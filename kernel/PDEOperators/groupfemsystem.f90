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
!#
!# Remark: The general internal data layout is as follows:
!# 
!# The scalar nodal vectors 1-6 contain the values 
!# for Pp, Pm, Qp, Qm, Rp, Rm in this order.
!#
!# The scalar edgewise vectors contain the following data
!#
!# 1 | alpha
!# 2 | total raw antidiffusive flux (explicit+implicit parts)
!# 3 | explicit part of the raw antidiffusive flux
!# 4 | constraining flux for prelimiting/implicit flux correction
!# </purpose>
!##############################################################################

module groupfemsystem

  use afcstabilisation
  use basicgeometry
  use collection
  use fsystem
  use genoutput
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

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of nodes to handle simultaneously when building matrices
#ifdef GFSYS_NEQSIM
  integer, parameter, public :: GFSYS_NEQSIM = GFSYS_NEQSIM
#else
  integer, parameter, public :: GFSYS_NEQSIM = 128
#endif

  ! Number of edges to handle simultaneously when building matrices
#ifdef GFSYS_NEDGESIM
  integer, parameter, public :: GFSYS_NEDGESIM = GFSYS_NEDGESIM
#else
  integer, parameter, public :: GFSYS_NEDGESIM = 64
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

    ! local variables
    integer, dimension(2) :: Isize


    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and. &
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
    if (rmatrixBlockTemplate%nblocksPerCol .ne. &
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
    rafcstab%NEDGE = int(0.5*(rmatrixBlockTemplate%RmatrixBlock(1,1)%NA-&
                              rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ))


    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationBlock', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)


    case (AFCSTAB_FEMTVD)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationBlock', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      allocate(rafcstab%p_rvectorPp)
      allocate(rafcstab%p_rvectorPm)
      allocate(rafcstab%p_rvectorQp)
      allocate(rafcstab%p_rvectorQm)
      allocate(rafcstab%p_rvectorRp)
      allocate(rafcstab%p_rvectorRm)

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      

    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationBlock', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      allocate(rafcstab%p_rvectorPp)
      allocate(rafcstab%p_rvectorPm)
      allocate(rafcstab%p_rvectorQp)
      allocate(rafcstab%p_rvectorQm)
      allocate(rafcstab%p_rvectorRp)
      allocate(rafcstab%p_rvectorRm)

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)

      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      allocate(rafcstab%p_rvectorAlpha)
      allocate(rafcstab%p_rvectorFlux0)
      allocate(rafcstab%p_rvectorFlux)

      call lsyssc_createVector(rafcstab%p_rvectorAlpha, rafcstab%NEDGE,&
          1, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorFlux0, rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorFlux,  rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ST_DOUBLE)

      ! We need the edgewise vector for the prelimited fluxes
      if (rafcstab%bprelimiting .or.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        allocate(rafcstab%p_rvectorFluxPrel)
        call lsyssc_createVector(rafcstab%p_rvectorFluxPrel,&
            rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      end if

      ! We need the nodal block vector for the low-order predictor
      allocate(rafcstab%p_rvectorPredictor)
      if (present(rblockDiscretisation)) then
        call lsysbl_createVectorBlock(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., ST_DOUBLE)
      else
        call lsysbl_createVectorBlock(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end if


    case (AFCSTAB_FEMFCT_LINEARISED,&
          AFCSTAB_FEMFCT_MASS)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationBlock', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      allocate(rafcstab%p_rvectorPp)
      allocate(rafcstab%p_rvectorPm)
      allocate(rafcstab%p_rvectorQp)
      allocate(rafcstab%p_rvectorQm)
      allocate(rafcstab%p_rvectorRp)
      allocate(rafcstab%p_rvectorRm)

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      allocate(rafcstab%p_rvectorAlpha)
      allocate(rafcstab%p_rvectorFlux)

      call lsyssc_createVector(rafcstab%p_rvectorAlpha, rafcstab%NEDGE,&
          1, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorFlux,  rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ST_DOUBLE)

      ! We need the edgewise vector if raw antidiffusive fluxes should be prelimited
      if (rafcstab%bprelimiting) then
        allocate(rafcstab%p_rvectorFluxPrel)
        call lsyssc_createVector(rafcstab%p_rvectorFluxPrel,&
            rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      end if

      ! We need the nodal block vector for the low-order predictor
      allocate(rafcstab%p_rvectorPredictor)
      if (present(rblockDiscretisation)) then
        call lsysbl_createVectorBlock(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., ST_DOUBLE)
      else
        call lsysbl_createVectorBlock(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
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
    integer, dimension(2) :: Isize


    ! Set atomic data
    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rmatrixTemplate%NVAR
    end if
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = int(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ))


    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationScalar', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)


    case (AFCSTAB_FEMTVD)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationScalar', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      allocate(rafcstab%p_rvectorPp)
      allocate(rafcstab%p_rvectorPm)
      allocate(rafcstab%p_rvectorQp)
      allocate(rafcstab%p_rvectorQm)
      allocate(rafcstab%p_rvectorRp)
      allocate(rafcstab%p_rvectorRm)

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)

      
    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationScalar', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      allocate(rafcstab%p_rvectorPp)
      allocate(rafcstab%p_rvectorPm)
      allocate(rafcstab%p_rvectorQp)
      allocate(rafcstab%p_rvectorQm)
      allocate(rafcstab%p_rvectorRp)
      allocate(rafcstab%p_rvectorRm)

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVAR, .false., ST_DOUBLE)

      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      allocate(rafcstab%p_rvectorAlpha)
      allocate(rafcstab%p_rvectorFlux0)
      allocate(rafcstab%p_rvectorFlux)

      call lsyssc_createVector(rafcstab%p_rvectorAlpha, rafcstab%NEDGE,&
          1, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorFlux0, rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorFlux,  rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ST_DOUBLE)

      ! We need the edgewise vector for the prelimited fluxes
      if (rafcstab%bprelimiting .or.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        allocate(rafcstab%p_rvectorFluxPrel)
        call lsyssc_createVector(rafcstab%p_rvectorFluxPrel,&
            rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
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
      

    case (AFCSTAB_FEMFCT_LINEARISED,&
          AFCSTAB_FEMFCT_MASS)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationScalar', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      allocate(rafcstab%p_rvectorPp)
      allocate(rafcstab%p_rvectorPm)
      allocate(rafcstab%p_rvectorQp)
      allocate(rafcstab%p_rvectorQm)
      allocate(rafcstab%p_rvectorRp)
      allocate(rafcstab%p_rvectorRm)

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ST_DOUBLE)
      
      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      allocate(rafcstab%p_rvectorAlpha)
      allocate(rafcstab%p_rvectorFlux)

      call lsyssc_createVector(rafcstab%p_rvectorAlpha, rafcstab%NEDGE,&
          1, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%p_rvectorFlux,  rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ST_DOUBLE)

      ! We need the edgewise vector if raw antidiffusive fluxes should be prelimited
      if (rafcstab%bprelimiting) then
        allocate(rafcstab%p_rvectorFluxPrel)
        call lsyssc_createVector(rafcstab%p_rvectorFluxPrel,&
            rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
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
      fcb_calcMatrixDiagonal_sim, fcb_calcMatrix_sim, dscale,&
      bclear, rdivMatrix, rcollection)

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
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The divergence operator
    type(t_matrixBlock), intent(inout) :: rdivMatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_array), dimension(rx%nblocks,rx%nblocks)  :: rarray
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal
    logical :: bisFullMatrix


    ! Check if block vector contains only one block and if
    ! global operator is stored in interleave format.
    if ((rx%nblocks .eq. 1) .and.&
        (rdivMatrix%nblocksPerCol .eq. 1) .and. &
        (rdivMatrix%nblocksPerRow .eq. 1)) then
      call gfsys_buildDivOperatorScalar(rafcstab, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagonal_sim, fcb_calcMatrix_sim,&
          dscale, bclear, rdivMatrix%RmatrixBlock(1,1), rcollection)
      return
    end if

    ! Check if block matrix exhibits group structure
    if (rdivMatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsysbl_clearMatrix(rdivMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DmatCoeffAtNode(rafcstab, p_DmatrixCoeffsAtNode)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call afcstab_getbase_array(rdivMatrix, rarray, bisFullMatrix)
    call lsysbl_getbase_double(rx, p_Dx)
    
    ! What kind of matrix are we?
    select case(rdivMatrix%RmatrixBlock(1,1)%cmatrixFormat)
    case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 7 and 9
      !-------------------------------------------------------------------------

      ! Set diagonal pointer
      if (rdivMatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX7) then
        call lsyssc_getbase_Kld(rdivMatrix%RmatrixBlock(1,1), p_Kdiagonal)
      else
        call lsyssc_getbase_Kdiagonal(rdivMatrix%RmatrixBlock(1,1), p_Kdiagonal)
      end if

      ! What type of matrix are we?
      if (bisFullMatrix) then
        
        call doOperatorMat79(p_Kdiagonal, p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rx%nblocks,&
            p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, rarray)
     
      else   ! bisFullMatrix == no

        call doOperatorMat79Diag(p_Kdiagonal, p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rx%nblocks,&
            p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, rarray)

      end if   ! bisFullMatrix

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble block-diagonal divergence operator K
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag(Kdiagonal, IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, K)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE, NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: K

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IverticesAtNode
      
      ! local variables
      integer :: idx,IEQset,IEQmax,IEDGEset,IEDGEmax
      integer :: i,ii,jj,ij,ji,iedge,ivar

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IverticesAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(NVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
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
        call fcb_calcMatrixDiagonal_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IverticesAtNode(:,1:IEQmax-IEQset+1), dscale,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)
        
        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        do idx = 1, IEQmax-IEQset+1

          ! Get position of diagonal entry
          ii = IverticesAtNode(2,idx)

          ! Update the diagonal coefficient
          do ivar = 1, NVAR
            K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)+&
                DcoefficientsAtNode(ivar,1,idx)
          end do
        end do       
      end do
      
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

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute off-diagonal entries
        call fcb_calcMatrix_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global matrix
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of diagonal entries
          ii = Kdiagonal(IverticesAtEdge(1,iedge))
          jj = Kdiagonal(IverticesAtEdge(2,iedge))

          ! Get position of off-diagonal entries
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

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
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)

    end subroutine doOperatorMat79Diag

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, K)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE, NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: K

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IverticesAtNode
      
      ! local variables
      integer :: idx,IEQset,IEQmax,IEDGEset,IEDGEmax
      integer :: i,ii,jj,ij,ji,iedge,ivar,jvar,ijpos

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IverticesAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(NVAR*NVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
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
        call fcb_calcMatrixDiagonal_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IverticesAtNode(:,1:IEQmax-IEQset+1), dscale,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
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
      end do

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

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute off-diagonal entries
        call fcb_calcMatrix_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of diagonal entries
          ii = Kdiagonal(IverticesAtEdge(1,iedge))
          jj = Kdiagonal(IverticesAtEdge(2,iedge))

          ! Get position of off-diagonal entries
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

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
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)

    end subroutine doOperatorMat79

  end subroutine gfsys_buildDivOperatorBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorScalar(rafcstab, rx,&
      fcb_calcMatrixDiagonal_sim, fcb_calcMatrix_sim, dscale,&
      bclear, rdivMatrix, rcollection)

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
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The divergence operator
    type(t_matrixScalar), intent(inout) :: rdivMatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DivOp,p_Dx
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rdivMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DmatCoeffAtNode(rafcstab, p_DmatrixCoeffsAtNode)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsyssc_getbase_double(rdivMatrix, p_DivOp)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! What kind of matrix are we?
    select case(rdivMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7INTL, LSYSSC_MATRIX9INTL)
      !-------------------------------------------------------------------------
      ! Matrix format 7 and 9 interleaved
      !-------------------------------------------------------------------------

      ! Set diagonal pointer
      if (rdivMatrix%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
        call lsyssc_getbase_Kld(rdivMatrix, p_Kdiagonal)
      else
        call lsyssc_getbase_Kdiagonal(rdivMatrix, p_Kdiagonal)
      end if

      ! What type of matrix are we?
      select case(rdivMatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIX1)
        call doOperatorMat79(p_Kdiagonal, p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rdivMatrix%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
            p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_DivOp)
        
      case (LSYSSC_MATRIXD)
        call doOperatorMat79(p_Kdiagonal, p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rdivMatrix%NA, rx%NVAR, rx%NVAR,&
            p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_DivOp)

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

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence operator K
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79(Kdiagonal, IverticesAtEdge, NEDGE, NEQ, NA,&
        NVAR, MVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, K)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: K

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IverticesAtNode
      
      ! local variables
      integer :: idx,IEQset,IEQmax,IEDGEset,IEDGEmax
      integer :: i,ii,jj,ij,ji,iedge

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IverticesAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(MVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
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
        call fcb_calcMatrixDiagonal_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IverticesAtNode(:,1:IEQmax-IEQset+1), dscale,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        do idx = 1, IEQmax-IEQset+1

          ! Get position of diagonal entry
          ii = IverticesAtNode(2,idx)

          ! Update the diagonal coefficient
          K(:,ii) = K(:,ii) + DcoefficientsAtNode(:,1,idx)
        end do
      end do

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

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
        end do

        ! Use callback function to compute off-diagonal entries
        call fcb_calcMatrix_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
        
        ! Loop through all edges in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of diagonal entries
          ii = Kdiagonal(IverticesAtEdge(1,iedge))
          jj = Kdiagonal(IverticesAtEdge(2,iedge))

          ! Get position of off-diagonal entries
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Update the global operator
          K(:,ii) = K(:,ii) - DcoefficientsAtEdge(:,1,idx)
          K(:,jj) = K(:,jj) - DcoefficientsAtEdge(:,1,idx)
          K(:,ij) = K(:,ij) + DcoefficientsAtEdge(:,2,idx) + DcoefficientsAtEdge(:,1,idx) 
          K(:,ji) = K(:,ji) + DcoefficientsAtEdge(:,3,idx) + DcoefficientsAtEdge(:,1,idx) 
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      
    end subroutine doOperatorMat79

  end subroutine gfsys_buildDivOperatorScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVectorBlock(rafcstab, rx,&
      fcb_calcFlux_sim, dscale, bclear, ry, rcollection)

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

    ! callback functions to compute local fluxes
    include 'intf_gfsyscallback.inc'
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
    integer, dimension(:,:), pointer :: p_IverticesAtEdge


    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildDivVectorScalar(rafcstab, rx%RvectorBlock(1),&
          fcb_calcFlux_sim, dscale, bclear, ry%RvectorBlock(1), rcollection)
      return
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorBlock')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)
    
    ! Assemble the divergence vector
    call doDivVector(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
        rx%nblocks, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dy)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector
    
    subroutine doDivVector(IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtEdge, Dx, dscale, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
          ! Update the global vector
          Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
          Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      
    end subroutine doDivVector
    
  end subroutine gfsys_buildDivVectorBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVectorScalar(rafcstab, rx,&
      fcb_calcFlux_sim, dscale, bclear, ry, rcollection)

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
    include 'intf_gfsyscallback.inc'
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
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorScalar')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

    ! Assemble the divergence vector
    call doDivVector(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
        rx%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dy)
    
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector

    subroutine doDivVector(IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtEdge, Dx, dscale, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
        end do

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Update the global vector
          Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
          Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      
    end subroutine doDivVector

  end subroutine gfsys_buildDivVectorScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecTVDBlock(RcoeffMatrices, rafcstab, rx,&
      fcb_calcFlux_sim, fcb_calcCharacteristics_sim,&
      dscale, bclear, ry, rcollection)

!<description>
    ! This subroutine assembles the divergence vector for FEM-TVD schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
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
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim


    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildDivVecTVDScalar(RcoeffMatrices, rafcstab,&
          rx%RvectorBlock(1), fcb_calcFlux_sim,&
          fcb_calcCharacteristics_sim, dscale, bclear,&
          ry%RvectorBlock(1), rcollection)
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)


    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDBlock')
      call sys_halt()
    end select

    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 7 and 9
      !-------------------------------------------------------------------------

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doLimitTVDMat79_1D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
            p_DcoeffX, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
      case (NDIM2D)
        call doLimitTVDMat79_2D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
            p_DcoeffX, p_DcoeffY, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
      case (NDIM3D)
        call doLimitTVDMat79_3D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
      end select

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDBlock')
      call sys_halt()
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax
      

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DmatrixCoeffsAtEdge(1,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))

      ! Clear P's and Q's (X-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Assemble high-order Galerkin fluxes
        !-----------------------------------------------------------------------

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(1, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)
      
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (X-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the X-direction
        call fcb_calcCharacteristics_sim(XDir1D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DmatrixCoeffsAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)

    end subroutine doLimitTVDMat79_1D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax

      
      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DmatrixCoeffsAtEdge(2,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))

      ! Clear P's and Q's (X-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Assemble high-order Galerkin fluxes
        !-----------------------------------------------------------------------

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1
          
          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(1, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do
      
      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------
      
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (X-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the X-direction
        call fcb_calcCharacteristics_sim(XDir2D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)

        !-----------------------------------------------------------------------
        ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
        !-----------------------------------------------------------------------
        
        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir2D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(2, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (Y-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir2D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)
      end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DmatrixCoeffsAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)

    end subroutine doLimitTVDMat79_2D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      !input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DmatrixCoeffsAtEdge(3,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (X-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Assemble high-order Galerkin fluxes
        !-----------------------------------------------------------------------

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1
          
          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(1, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do
      
      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------
      
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (X-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the X-direction
        call fcb_calcCharacteristics_sim(XDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)

        !-----------------------------------------------------------------------
        ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
        !-----------------------------------------------------------------------
        
        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(2, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (Y-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)

        !-----------------------------------------------------------------------
        ! Assemble artificial viscosities and antidiffusive fluxes (Z-direction)
        !-----------------------------------------------------------------------
        
        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Z-direction
        call fcb_calcCharacteristics_sim(ZDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)
        
        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(3, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Z-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx)         = Dx(IverticesAtEdge(2,iedge),:)
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (Z-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Z-direction
        call fcb_calcCharacteristics_sim(ZDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(3, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DmatrixCoeffsAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)

    end subroutine doLimitTVDMat79_3D


    !**************************************************************
    ! Assemble the upper and lower bounds Q and the sums of
    ! antidiffusive contributions P for a given set of edges

    pure subroutine doBoundsAndIncrements_sim(idirection, NVAR,&
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IverticesAtEdge, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dpp,Dpm,Dqp,Dqm
      
      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      integer :: idx,ivar,i,j
      
      ! Loop over all edges in the set
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
            i = IverticesAtEdge(2,idx)
            j = IverticesAtEdge(1,idx)
            Dflux(ivar) = -Dflux(ivar)
          else
            i = IverticesAtEdge(1,idx)
            j = IverticesAtEdge(2,idx)
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

    end subroutine doBoundsAndIncrements_sim

    
    !**************************************************************
    ! Limit the antidiffusive fluxes and apply them to the vector
    
    pure subroutine doLimitADFluxes_sim(idirection, NVAR, dscale, Drp, Drm,&
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IverticesAtEdge, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      real(DP), dimension(:,:), intent(in) :: DrighteigenvectorsAtEdge,Drp,Drm
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dy
      
      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      real(DP) :: daux
      integer :: idx,ivar,jvar,i,j
      
      ! Loop over all edges in the set
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
        i = IverticesAtEdge(1,idx)
        j = IverticesAtEdge(2,idx)

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

    end subroutine doLimitADFluxes_sim

  end subroutine gfsys_buildDivVecTVDBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecTVDScalar(RcoeffMatrices, rafcstab, rx,&
      fcb_calcFlux_sim, fcb_calcCharacteristics_sim,&
      dscale, bclear, ry, rcollection)

!<description>
    ! This subroutine assembles the divergence vector for FEM-TVD schemes
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
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
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim


    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDScalar')
      call sys_halt()
    end select

    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 7 and 9
      !-------------------------------------------------------------------------

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doLimitTVDMat79_1D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%NVAR,&
            p_DcoeffX, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
      case (NDIM2D)
        call doLimitTVDMat79_2D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%NVAR,&
            p_DcoeffX, p_DcoeffY, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
      case (NDIM3D)
        call doLimitTVDMat79_3D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%NVAR,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
      end select

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDScalar')
      call sys_halt()
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameter
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax
      

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DmatrixCoeffsAtEdge(1,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (X-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Assemble high-order Galerkin fluxes
        !-----------------------------------------------------------------------

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(1, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do
      
      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)
      
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (X-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the X-direction
        call fcb_calcCharacteristics_sim(XDir1D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DmatrixCoeffsAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)

    end subroutine doLimitTVDMat79_1D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax

      
      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DmatrixCoeffsAtEdge(2,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))

      ! Clear P's and Q's (X-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
        end do
        
        !-----------------------------------------------------------------------
        ! Assemble high-order Galerkin fluxes
        !-----------------------------------------------------------------------

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1
          
          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(1, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do
      
      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------
      
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (X-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the X-direction
        call fcb_calcCharacteristics_sim(XDir2D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)

        !-----------------------------------------------------------------------
        ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
        !-----------------------------------------------------------------------
        
        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir2D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(2, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (Y-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir2D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)
      end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DmatrixCoeffsAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)

    end subroutine doLimitTVDMat79_2D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax

      
      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DmatrixCoeffsAtEdge(3,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (X-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Assemble high-order Galerkin fluxes
        !-----------------------------------------------------------------------

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1
          
          ! Get actual edge number
          iedge = idx+IEDGEset-1
          
          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(1, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do
      
      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------
      
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do
      
        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (X-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the X-direction
        call fcb_calcCharacteristics_sim(XDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)

        !-----------------------------------------------------------------------
        ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
        !-----------------------------------------------------------------------
        
        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(2, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (Y-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Y-direction
        call fcb_calcCharacteristics_sim(YDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)

        !-----------------------------------------------------------------------
        ! Assemble artificial viscosities and antidiffusive fluxes (Z-direction)
        !-----------------------------------------------------------------------
        
        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Z-direction
        call fcb_calcCharacteristics_sim(ZDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)
        
        ! Assemble the upper and lower bounds Q and the sums of
        ! antidiffusive contributions P for the set of edges
        call doBoundsAndIncrements_sim(3, NVAR,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax),&
            Dpp, Dpm, Dqp, Dqm)
      end do

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Z-direction)
      !-------------------------------------------------------------------------

      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx)         = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx)         = Dx(:,IverticesAtEdge(2,iedge))
          DmatrixCoeffsAtEdge(1,1,idx) = DcoeffX(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(2,1,idx) = DcoeffY(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(3,1,idx) = DcoeffZ(IverticesAtEdge(3,iedge))
          DmatrixCoeffsAtEdge(1,2,idx) = DcoeffX(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(2,2,idx) = DcoeffY(IverticesAtEdge(4,iedge))
          DmatrixCoeffsAtEdge(3,2,idx) = DcoeffZ(IverticesAtEdge(4,iedge))
        end do

        !-----------------------------------------------------------------------
        ! Apply artificial viscosities and limited antidiffusion (Z-direction)
        !-----------------------------------------------------------------------

        ! Use callback function to compute the characteristic variables
        ! and corresponding eigenvalues along the Z-direction
        call fcb_calcCharacteristics_sim(ZDir3D,&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection=rcollection)

        ! Apply limited characteristic fluxes to global vector
        call doLimitADFluxes_sim(3, NVAR, dscale, Drp, Drm,&
            DmatrixCoeffsAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), Dy)
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DmatrixCoeffsAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      
    end subroutine doLimitTVDMat79_3D

    !**************************************************************
    ! Assemble the upper and lower bounds Q and the sums of
    ! antidiffusive contributions P for a given set of edges

    pure subroutine doBoundsAndIncrements_sim(idirection, NVAR,&
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IverticesAtEdge, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dpp,Dpm,Dqp,Dqm
      
      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      integer :: idx,ivar,i,j
      
      ! Loop over all edges in the set
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
            i = IverticesAtEdge(2,idx)
            j = IverticesAtEdge(1,idx)
            Dflux(ivar) = -Dflux(ivar)
          else
            i = IverticesAtEdge(1,idx)
            j = IverticesAtEdge(2,idx)
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

    end subroutine doBoundsAndIncrements_sim

    !**************************************************************
    ! Limit the antidiffusive fluxes and apply them to the vector
    
    pure subroutine doLimitADFluxes_sim(idirection, NVAR, dscale, Drp, Drm,&
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IverticesAtEdge, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      real(DP), dimension(:,:), intent(in) :: DrighteigenvectorsAtEdge,Drp,Drm
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dy
      
      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      real(DP) :: daux
      integer :: idx,ivar,jvar,i,j
      
      ! Loop over all edges in the set
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
        i = IverticesAtEdge(1,idx)
        j = IverticesAtEdge(2,idx)

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

    end subroutine doLimitADFluxes_sim

  end subroutine gfsys_buildDivVecTVDScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecFCTBlock(rlumpedMassMatrix,&
      rafcstab, rx, dscale, bclear, ioperationSpec, ry, NVARtransformed, &
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
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

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
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation_sim
    optional :: fcb_calcDiffTransformation_sim

    ! OPTIONAL: callback functions to overwrite the standard operations
    include 'intf_groupfemcallback.inc'
    optional :: fcb_calcADIncrements
    optional :: fcb_calcBounds
    optional :: fcb_limitNodal
    optional :: fcb_limitEdgewise
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
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: nvariable

    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1)) then
      call gfsys_buildDivVecFCTScalar(rlumpedMassMatrix,&
          rafcstab, rx%RvectorBlock(1), dscale, bclear,&
          ioperationSpec, ry%RvectorBlock(1), NVARtransformed,&
          fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
          fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
          fcb_limitEdgewise, fcb_calcCorrection, rcollection)
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBLock')
      call sys_halt()
    end if

    !---------------------------------------------------------------------------
    ! The nonlinear FEM-FCT algorithm is split into the following
    ! steps which can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Compute the antidiffusive increments (Pp, Pm)
    !
    ! 3) Compute the local solution bounds (Qp, Qm).
    !
    ! 3) Compute the nodal correction factors (Rp, Rm).
    !
    ! 4) Apply the limited antidifusive fluxes to the divergence
    !
    !    Step 4) may be split into the following substeps
    !
    !    4.1) Compute the edgewise correction factors based on the pre-
    !         computed raw-antidiffusive fluxes.
    !
    !    4.2) Compute the raw antidiffusive fluxes for a different set of
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

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if

    if (iand(ioperationSpec, AFCSTAB_FCTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute sums of antidiffusive increments
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
      call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux)
      else
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      end if

      ! Compute sums of antidiffusive increments
      if (rafcstab%bprelimiting) then
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux0)

        if (present(fcb_calcADIncrements)) then
          ! User-supplied callback routine
          call fcb_calcADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dflux, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim, p_Dflux0,&
              rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doPreADIncrementsTransformed(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
        else
          ! Standard routine without flux transformation
          call doPreADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
        end if
        
      else

        if (present(fcb_calcADIncrements)) then
          ! User-supplied callback routine
          call fcb_calcADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dflux, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim,&
              rcollection=rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doADIncrementsTransformed(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              nvariable, p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        else
          ! Standard routine without flux transformation
          call doADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        end if
      end if
        
      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute local bounds
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
      call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

       ! Compute local bounds
      if (present(fcb_calcBounds)) then
        ! User-supplied callback routine
        call fcb_calcBounds(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dqp, p_Dqm,&
            fcb_calcDiffTransformation_sim, rcollection)
      elseif (present(fcb_calcDiffTransformation_sim)) then
        ! Standard routine with difference transformation
        call doBoundsTransformed(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            nvariable, p_Dx, p_Dqp, p_Dqm)
        ! Standard routine without difference transformation
      else
        call doBounds(p_IverticesAtEdge, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dqp, p_Dqm)
      end if
      
      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_BOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_BOUNDS)       .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
      call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
      call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
      call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
      call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
      call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
      call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, nvariable, dscale,&
            p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection)
      elseif (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        ! Standard routine without constraints
        call doLimitNodal(rafcstab%NEQ, nvariable,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        ! Standard routine with constraints
        call doLimitNodalConstrained(rafcstab%NEQ, nvariable,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
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
        call output_line('Stabilisation does not provides nodal correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
      call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux0)

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dflux, p_Drp, p_Drm,&
              p_Dalpha, fcb_calcFluxTransformation_sim, p_Dflux0, rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseConstrainedTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, nvariable, p_Dx,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseConstrained(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      else

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dflux, p_Drp, p_Drm,&
              p_Dalpha, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, nvariable, p_Dx,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewise(p_IverticesAtEdge,&
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
        call output_line('Stabilisation does not provides edgewise correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(ry, p_Dy)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Clear divergence vector?
      if (bclear) call lsysbl_clearVector(ry)

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then
        call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
        
        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IverticesAtEdge, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, dscale, p_Dalpha, p_Dflux,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dy, p_ML,rcollection)
        else
          ! Standard routine
          call doCorrectScaleByMass(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
        end if
        
      else

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IverticesAtEdge, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, dscale, p_Dalpha, p_Dflux,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dy, rcollection=rcollection)
        else
          ! Standard routine
          call doCorrect(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              dscale, p_Dalpha, p_Dflux, p_Dy)
        end if
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j


      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Apply multiplicative correction factor
        F_ij = Dalpha(iedge) * Dflux(:,iedge)

        ! Compute the sums of antidiffusive increments
        Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
        Dpp(:,j) = Dpp(:,j)+max(0.0_DP,-F_ij)
        Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
        Dpm(:,j) = Dpm(:,j)+min(0.0_DP,-F_ij)
      end do

    end subroutine doADIncrements

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-
    ! defined set of variables prior to computing the sums

    subroutine doADIncrementsTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
      
      ! local variables
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
          DfluxesAtEdge(:,idx) = Dalpha(iedge)*Dflux(:,iedge)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1), &
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vectors
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Compute the sums of positive/negative antidiffusive increments
          Dpp(:,i) = Dpp(:,i)+max(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
          Dpp(:,j) = Dpp(:,j)+max(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
          Dpm(:,i) = Dpm(:,i)+min(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
          Dpm(:,j) = Dpm(:,j)+min(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      deallocate(DtransformedFluxesAtEdge)

    end subroutine doADIncrementsTransformed

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and with prelimiting

    subroutine doPreADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Dflux0, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      real(DP) :: alpha_ij
      integer :: iedge,i,j


      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Apply multiplicative correction factor
        F_ij = Dalpha(iedge) * Dflux(:,iedge)

        ! MinMod prelimiting
        alpha_ij = minval(mprim_minmod3(F_ij, Dflux0(:,iedge), F_ij))

        ! Synchronisation of correction factors
        Dalpha(iedge) = Dalpha(iedge) * alpha_ij

        ! Update the raw antidiffusive flux
        F_ij = alpha_ij * F_ij

        ! Compute the sums of antidiffusive increments
        Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
        Dpp(:,j) = Dpp(:,j)+max(0.0_DP,-F_ij)
        Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
        Dpm(:,j) = Dpm(:,j)+min(0.0_DP,-F_ij)
      end do

    end subroutine doPreADIncrements

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-
    ! defined set of variables prior to computing the sums
    ! Perform minmod prelimiting of the raw antidiffusive fluxes

    subroutine doPreADIncrementsTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Dflux0,&
        Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(NVARtransformed,NEQ), intent(inout) :: Dpp,Dpm

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedPrelFluxesAtEdge
      
      ! local variables
      real(DP), dimension(NVAR) :: F_ij,F_ji
      real(DP) :: alpha_ij,alpha_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      allocate(DtransformedPrelFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
          DfluxesAtEdge(:,idx) = Dalpha(iedge)*Dflux(:,iedge)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1), &
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Use callback function to compute transformed fluxes
        ! for the explicit part for prelimiting
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux0(:,IEDGEset:IEDGEmax),&
            DtransformedPrelFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! MinMod prelimiting
          alpha_ij = minval(mprim_minmod3(&
                            DtransformedFluxesAtEdge(:,1,idx),&
                            DtransformedPrelFluxesAtEdge(:,1,idx), F_ij))
          alpha_ji = minval(mprim_minmod3(&
                            DtransformedFluxesAtEdge(:,2,idx),&
                            DtransformedPrelFluxesAtEdge(:,2,idx), F_ji))
          
          ! Synchronisation of correction factors TODO!!!!
          Dalpha(iedge) = Dalpha(iedge) * alpha_ij
          
          ! Update the raw antidiffusive fluxes
          F_ij = alpha_ij * F_ij
          F_ji = alpha_ij * F_ji

          ! Compute the sums of positive/negative antidiffusive increments
          Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
          Dpp(:,j) = Dpp(:,j)+max(0.0_DP, F_ji)
          Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
          Dpm(:,j) = Dpm(:,j)+min(0.0_DP, F_ji)
        end do
      end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      deallocate(DtransformedPrelFluxesAtEdge)

    end subroutine doPreADIncrementsTransformed

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! without transformation

    subroutine doBounds(IverticesAtEdge, NEDGE, NEQ, NVAR, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Diff
      integer :: iedge,i,j


      ! Clear Q`s
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute solution difference
        Diff = Dx(j,:)-Dx(i,:)

        ! Compute the distance to a local extremum
        ! of the predicted solution
        Dqp(:,i) = max(Dqp(:,i), Diff)
        Dqp(:,j) = max(Dqp(:,j),-Diff)
        Dqm(:,i) = min(Dqm(:,i), Diff)
        Dqm(:,j) = min(Dqm(:,j),-Diff)
      end do

    end subroutine doBounds

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! which is transformed to a user-defined set of variables

    subroutine doBoundsTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(inout) :: Dqp,Dqm

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DtransformedDataAtEdge

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
      

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedDataAtEdge(NVARtransformed,GFSYS_NEDGESIM))
      
      ! Clear Q`s
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute transformed differences
        call fcb_calcDiffTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DtransformedDataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
        
          ! Compute the distance to a local extremum of the predicted solution
          Dqp(:,i) = max(Dqp(:,i), DtransformedDataAtEdge(:,idx))
          Dqp(:,j) = max(Dqp(:,j),-DtransformedDataAtEdge(:,idx))
          Dqm(:,i) = min(Dqm(:,i), DtransformedDataAtEdge(:,idx))
          Dqm(:,j) = min(Dqm(:,j),-DtransformedDataAtEdge(:,idx))
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedDataAtEdge)

    end subroutine doBoundsTransformed

    !**************************************************************
    ! Compute nodal correction factors without constraints

    subroutine doLimitNodal(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq


      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpp(:,ieq) .gt. AFCSTAB_EPSABS)
          Drp(:,ieq) = ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq))
        elsewhere
          Drp(:,ieq) = 1.0_DP
        end where
      end do

      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpm(:,ieq) .lt. -AFCSTAB_EPSABS)
          Drm(:,ieq) = ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq))
        elsewhere
          Drm(:,ieq) = 1.0_DP
        end where
      end do

    end subroutine doLimitNodal

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrained(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq


      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpp(:,ieq) .gt. AFCSTAB_EPSABS)
          Drp(:,ieq) = min(1.0_DP, ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq)))
        elsewhere
          Drp(:,ieq) = 1.0_DP
        end where
      end do

      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpm(:,ieq) .lt. -AFCSTAB_EPSABS)
          Drm(:,ieq) = min(1.0_DP, ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq)))
        elsewhere
          Drm(:,ieq) = 1.0_DP
        end where
      end do

    end subroutine doLimitNodalConstrained

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewise(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F_ij,R_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Get precomputed raw antidiffusive fluxes
        F_ij = Dflux(:,iedge)

        ! Compute nodal correction factors
        where (F_ij .ge. 0.0_DP)
          R_ij = min(Drp(:,i),Drm(:,j))
        elsewhere
          R_ij = min(Drp(:,j),Drm(:,i))
        end where

!!$       REMARK: Numerical test demonstrate that this modification has no
!!$               significant influence on the order-of-accuracy
!!$
!!$        where (F_ij .gt. AFCSTAB_EPSABS)
!!$          R_ij = min(Drp(:,i),Drm(:,j))
!!$        elsewhere (F_ij .lt. -AFCSTAB_EPSABS)
!!$          R_ij = min(Drp(:,j),Drm(:,i))
!!$        elsewhere
!!$          R_ij = 1.0_DP
!!$        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(R_ij)
      end do

    end subroutine doLimitEdgewise

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    ! which are transformed to a user-defined set of variables
    ! priori to computing the correction factors

    subroutine doLimitEdgewiseTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Drp, Drm, Dalpha)

      ! input  parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux(:,IEDGEset:IEDGEmax),&
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Compute nodal correction factors
          R_ij = merge(Drp(:,i), Drm(:,i),&
                       DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
          R_ji = merge(Drp(:,j), Drm(:,j),&
                       DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)

!!$       REMARK: Numerical test demonstrate that this modification has no
!!$               significant influence on the order-of-accuracy
!!$
!!$          where (DtransformedFluxesAtEdge(:,1,idx) .gt. AFCSTAB_EPSABS)
!!$            R_ij = Drp(:,i)
!!$          elsewhere (DtransformedFluxesAtEdge(:,1,idx) .lt. -AFCSTAB_EPSABS)
!!$            R_ij = Drm(:,i)
!!$          elsewhere
!!$            R_ij = 1.0_DP
!!$          end where
!!$
!!$          where (DtransformedFluxesAtEdge(:,2,idx) .gt. AFCSTAB_EPSABS)
!!$            R_ji = Drp(:,j)
!!$          elsewhere (DtransformedFluxesAtEdge(:,2,idx) .lt. -AFCSTAB_EPSABS)
!!$            R_ji = Drm(:,j)
!!$          elsewhere
!!$            R_ij = 1.0_DP
!!$          end where

          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxesAtEdge)
     
    end subroutine doLimitEdgewiseTransformed

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrained(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F1_ij,F2_ij,R_ij
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

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

    end subroutine doLimitEdgewiseConstrained

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes which are transformed
    ! to a user-defined set of variables priori to computing the
    ! correction factors

    subroutine doLimitEdgewiseConstrainedTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux1, Dflux2, Drp, Drm, Dalpha)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes1AtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes2AtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes1AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes2AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux1(:,IEDGEset:IEDGEmax),&
            DtransformedFluxes1AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux2(:,IEDGEset:IEDGEmax),&
            DtransformedFluxes2AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1
          
          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

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
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxes1AtEdge)
      deallocate(DtransformedFluxes2AtEdge)
      
    end subroutine doLimitEdgewiseConstrainedTransformed

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrect(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, dscale, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Correct antidiffusive flux
        F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

        ! Apply limited antidiffusive fluxes
        Dy(i,:) = Dy(i,:) + F_ij
        Dy(j,:) = Dy(j,:) - F_ij
      end do

    end subroutine doCorrect

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

    subroutine doCorrectScaleByMass(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, dscale, ML, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha,ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Correct antidiffusive flux
        F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

        ! Apply limited antidiffusive fluxes
        Dy(i,:) = Dy(i,:) + F_ij/ML(i)
        Dy(j,:) = Dy(j,:) - F_ij/ML(j)
      end do
    end subroutine doCorrectScaleByMass

  end subroutine gfsys_buildDivVecFCTBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecFCTScalar(rlumpedMassMatrix,&
      rafcstab, rx, dscale, bclear, ioperationSpec, ry, NVARtransformed,&
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
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

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
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation_sim
    optional :: fcb_calcDiffTransformation_sim

    ! OPTIONAL: callback functions to overwrite the standard operations
    include 'intf_groupfemcallback.inc'
    optional :: fcb_calcADIncrements
    optional :: fcb_calcBounds
    optional :: fcb_limitNodal
    optional :: fcb_limitEdgewise
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
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: nvariable

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
      call sys_halt()
    end if

    !---------------------------------------------------------------------------
    ! The nonlinear FEM-FCT algorithm is split into the following
    ! steps which can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Compute the antidiffusive increments (Pp, Pm)
    !
    ! 3) Compute the local solution bounds (Qp, Qm).
    !
    ! 3) Compute the nodal correction factors (Rp, Rm).
    !
    ! 4) Apply the limited antidifusive fluxes to the divergence
    !
    !    Step 4) may be split into the following substeps
    !
    !    4.1) Compute the edgewise correction factors based on the pre-
    !         computed raw-antidiffusive fluxes.
    !
    !    4.2) Compute the raw antidiffusive fluxes for a different set of
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

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute sums of antidiffusive increments
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
      call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux)
      else
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      end if

      ! Compute sums of antidiffusive increments
      if (rafcstab%bprelimiting) then
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux0)

        if (present(fcb_calcADIncrements)) then
          ! User-supplied callback routine
          call fcb_calcADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim, p_Dflux0,&
              rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doPreADIncrementsTransformed(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
        else
          ! Standard routine without flux transformation
          call doPreADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
        end if
        
      else
        
        if (present(fcb_calcADIncrements)) then
          ! User-supplied callback routine
          call fcb_calcADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim,&
              rcollection=rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doADIncrementsTransformed(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              nvariable, p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        else
          ! Standard routine without flux transformation
          call doADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        end if
      end if

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute local bounds
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
      call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute local bounds
      if (present(fcb_calcBounds)) then
        ! User-supplied callback routine
        call fcb_calcBounds(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
            rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dqp, p_Dqm,&
            fcb_calcDiffTransformation_sim, rcollection)
      elseif (present(fcb_calcDiffTransformation_sim)) then
        ! Standard routine with difference transformation
        call doBoundsTransformed(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            nvariable, p_Dx, p_Dqp, p_Dqm)
        ! Standard routine without difference transformation
      else
        call doBounds(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            p_Dx, p_Dqp, p_Dqm)
      end if

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_BOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_BOUNDS)       .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
      call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
      call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
      call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
      call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
      call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
      call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, nvariable, dscale,&
            p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection)
      elseif (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        ! Standard routine without constraints
        call doLimitNodal(rafcstab%NEQ, nvariable,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        ! Standard routine with constraints
        call doLimitNodalConstrained(rafcstab%NEQ, nvariable,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
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
        call output_line('Stabilisation does not provides nodal correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
      call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux0)

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Drp, p_Drm,&
              p_Dalpha, fcb_calcFluxTransformation_sim, p_Dflux0, rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseConstrainedTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, nvariable, p_Dx,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseConstrained(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      else

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Drp, p_Drm,&
              p_Dalpha, fcb_calcFluxTransformation_sim,&
              rcollection=rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, nvariable, p_Dx,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewise(p_IverticesAtEdge,&
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
        call output_line('Stabilisation does not provides edgewise correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(ry, p_Dy)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Clear divergence vector?
      if (bclear) call lsyssc_clearVector(ry)

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then
        call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
        
        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IverticesAtEdge, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, dscale, p_Dalpha, p_Dflux,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dy, p_ML,rcollection)
        else
          ! Standard routine
          call doCorrectScaleByMass(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
        end if
        
      else

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IverticesAtEdge, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, dscale, p_Dalpha, p_Dflux,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dy, rcollection=rcollection)
        else
          ! Standard routine
          call doCorrect(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              dscale, p_Dalpha, p_Dflux, p_Dy)
        end if
      end if
    end if
    
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j


      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Apply multiplicative correction factor
        F_ij = Dalpha(iedge) * Dflux(:,iedge)

        ! Compute the sums of antidiffusive increments
        Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
        Dpp(:,j) = Dpp(:,j)+max(0.0_DP,-F_ij)
        Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
        Dpm(:,j) = Dpm(:,j)+min(0.0_DP,-F_ij)
      end do

    end subroutine doADIncrements

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-
    ! defined set of variables prior to computing the sums

    subroutine doADIncrementsTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(inout) :: Dpp,Dpm

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
      
      ! local variables
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
          DfluxesAtEdge(:,idx) = Dalpha(iedge)*Dflux(:,iedge)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1), &
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vectors
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Compute the sums of positive/negative antidiffusive increments
          Dpp(:,i) = Dpp(:,i)+max(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
          Dpp(:,j) = Dpp(:,j)+max(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
          Dpm(:,i) = Dpm(:,i)+min(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
          Dpm(:,j) = Dpm(:,j)+min(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      deallocate(DtransformedFluxesAtEdge)

    end subroutine doADIncrementsTransformed

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and with prelimiting

    subroutine doPreADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Dflux0, Dalpha, Dpp, Dpm)
      
      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      real(DP) :: alpha_ij
      integer :: iedge,i,j


      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Apply multiplicative correction factor
        F_ij = Dalpha(iedge) * Dflux(:,iedge)

        ! MinMod prelimiting
        alpha_ij = minval(mprim_minmod3(F_ij, Dflux0(:,iedge), F_ij))

        ! Synchronisation of correction factors
        Dalpha(iedge) = Dalpha(iedge) * alpha_ij

        ! Update the raw antidiffusive flux
        F_ij = alpha_ij * F_ij

        ! Compute the sums of antidiffusive increments
        Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
        Dpp(:,j) = Dpp(:,j)+max(0.0_DP,-F_ij)
        Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
        Dpm(:,j) = Dpm(:,j)+min(0.0_DP,-F_ij)
      end do

    end subroutine doPreADIncrements

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-
    ! defined set of variables prior to computing the sums
    ! Perform minmod prelimiting of the raw antidiffusive fluxes

    subroutine doPreADIncrementsTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Dflux0,&
        Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/outpu parameters
      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(NVARtransformed,NEQ), intent(inout) :: Dpp,Dpm

      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedPrelFluxesAtEdge
      
      ! local variables
      real(DP), dimension(NVAR) :: F_ij,F_ji
      real(DP) :: alpha_ij,alpha_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      allocate(DtransformedPrelFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Clear P`s
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
          DfluxesAtEdge(:,idx) = Dalpha(iedge)*Dflux(:,iedge)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1), &
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Use callback function to compute transformed fluxes
        ! for the explicit part for prelimiting
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux0(:,IEDGEset:IEDGEmax),&
            DtransformedPrelFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! MinMod prelimiting
          alpha_ij = minval(mprim_minmod3(&
                            DtransformedFluxesAtEdge(:,1,idx),&
                            DtransformedPrelFluxesAtEdge(:,1,idx), F_ij))
          alpha_ji = minval(mprim_minmod3(&
                            DtransformedFluxesAtEdge(:,2,idx),&
                            DtransformedPrelFluxesAtEdge(:,2,idx), F_ji))
          
          ! Synchronisation of correction factors TODO!!!!
          Dalpha(iedge) = Dalpha(iedge) * alpha_ij
          
          ! Update the raw antidiffusive fluxes
          F_ij = alpha_ij * F_ij
          F_ji = alpha_ij * F_ji

          ! Compute the sums of positive/negative antidiffusive increments
          Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
          Dpp(:,j) = Dpp(:,j)+max(0.0_DP, F_ji)
          Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
          Dpm(:,j) = Dpm(:,j)+min(0.0_DP, F_ji)
        end do
      end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      deallocate(DtransformedPrelFluxesAtEdge)

    end subroutine doPreADIncrementsTransformed
    
    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! without transformation

    subroutine doBounds(IverticesAtEdge, NEDGE, NEQ, NVAR, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Diff
      integer :: iedge,i,j


      ! Clear Q`s
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute solution difference
        Diff = Dx(:,j)-Dx(:,i)

        ! Compute the distance to a local extremum
        ! of the predicted solution
        Dqp(:,i) = max(Dqp(:,i), Diff)
        Dqp(:,j) = max(Dqp(:,j),-Diff)
        Dqm(:,i) = min(Dqm(:,i), Diff)
        Dqm(:,j) = min(Dqm(:,j),-Diff)
      end do

    end subroutine doBounds

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! which is transformed to a user-defined set of variables

    subroutine doBoundsTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(inout) :: Dqp,Dqm

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DtransformedDataAtEdge

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
      

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedDataAtEdge(NVARtransformed,GFSYS_NEDGESIM))
      
      ! Clear Q`s
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
        end do

        ! Use callback function to compute transformed differences
        call fcb_calcDiffTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            DtransformedDataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
        
          ! Compute the distance to a local extremum of the predicted solution
          Dqp(:,i) = max(Dqp(:,i), DtransformedDataAtEdge(:,idx))
          Dqp(:,j) = max(Dqp(:,j),-DtransformedDataAtEdge(:,idx))
          Dqm(:,i) = min(Dqm(:,i), DtransformedDataAtEdge(:,idx))
          Dqm(:,j) = min(Dqm(:,j),-DtransformedDataAtEdge(:,idx))
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedDataAtEdge)

    end subroutine doBoundsTransformed

    !**************************************************************
    ! Compute nodal correction factors without constraints

    subroutine doLimitNodal(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq


      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpp(:,ieq) .gt. AFCSTAB_EPSABS)
          Drp(:,ieq) = ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq))
        elsewhere
          Drp(:,ieq) = 1.0_DP
        end where
      end do

      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpm(:,ieq) .lt. -AFCSTAB_EPSABS)
          Drm(:,ieq) = ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq))
        elsewhere
          Drm(:,ieq) = 1.0_DP
        end where
      end do

    end subroutine doLimitNodal

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrained(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq


      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpp(:,ieq) .gt. AFCSTAB_EPSABS)
          Drp(:,ieq) = min(1.0_DP, ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq)))
        elsewhere
          Drp(:,ieq) = 1.0_DP
        end where
      end do

      ! Loop over all vertices
      do ieq = 1, NEQ
        where (dscale*Dpm(:,ieq) .lt. -AFCSTAB_EPSABS)
          Drm(:,ieq) = min(1.0_DP, ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq)))
        elsewhere
          Drm(:,ieq) = 1.0_DP
        end where
      end do

    end subroutine doLimitNodalConstrained

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewise(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F_ij,R_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Get precomputed raw antidiffusive fluxes
        F_ij = Dflux(:,iedge)

        ! Compute nodal correction factors
        where (F_ij .ge. 0.0_DP)
          R_ij = min(Drp(:,i),Drm(:,j))
        elsewhere
          R_ij = min(Drp(:,j),Drm(:,i))
        end where

!!$       REMARK: Numerical test demonstrate that this modification has no
!!$               significant influence on the order-of-accuracy
!!$
!!$        where (F_ij .gt. AFCSTAB_EPSABS)
!!$          R_ij = min(Drp(:,i),Drm(:,j))
!!$        elsewhere (F_ij .lt. -AFCSTAB_EPSABS)
!!$          R_ij = min(Drp(:,j),Drm(:,i))
!!$        elsewhere
!!$          R_ij = 1.0_DP
!!$        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(R_ij)
      end do

    end subroutine doLimitEdgewise

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    ! which are transformed to a user-defined set of variables
    ! priori to computing the correction factors

    subroutine doLimitEdgewiseTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux(:,IEDGEset:IEDGEmax),&
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Compute nodal correction factors
          R_ij = merge(Drp(:,i), Drm(:,i),&
                       DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
          R_ji = merge(Drp(:,j), Drm(:,j),&
                       DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)

!!$       REMARK: Numerical test demonstrate that this modification has no
!!$               significant influence on the order-of-accuracy
!!$
!!$          where (DtransformedFluxesAtEdge(:,1,idx) .gt. AFCSTAB_EPSABS)
!!$            R_ij = Drp(:,i)
!!$          elsewhere (DtransformedFluxesAtEdge(:,1,idx) .lt. -AFCSTAB_EPSABS)
!!$            R_ij = Drm(:,i)
!!$          elsewhere
!!$            R_ij = 1.0_DP
!!$          end where
!!$
!!$          where (DtransformedFluxesAtEdge(:,2,idx) .gt. AFCSTAB_EPSABS)
!!$            R_ji = Drp(:,j)
!!$          elsewhere (DtransformedFluxesAtEdge(:,2,idx) .lt. -AFCSTAB_EPSABS)
!!$            R_ji = Drm(:,j)
!!$          elsewhere
!!$            R_ij = 1.0_DP
!!$          end where

          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
        end do
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      
    end subroutine doLimitEdgewiseTransformed

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrained(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F1_ij,F2_ij,R_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

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

    end subroutine doLimitEdgewiseConstrained

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes which are transformed
    ! to a user-defined set of variables priori to computing the
    ! correction factors

    subroutine doLimitEdgewiseConstrainedTransformed(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes1AtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes2AtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes1AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes2AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      
      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux1(:,IEDGEset:IEDGEmax),&
            DtransformedFluxes1AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
            Dflux2(:,IEDGEset:IEDGEmax),&
            DtransformedFluxes2AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1
          
          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

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
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxes1AtEdge)
      deallocate(DtransformedFluxes2AtEdge)
      
    end subroutine doLimitEdgewiseConstrainedTransformed

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrect(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, dscale, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Correct antidiffusive flux
        F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

        ! Apply limited antidiffusive fluxes
        Dy(:,i) = Dy(:,i) + F_ij
        Dy(:,j) = Dy(:,j) - F_ij
      end do
      
    end subroutine doCorrect

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

    subroutine doCorrectScaleByMass(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, dscale, ML, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha,ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Correct antidiffusive flux
        F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

        ! Apply limited antidiffusive fluxes
        Dy(:,i) = Dy(:,i) + F_ij/ML(i)
        Dy(:,j) = Dy(:,j) - F_ij/ML(j)
      end do

    end subroutine doCorrectScaleByMass

  end subroutine gfsys_buildDivVecFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildFluxFCTBlock(rafcstab, rx, fcb_calcFlux2_sim,&
      theta, tstep, dscale, binit, rmatrix, ry, rcollection)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for FEM-FCT schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
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
    ! TRUE  : assemble the initial antidiffusive flux
    ! FALSE : assemble the antidiffusive flux using some initial values
    logical, intent(in) :: binit

    ! callback functions to compute antidiffusive fluxes
    include 'intf_gfsyscallback.inc'

    ! OPTIONAL: mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative
    type(t_vectorBlock), intent(in), optional :: ry
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dflux0,p_Dflux,p_Dalpha
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge


    ! Check if block vector contains only one block
    if (rx%nblocks .eq. 1) then
      if (present(ry)) then
        call gfsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            fcb_calcFlux2_sim, theta, tstep, dscale, binit,&
            rmatrix, ry%RvectorBlock(1), rcollection)
      else
        call gfsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            fcb_calcFlux2_sim, theta, tstep, dscale, binit,&
            rmatrix, rcollection=rcollection)
      end if
      return
    end if
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0) .and.&
        (rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_FEMFCT_MASS)) then
      call output_line('Stabilisation does not provide edge-based data structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix (if present)
    if (present(rmatrix)) then
      if ((rafcstab%NEQ .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE .ne. int(0.5*(rmatrix%NA-rmatrix%NEQ),I32)) .or.&
          (rafcstab%cmatrixFormat .ne. rmatrix%cmatrixFormat)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
        call sys_halt()
      end if
    end if

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsysbl_getbase_double(rx, p_Dx)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

      if (rafcstab%bprelimiting) then
        call output_line('Prelimiting for classical FEM-FCT has not been implemented yet!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
        call sys_halt()
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

      ! Set pointers to antidiffusive fluxes
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)

      ! Check if the amount of rejected antidiffusion should be
      ! included in the initial raw antidiffusive fluxes
      if (.not.binit .and.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_ITERATIVE) then

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
          call output_line('Stabilisation does not provide correction factors',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
          call sys_halt()
        end if

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
          call sys_halt()
        end if

        ! Set pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

        ! Subtract amount of rejected antidiffusion
        call doCombineFluxes(rafcstab%NVAR, rafcstab%NEDGE,&
            -1.0_DP, p_Dalpha, p_Dflux, p_Dflux0)
      end if

      ! Do we have a consistent mass matrix?
      if (present(rmatrix) .and. present(ry)) then

        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------
        
        ! Set pointers
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        call lsysbl_getbase_double(ry, p_Dy)

        ! Are we in the first step?
        if (binit) then
          ! Assemble total raw-antidiffusive fluxes for first step
          ! (without the contribution of the time derivative)
          call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)

          ! Assemble explicit part of raw-antidiffusive fluxes
          call lalg_vectorLinearComb(p_Dflux, p_Dflux0, 1.0_DP-theta, 0.0_DP)
          
          ! Apply mass antidiffusion to explicit fluxes
          call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_Dmatrix, p_Dy, -dscale/tstep, p_Dflux0)
        else
          ! Assemble implicit part of raw-antidiffusive fluxes
          call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, theta*dscale, p_Dflux)

          ! Apply mass antidiffusion to implicit fluxes
          call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_Dmatrix, p_Dy, dscale/tstep, p_Dflux)
          
          ! Assemble total raw-antidiffusive fluxes
          call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
        end if

      else

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Assemble raw-antidiffusive fluxes
        call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)
        
        ! Combine explicit and implicit fluxes
        if (binit) then
          call lalg_copyVector(p_Dflux, p_Dflux0)
        elseif (1.0_DP-theta .gt. SYS_EPSREAL) then
          call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP-theta, theta)
        elseif (theta .gt. SYS_EPSREAL) then
          call lalg_scaleVector(p_Dflux, theta)
        else
          call lalg_clearVector(p_Dflux)
        end if

      end if
      
      ! Do we have to store the initial fluxes separately?
      if (binit .and.&
          (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT))&
          call lsyssc_copyVector(rafcstab%p_rvectorFlux,&
          rafcstab%p_rvectorFluxPrel)

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_FEMFCT_LINEARISED)

      !-------------------------------------------------------------------------
      ! Linearised FEM-FCT algorithm
      !-------------------------------------------------------------------------

      ! Do we have to use the consistent mass matrix?
      if (present(rmatrix) .and. present(ry)) then

        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        call lsysbl_getbase_double(ry, p_Dy)

        ! Assemble raw-antidiffusive fluxes
        call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)

        ! Apply mass antidiffusion to antidiffusive fluxes
        call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_Dmatrix, p_Dy, dscale, p_Dflux)
      end if

      if (.not.(present(rmatrix)) .or. rafcstab%bprelimiting) then

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointer
        if (rafcstab%bprelimiting) then
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux)
        else
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        end if

        ! Assemble raw-antidiffusive fluxes
        call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)
        
        ! Prelimiting is only necessary if the consistent mass matrix
        ! is built into the raw antidiffusive fluxes. However, if the
        ! switch for prelimiting was not set to .false. and no mass
        ! antidiffusion is built into the fluxes we can simply copy it
        if (.not.(present(rmatrix)) .and. rafcstab%bprelimiting)&
            call lsyssc_copyVector(&
            rafcstab%p_rvectorFluxPrel, rafcstab%p_rvectorFlux)
      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_FEMFCT_MASS)

      !-------------------------------------------------------------------------
      ! FEM-FCT algorithm for mass antidiffusion
      !-------------------------------------------------------------------------

      if (present(rmatrix)) then

        ! Set pointers
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)

        ! Clear vector and assemble antidiffusive fluxes
        call lalg_clearVector(p_Dflux)
        call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_Dmatrix, p_Dx, dscale, p_Dflux)

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
    ! Assemble raw antidiffusive fluxes (without
    ! contribution of the consistent mass matrix)

    subroutine doFluxes(IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtEdge, Dx, dscale, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      
      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
        end do

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux2_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            Dflux(:,IEDGEset:IEDGEmax), rcollection)
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)

    end subroutine doFluxes

    !**************************************************************
    ! Assemble raw antidiffusive mass fluxes

#ifndef USE_OPENMP
    pure &
#endif

    subroutine doMassFluxes(IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixData, Dx, dscale, Dflux)
      
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DmatrixData
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      integer :: iedge,ij,i,j
      
      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)

        ! Compute the raw antidiffusives fluxes
        Dflux(:,iedge) = Dflux(:,iedge) + dscale*DmatrixData(ij)*(Dx(i,:)-Dx(j,:))
      end do

    end subroutine doMassFluxes

    !**************************************************************
    ! Combine two fluxes: flux2 := flux2+dscale*alpha*flux2

#ifndef USE_OPENMP
    pure &
#endif

    subroutine doCombineFluxes(NVAR, NEDGE, dscale, Dalpha, Dflux1, Dflux2)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR, NEDGE

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux2

      ! local variables
      integer :: iedge

      do iedge = 1, NEDGE
        Dflux2(:,iedge) = Dflux2(:,iedge) +&
            dscale * Dalpha(iedge) * Dflux1(:,iedge)
      end do

    end subroutine doCombineFluxes

  end subroutine gfsys_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildFluxFCTScalar(rafcstab, rx, fcb_calcFlux2_sim,&
      theta, tstep, dscale, binit, rmatrix, ry, rcollection)

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
    ! TRUE  : assemble the initial antidiffusive flux
    ! FALSE : assemble the antidiffusive flux using some initial values
    logical, intent(in) :: binit

    ! callback functions to compute antidiffusive fluxes
    include 'intf_gfsyscallback.inc'

    ! OPTIONAL: mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative
    type(t_vectorScalar), intent(in), optional :: ry
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dflux0,p_Dflux,p_Dalpha
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge


    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0) .and.&
        (rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_FEMFCT_MASS)) then
      call output_line('Stabilisation does not provide edge-based data structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix (if present)
    if (present(rmatrix)) then
      if ((rafcstab%NEQ .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE .ne. int(0.5*(rmatrix%NA-rmatrix%NEQ),I32)) .or.&
          (rafcstab%cmatrixFormat .ne. rmatrix%cmatrixFormat)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
        call sys_halt()
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

      if (rafcstab%bprelimiting) then
        call output_line('Prelimiting for classical FEM-FCT has not been implemented yet!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
        call sys_halt()
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

      ! Get pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
      
      ! Check if the amount of rejected antidiffusion should be
      ! included in the initial raw antidiffusive fluxes
      if (.not.binit .and.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_ITERATIVE) then

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
          call output_line('Stabilisation does not provide correction factors',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
          call sys_halt()
        end if

        ! Set pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

        ! Subtract amount of rejected antidiffusion
        call doCombineFluxes(rafcstab%NVAR, rafcstab%NEDGE,&
            -1.0_DP, p_Dalpha, p_Dflux, p_Dflux0)
      end if

      ! Do we have to use the consistent mass matrix?
      if (present(rmatrix) .and. present(ry)) then
        
        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------
        
        ! Set pointer for consistent mass matrix
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        call lsyssc_getbase_double(ry, p_Dy)
        
        ! Are we in the first step?
        if (binit) then
          ! Assemble total raw-antidiffusive fluxes for first step
          ! (without the contribution of the time derivative)
          call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)

          ! Assemble explicit part of raw-antidiffusive fluxes
          call lalg_vectorLinearComb(p_Dflux, p_Dflux0, 1.0_DP-theta, 0.0_DP)

          ! Apply mass antidiffusion to explicit fluxes
          call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_Dmatrix, p_Dy, -dscale/tstep, p_Dflux0)
        else
          ! Assemble implicit part of raw-antidiffusive fluxes
          call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, theta*dscale, p_Dflux)

          ! Apply mass antidiffusion to implicit fluxes
          call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_Dmatrix, p_Dy, dscale/tstep, p_Dflux)

          ! Assemble total raw-antidiffusive fluxes
          call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
        end if
        
      else

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Assemble raw-antidiffusive fluxes
        call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)
        
        ! Combine explicit and implicit fluxes
        if (binit) then
          call lalg_copyVector(p_Dflux, p_Dflux0)
        elseif (1.0_DP-theta .gt. SYS_EPSREAL) then
          call lalg_vectorLinearComb(&
              p_Dflux0, p_Dflux, 1.0_DP-theta, theta)
        elseif (theta .gt. SYS_EPSREAL) then
          call lalg_scaleVector(p_Dflux, theta)
        else
          call lalg_clearVector(p_Dflux)
        end if

      end if

      ! Do we have to store the initial fluxes separately?
      if (binit .and.&
          (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT))&
          call lsyssc_copyVector(rafcstab%p_rvectorFlux,&
          rafcstab%p_rvectorFluxPrel)
      
      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_FEMFCT_LINEARISED)

      !-------------------------------------------------------------------------
      ! Linearised FEM-FCT algorithm
      !-------------------------------------------------------------------------

      ! Do we have to use the consistent mass matrix?
      if (present(rmatrix) .and. present(ry)) then
        
        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        call lsyssc_getbase_double(ry, p_Dy)

        ! Assemble raw-antidiffusive fluxes
        call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)

        ! Apply mass antidiffusion to antidiffusive fluxes
        call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_Dmatrix, p_Dy, dscale, p_Dflux)
      end if

      if (.not.(present(rmatrix)) .or. rafcstab%bprelimiting) then

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointer
        if (rafcstab%bprelimiting) then
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux)
        else
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        end if

        ! Assemble raw-antidiffusive fluxes
        call doFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dflux)
        
        ! Prelimiting is only necessary if the consistent mass matrix
        ! is built into the raw antidiffusive fluxes. However, if the
        ! switch for prelimiting was not set to .false. and no mass
        ! antidiffusion is built into the fluxes we can simply copy it
        if (.not.(present(rmatrix)) .and. rafcstab%bprelimiting)&
            call lsyssc_copyVector(&
            rafcstab%p_rvectorFluxPrel, rafcstab%p_rvectorFlux)
      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)

      
    case (AFCSTAB_FEMFCT_MASS)

      !-------------------------------------------------------------------------
      ! FEM-FCT algorithm for mass antidiffusion
      !-------------------------------------------------------------------------

      if (present(rmatrix)) then

        ! Set pointers
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)

        ! Clear vector and assemble antidiffusive fluxes
        call lalg_clearVector(p_Dflux)
        call doMassFluxes(p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_Dmatrix, p_Dx, dscale, p_Dflux)
        
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
    ! Assemble raw antidiffusive fluxes (without
    ! contribution of the consistent mass matrix)

    subroutine doFluxes(IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtEdge, Dx, dscale, Dflux)
      
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      
      ! local variables
      integer :: i,j,idx,iedge,IEDGEset,IEDGEmax


      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edges
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
          DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
        end do

        ! Use callback function to compute internodal fluxes
        call fcb_calcFlux2_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
            IverticesAtEdge(:,IEDGEset:IEDGEmax), dscale,&
            Dflux(:,IEDGEset:IEDGEmax), rcollection)
      end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)

    end subroutine doFluxes

    !**************************************************************
    ! Assemble raw antidiffusive mass fluxes

#ifndef USE_OPENMP
    pure &
#endif

    subroutine doMassFluxes(IverticesAtEdge, NEDGE, NEQ, NVAR,&
        DmatrixData, Dx,dscale, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DmatrixData
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      integer :: iedge,ij,i,j
      
      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)

        ! Compute the raw antidiffusives fluxes
        Dflux(:,iedge) = Dflux(:,iedge) + dscale*DmatrixData(ij)*(Dx(:,i)-Dx(:,j))
      end do

    end subroutine doMassFluxes

    !**************************************************************
    ! Combine two fluxes: flux2 := flux2+dscale*alpha*flux2

#ifndef USE_OPENMP
    pure &
#endif

    subroutine doCombineFluxes(NVAR, NEDGE, dscale, Dalpha, Dflux1, Dflux2)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR, NEDGE

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux2

      ! local variables
      integer :: iedge

      do iedge = 1, NEDGE
        Dflux2(:,iedge) = Dflux2(:,iedge) +&
            dscale * Dalpha(iedge) * Dflux1(:,iedge)
      end do

    end subroutine doCombineFluxes

  end subroutine gfsys_buildFluxFCTScalar

end module groupfemsystem

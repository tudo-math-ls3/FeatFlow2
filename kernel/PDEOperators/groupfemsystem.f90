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
!# gfsc_buildDivVectorXXX. They can be used to update/initialise the
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
!# 2.) gfsys_isMatrixCompatible = gfsys_isMatrixCompatibleScalar /
!#                                gfsys_isMatrixCompatibleBlock
!#     -> checks wether a matrix and a stabilisation structure are compatible
!#
!# 3.) gfsys_isVectorCompatible = gfsys_isVectorCompatibleScalar /
!#                                gfsys_isVectorCompatibleBlock
!#     -> checks wether a vector and a stabilisation structure are compatible
!#
!# 4.) gfsys_buildDivOperator = gfsys_buildDivOperatorScalar /
!#                              gfsys_buildDivOperatorBlock
!#     -> assembles the global operator that results from the discretisation
!#        of the divergence term $div(F)$ by means of the Galerkin method
!#        and some sort of artificial dissipation (if required)
!#
!# 5.) gfsys_buildDivVector = gfsys_buildVecScalar /
!#                            gfsys_buildVecBlock
!#     -> assembles the divergence vector
!#
!# 6.) gfsys_buildDivVectorTVD = gfsys_buildVecTVDScalar /
!#                               gfsys_buildVecTVDBlock
!#     -> assembles the divergence term for FEM-TVD stabilisation
!#
!# 7.) gfsys_buildDivVectorFCT = gfsys_buildVecFCTScalar /
!#                               gfsys_buildVecFCTBlock
!#     -> assembles the divergence term for FEM-FCT stabilisation
!#
!# 8.) gfsys_buildFluxFCT = gfsys_buildFluxFCTScalar /
!#                          gfsys_buildFluxFCTBlock
!#     -> assembles the raw antidiffusive flux for FEM-FCT stabilisation
!#
!# The following internal routines are available:
!#
!# 1.) gfsys_getbase_double
!#     -> assigns the pointers of an array to a given block matrix
!#
!# 2.) gfsys_getbase_single
!#     -> assigns the pointers of an array to a given block matrix
!#
!# 3.) gfsys_getbase_int
!#     -> assigns the pointers of an array to a given block matrix
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
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use storage
  use triangulation

  implicit none

  private

  public :: gfsys_initStabilisation
  public :: gfsys_isMatrixCompatible
  public :: gfsys_isVectorCompatible
  public :: gfsys_buildDivOperator
  public :: gfsys_buildDivVector
  public :: gfsys_buildDivVectorTVD
  public :: gfsys_buildDivVectorFCT
  public :: gfsys_buildFluxFCT

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

!<types>
!<typeblock>

  ! An auxiliary array that is used to store the content of all
  ! scalar submatrices simultaneously. Note that this derived type
  ! is PRIVATE and cannot be accessed from outside of this module.
  private :: t_array
  type t_array

    ! Pointer to the double-valued matrix data
    real(DP), dimension(:), pointer :: Da

    ! Pointer to the single-valued matrix data
    real(SP), dimension(:), pointer :: Fa

    ! Pointer to the integer-valued matrix data
    integer, dimension(:), pointer :: Ia
  end type t_array
!</typeblock>
!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface gfsys_initStabilisation
    module procedure gfsys_initStabilisationScalar
    module procedure gfsys_initStabilisationBlock
  end interface

  interface gfsys_isMatrixCompatible
    module procedure gfsys_isMatrixCompatibleScalar
    module procedure gfsys_isMatrixCompatibleBlock
  end interface

  interface gfsys_isVectorCompatible
    module procedure gfsys_isVectorCompatibleScalar
    module procedure gfsys_isVectorCompatibleBlock
  end interface

  interface gfsys_buildDivOperator
     module procedure gfsys_buildDivOperatorScalar
     module procedure gfsys_buildDivOperatorBlock
  end interface

  interface gfsys_buildDivVector
    module procedure gfsys_buildVecScalar
    module procedure gfsys_buildVecBlock
  end interface

  interface gfsys_buildDivVectorTVD
    module procedure gfsys_buildVecTVDScalar
    module procedure gfsys_buildVecTVDBlock
  end interface

  interface gfsys_buildDivVectorFCT
    module procedure gfsys_buildVecFCTScalar
    module procedure gfsys_buildVecFCTBlock
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
      rafcstab, NVARtransformed)

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
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: i


    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and. &
        (rmatrixBlockTemplate%nblocksPerRow .eq. 1)) then
      call gfsys_initStabilisationScalar(&
          rmatrixBlockTemplate%RmatrixBlock(1,1),&
          rafcstab, NVARtransformed)
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

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end do


    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationBlock', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end do

      ! We need 3 edgewise vectors at most
      allocate(rafcstab%RedgeVectors(3))

      ! We need 1 edgewise vector for the correction factors
      call lsyssc_createVector(rafcstab%RedgeVectors(1),&
          rafcstab%NEDGE, 1, .false., ST_DOUBLE)

      ! We need 2 edgewise vectors for the raw antidiffusive fluxes
      call lsyssc_createVector(rafcstab%RedgeVectors(2),&
          rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%RedgeVectors(3),&
          rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      

    case (AFCSTAB_FEMFCT_LINEARISED)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationBlock', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end do

      ! We need 4 edgewise vectors at most
      allocate(rafcstab%RedgeVectors(4))
      
      ! We need 1 edgewise vector for the correction factors
      call lsyssc_createVector(rafcstab%RedgeVectors(1),&
          rafcstab%NEDGE, 1, .false., ST_DOUBLE)
      
      ! We need 1 edgewise vector for the raw antidiffusive fluxes.
      call lsyssc_createVector(rafcstab%RedgeVectors(2),&
          rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      
      ! If the raw antidiffusive fluxes should be prelimited
      ! then 1 additional edgewie vector is required
      if (rafcstab%bprelimiting) then
        call lsyssc_createVector(rafcstab%RedgeVectors(4),&
            rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      end if


    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end select

    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED

  end subroutine gfsys_initStabilisationBlock

  ! *****************************************************************************

  subroutine gfsys_initStabilisationScalar(rmatrixTemplate, rafcstab,&
      NVARtransformed)

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
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: i


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

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end do


    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationScalar', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end do
      
      ! We need 3 edgewise vectors at most
      allocate(rafcstab%RedgeVectors(3))

      ! We need 1 edgewise vector for the correction factors
      call lsyssc_createVector(rafcstab%RedgeVectors(1),&
          rafcstab%NEDGE, 1, .false., ST_DOUBLE)

      ! We need 2 edgewise vectors for the raw antidiffusive fluxes
      call lsyssc_createVector(rafcstab%RedgeVectors(2),&
          rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      call lsyssc_createVector(rafcstab%RedgeVectors(3),&
          rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)


    case (AFCSTAB_FEMFCT_LINEARISED)

      ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsys_initStabilisationScalar', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end do
      
      ! We need 4 edgewise vectors at most
      allocate(rafcstab%RedgeVectors(4))
      
      ! We need 1 edgewise vector for the correction factors
      call lsyssc_createVector(rafcstab%RedgeVectors(1),&
          rafcstab%NEDGE, 1, .false., ST_DOUBLE)
      
      ! We need 1 edgewise vector for the raw antidiffusive fluxes.
      call lsyssc_createVector(rafcstab%RedgeVectors(2),&
          rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      
      ! If the raw antidiffusive fluxes should be prelimited
      ! then 1 additional edgewie vector is required
      if (rafcstab%bprelimiting) then
        call lsyssc_createVector(rafcstab%RedgeVectors(4),&
            rafcstab%NEDGE, rafcstab%NVAR, .false., ST_DOUBLE)
      end if
      

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationScalar')
      call sys_halt()
    end select

    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED

  end subroutine gfsys_initStabilisationScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_isMatrixCompatibleBlock(rafcstab, rmatrixBlock, bcompatible)

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
    if ((rmatrixBlock%nblocksPerCol .eq. 1) .and. &
        (rmatrixBlock%nblocksPerRow .eq. 1)) then
      call gfsys_isMatrixCompatibleScalar(rafcstab,&
          rmatrixBlock%RmatrixBlock(1,1), bcompatible)
      return
    end if

    ! Check that number of columns equans number of rows
    if (rmatrixBlock%nblocksPerCol .ne. &
        rmatrixBlock%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isMatrixCompatibleBlock')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrixBlock%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Block matrix must have group structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isMatrixCompatibleBlock')
        call sys_halt()
      end if
    end if

    ! Matrix/operator must have the same size
    if (rafcstab%NVAR  .ne. rmatrixBlock%nblocksPerCol .or.&
        rafcstab%NEQ   .ne. rmatrixBlock%RmatrixBlock(1,1)%NEQ  .or.&
        rafcstab%NEDGE .ne. int(0.5*(rmatrixBlock%RmatrixBlock(1,1)%NA-&
                                     rmatrixBlock%RmatrixBlock(1,1)%NEQ))) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfssy_isMatrixCompatibleBlock')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isMatrixCompatibleBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_isMatrixCompatibleScalar(rafcstab, rmatrix, bcompatible)

!<description>
    ! This subroutine checks whether a scalar matrix and a discrete
    ! stabilisation structure are compatible to each other,
    ! i.e. if they share the same structure, size and so on.
    !
    ! Note that the matrix is required as scalar matrix. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

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

    ! Matrix/operator must have the same size
    if (rafcstab%NEQ   .ne. rmatrix%NEQ  .or.&
        rafcstab%NVAR  .ne. rmatrix%NVAR .or.&
        rafcstab%NEDGE .ne. int(0.5*(rmatrix%NA-rmatrix%NEQ))) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfssy_isMatrixCompatibleScalar')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isMatrixCompatibleScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_isVectorCompatibleBlock(rafcstab, rvectorBlock, bcompatible)

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
      call gfsys_isVectorCompatibleScalar(rafcstab,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isVectorCompatibleBlock')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isVectorCompatibleBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_isVectorCompatibleScalar(rafcstab, rvector, bcompatible)

!<description>
    ! This subroutine checks whether a vector and a stabilisation
    ! structure are compatible to each other, i.e., share the same
    ! structure, size and so on.
    !
    ! Note that the vector is required as scalar vector. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! scalar vector
    type(t_vectorScalar), intent(in) :: rvector

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

    ! Matrix/operator must have the same size
    if (rafcstab%NEQ   .ne. rvector%NEQ .or.&
        rafcstab%NVAR  .ne. rvector%NVAR) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_isVectorCompatibleScalar')
        call sys_halt()
      end if
    end if
  end subroutine gfsys_isVectorCompatibleScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorBlock(RcoeffMatrices, rafcstab, rx,&
      fcb_calcMatrixDiagonal, fcb_calcMatrix, dscale, bclear, rdivMatrix)

!<description>
    ! This subroutine assembles the discrete divergence operator.
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! global transport operator
    type(t_matrixBlock), intent(inout) :: rdivMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_array), dimension(rx%nblocks,rx%nblocks)  :: rarray
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: ndim
    logical :: bisFullMatrix


    ! Check if block vector contains only one block and if
    ! global operator is stored in interleave format.
    if ((rx%nblocks .eq. 1) .and.&
        (rdivMatrix%nblocksPerCol .eq. 1) .and. &
        (rdivMatrix%nblocksPerRow .eq. 1)) then
      call gfsys_buildDivOperatorScalar(RcoeffMatrices, rafcstab,&
          rx%RvectorBlock(1), fcb_calcMatrixDiagonal, fcb_calcMatrix,&
          dscale, bclear, rdivMatrix%RmatrixBlock(1,1))
      return
    end if

    ! Check if block matrix exhibits group structure
    if (rdivMatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear matrix?
    if (bclear) call lsysbl_clearMatrix(rdivMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call gfsys_getbase_double(rdivMatrix, rarray, bisFullMatrix)
    call lsysbl_getbase_double(rx, p_Dx)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end select

    ! What kind of matrix are we?
    select case(RcoeffMatrices(1)%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------

      ! Set diagonal pointer
      call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kdiagonal)

      ! What type of matrix are we?
      if (bisFullMatrix) then

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_Dx, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat79_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat79_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, rarray)
        end select

      else   ! bisFullMatrix == no

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79Diag_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_Dx, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat79Diag_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat79Diag_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, rarray)
        end select

      end if   ! bisFullMatrix


    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set diagonal pointer
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)

      ! What type of matrix are we?
      if (bisFullMatrix) then

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_Dx, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat79_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat79_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, rarray)
        end select

      else   ! bisFullMatrix == no

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79Diag_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_Dx, dscale, rarray)
        case (NDIM2D)
          call doOperatorMat79Diag_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, rarray)
        case (NDIM3D)
          call doOperatorMat79Diag_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, rarray)
        end select

      end if   ! bisFullMatrix


    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag_1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE, NEQ,NVAR

      type(t_array), dimension(:,:), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: D_ij,K_ij,K_ji,Dx_i,Dx_j
      integer :: iedge,ii,ij,ji,jj,i,j,ivar


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij,Dx_i,ivar)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii)

        ! Get solution values at node
        Dx_i = Dx(i,:)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx_i, C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        do ivar = 1, NVAR
          K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + K_ij(ivar)
        end do
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Get solution values at node
        Dx_j = Dx(j,:)

        ! Compute matrices
        call fcb_calcMatrix(Dx_i, Dx_j, C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Assemble the global operator
        do ivar = 1, NVAR
          K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)              + D_ij(ivar)
          K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) + K_ij(ivar) - D_ij(ivar)
          K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + K_ji(ivar) - D_ij(ivar)
          K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)              + D_ij(ivar)
        end do
      end do
    end subroutine doOperatorMat79Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag_2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE, NEQ,NVAR

      type(t_array), dimension(:,:), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: D_ij,K_ij,K_ji,Dx_i,Dx_j
      integer :: iedge,ii,ij,ji,jj,i,j,ivar


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij,Dx_i,ivar)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii)

        ! Get solution values at node
        Dx_i = Dx(i,:)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx_i, C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        do ivar = 1, NVAR
          K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + K_ij(ivar)
        end do
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Get solution values at node
        Dx_j = Dx(j,:)

        ! Compute matrices
        call fcb_calcMatrix(Dx_i, Dx_j, C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Assemble the global operator
        do ivar = 1, NVAR
          K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)              + D_ij(ivar)
          K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) + K_ij(ivar) - D_ij(ivar)
          K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + K_ji(ivar) - D_ij(ivar)
          K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)              + D_ij(ivar)
        end do
      end do
    end subroutine doOperatorMat79Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag_3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NVAR

      type(t_array), dimension(:,:), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: D_ij,K_ij,K_ji,Dx_i,Dx_j
      integer :: iedge,ii,ij,ji,jj,i,j,ivar


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij,Dx_i,ivar)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii); C_ij(3) = DcoeffZ(ii)

        ! Get solution values at node
        Dx_i = Dx(i,:)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx_i, C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        do ivar = 1, NVAR
          K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii) + K_ij(ivar)
        end do
      end do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Get solution values at node
        Dx_j = Dx(j,:)

        ! Compute matrices
        call fcb_calcMatrix(Dx_i, Dx_j, C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Assemble the global operator
        do ivar = 1, NVAR
          K(ivar,ivar)%DA(ii) = K(ivar,ivar)%DA(ii)              + D_ij(ivar)
          K(ivar,ivar)%DA(ij) = K(ivar,ivar)%DA(ij) + K_ij(ivar) - D_ij(ivar)
          K(ivar,ivar)%DA(ji) = K(ivar,ivar)%DA(ji) + K_ji(ivar) - D_ij(ivar)
          K(ivar,ivar)%DA(jj) = K(ivar,ivar)%DA(jj)              + D_ij(ivar)
        end do
      end do
    end subroutine doOperatorMat79Diag_3D


    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79_1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NVAR

      type(t_array), dimension(:,:), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: D_ij,K_ij,K_ji
      real(DP), dimension(NVAR) :: Dx_i,Dx_j
      integer :: iedge,ii,ij,ji,jj,i,j,ivar,jvar,idx


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij,Dx_i,ivar)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii)

        ! Get solution values at node
        Dx_i = Dx(i,:)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx_i, C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        do ivar = 1, NVAR
          do jvar = 1, NVAR
            idx = NVAR*(ivar-1)+jvar
            K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + K_ij(idx)
          end do
        end do
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Get solution values at node
        Dx_j = Dx(j,:)

        ! Compute matrices
        call fcb_calcMatrix(Dx_i, Dx_j, C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        do ivar = 1, NVAR
          do jvar = 1, NVAR
            idx = NVAR*(ivar-1)+jvar
            K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)             + D_ij(idx)
            K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) + K_ij(idx) - D_ij(idx)
            K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + K_ji(idx) - D_ij(idx)
            K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)             + D_ij(idx)
          end do
        end do
      end do
    end subroutine doOperatorMat79_1D


    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79_2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NVAR

      type(t_array), dimension(:,:), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: D_ij,K_ij,K_ji
      real(DP), dimension(NVAR) :: Dx_i,Dx_j
      integer :: iedge,ii,ij,ji,jj,i,j,ivar,jvar,idx


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij,Dx_i,ivar)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ij)

        ! Get solution values at node
        Dx_i = Dx(i,:)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx_i, C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        do ivar = 1, NVAR
          do jvar = 1, NVAR
            idx = NVAR*(ivar-1)+jvar
            K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + K_ij(idx)
          end do
        end do
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Get solution values at node
        Dx_j = Dx(j,:)

        ! Compute matrices
        call fcb_calcMatrix(Dx_i, Dx_j, C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        do ivar = 1, NVAR
          do jvar = 1, NVAR
            idx = NVAR*(ivar-1)+jvar
            K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)             + D_ij(idx)
            K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) + K_ij(idx) - D_ij(idx)
            K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + K_ji(idx) - D_ij(idx)
            K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)             + D_ij(idx)
          end do
        end do
      end do
    end subroutine doOperatorMat79_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79_3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, K)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NVAR

      type(t_array), dimension(:,:), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: D_ij,K_ij,K_ji
      real(DP), dimension(NVAR) :: Dx_i,Dx_j
      integer :: iedge,ii,ij,ji,jj,i,j,ivar,jvar,idx


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij,Dx_i,ivar)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii); C_ij(3) = DcoeffZ(ii)

        ! Get solution values at node
        Dx_i = Dx(i,:)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx_i, C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        do ivar = 1, NVAR
          do jvar = 1, NVAR
            idx = NVAR*(ivar-1)+jvar
            K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii) + K_ij(idx)
          end do
        end do
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Get solution values at node
        Dx_j = Dx(j,:)

        ! Compute matrices
        call fcb_calcMatrix(Dx_i, Dx_j, C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        do ivar = 1, NVAR
          do jvar = 1, NVAR
            idx = NVAR*(ivar-1)+jvar
            K(jvar,ivar)%DA(ii) = K(jvar,ivar)%DA(ii)             + D_ij(idx)
            K(jvar,ivar)%DA(ij) = K(jvar,ivar)%DA(ij) + K_ij(idx) - D_ij(idx)
            K(jvar,ivar)%DA(ji) = K(jvar,ivar)%DA(ji) + K_ji(idx) - D_ij(idx)
            K(jvar,ivar)%DA(jj) = K(jvar,ivar)%DA(jj)             + D_ij(idx)
          end do
        end do
      end do
    end subroutine doOperatorMat79_3D

  end subroutine gfsys_buildDivOperatorBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorScalar(RcoeffMatrices, rafcstab, rx,&
      fcb_calcMatrixDiagonal, fcb_calcMatrix, dscale, bclear, rdivMatrix)

!<description>
    ! This subroutine assembles the discrete divergence operator.
    ! Note that this routine requires the scalar matrices/vectors
    ! are stored in the interleave format.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! scalar solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute local matrices
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! scalar transport operator
    type(t_matrixScalar), intent(inout) :: rdivMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_DivOp,p_Dx
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: ndim


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rdivMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsyssc_getbase_double(rdivMatrix, p_DivOp)
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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end select

    ! What kind of matrix are we?
    select case(rdivMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7INTL)
      !-------------------------------------------------------------------------
      ! Matrix format 7 interleaved
      !-------------------------------------------------------------------------

      ! Set diagonal pointer
      call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kdiagonal)

      ! What type of matrix are we?
      select case(rdivMatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIX1)

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_Dx, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat79_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat79_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, p_DivOp)
        end select


      case (LSYSSC_MATRIXD)

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79Diag_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_Dx, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat79Diag_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat79Diag_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, p_DivOp)
        end select


      case DEFAULT
        call output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
        call sys_halt()
      end select


    case (LSYSSC_MATRIX9INTL)
      !-------------------------------------------------------------------------
      ! Matrix format 9 interleaved
      !-------------------------------------------------------------------------

      ! Set diagonal pointer
      call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)

      ! What type of matrix are we?
      select case(rdivMatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIX1)

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_Dx, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat79_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat79_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, p_DivOp)
        end select


      case (LSYSSC_MATRIXD)

        select case(ndim)
        case (NDIM1D)
          call doOperatorMat79Diag_1D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_Dx, dscale, p_DivOp)
        case (NDIM2D)
          call doOperatorMat79Diag_2D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_Dx, dscale, p_DivOp)
        case (NDIM3D)
          call doOperatorMat79Diag_3D(p_Kdiagonal, p_IverticesAtEdge,&
              rafcstab%NEDGE, RcoeffMatrices(1)%NEQ,&
              RcoeffMatrices(1)%NA, rx%NVAR,&
              p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, p_DivOp)
        end select


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
    ! Assemble block-diagonal divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag_1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NA, NVAR, DcoeffX, Dx, dscale, K)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: D_ij,K_ij,K_ji
      integer :: iedge,ii,ij,ji,jj,i,j


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx(:,i), C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        K(:,ii) = K(:,ii) + K_ij
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute matrices
        call fcb_calcMatrix(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        K(:,ii) = K(:,ii)        + D_ij
        K(:,ij) = K(:,ij) + K_ij - D_ij
        K(:,ji) = K(:,ji) + K_ji - D_ij
        K(:,jj) = K(:,jj)        + D_ij
      end do
    end subroutine doOperatorMat79Diag_1D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag_2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NA, NVAR, DcoeffX, DcoeffY, Dx, dscale, K)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: D_ij,K_ij,K_ji
      integer :: iedge,ii,ij,ji,jj,i,j


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx(:,i), C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        K(:,ii) = K(:,ii) + K_ij
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute matrices
        call fcb_calcMatrix(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        K(:,ii) = K(:,ii)        + D_ij
        K(:,ij) = K(:,ij) + K_ij - D_ij
        K(:,ji) = K(:,ji) + K_ji - D_ij
        K(:,jj) = K(:,jj)        + D_ij
      end do
    end subroutine doOperatorMat79Diag_2D


    !**************************************************************
    ! Assemble block-diagonal divergence operator K in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag_3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NA, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, K)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR

      real(DP), dimension(NVAR,NA), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: D_ij,K_ij,K_ji
      integer :: iedge,ii,ij,ji,jj,i,j


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii); C_ij(3) = DcoeffZ(ii)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx(:,i), C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        K(:,ii) = K(:,ii) + K_ij
      end do
      !$omp end parallel do


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute matrices
        call fcb_calcMatrix(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        K(:,ii) = K(:,ii)        + D_ij
        K(:,ij) = K(:,ij) + K_ij - D_ij
        K(:,ji) = K(:,ji) + K_ji - D_ij
        K(:,jj) = K(:,jj)        + D_ij
      end do
    end subroutine doOperatorMat79Diag_3D


    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79_1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NA, NVAR, DcoeffX, Dx, dscale, K)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: D_ij,K_ij,K_ji
      integer :: iedge,ii,ij,ji,jj,i,j


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx(:,i), C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        K(:,ii) = K(:,ii) + K_ij
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute matrices
        call fcb_calcMatrix(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        K(:,ii) = K(:,ii)        + D_ij
        K(:,ij) = K(:,ij) + K_ij - D_ij
        K(:,ji) = K(:,ji) + K_ji - D_ij
        K(:,jj) = K(:,jj)        + D_ij
      end do
    end subroutine doOperatorMat79_1D


    !**************************************************************
    ! Assemble divergence operator K in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79_2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NA, NVAR, DcoeffX, DcoeffY, Dx, dscale, K)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: D_ij,K_ij,K_ji
      integer :: iedge,ii,ij,ji,jj,i,j


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx(:,i), C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        K(:,ii) = K(:,ii) + K_ij
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute matrices
        call fcb_calcMatrix(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        K(:,ii) = K(:,ii)        + D_ij
        K(:,ij) = K(:,ij) + K_ij - D_ij
        K(:,ji) = K(:,ji) + K_ji - D_ij
        K(:,jj) = K(:,jj)        + D_ij
      end do
    end subroutine doOperatorMat79_2D


    !**************************************************************
    ! Assemble divergence operator K in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79_3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, NA, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, K)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ,NA,NVAR

      real(DP), dimension(NVAR*NVAR,NA), intent(inout) :: K

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: D_ij,K_ij,K_ji
      integer :: iedge,ii,ij,ji,jj,i,j


      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,K_ij)
      do i = 1, NEQ

        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ii); C_ij(2) = DcoeffY(ii); C_ij(3) = DcoeffZ(ii)

        ! Compute matrix
        call fcb_calcMatrixDiagonal(Dx(:,i), C_ij, i, dscale, K_ij)

        ! Apply matrix to the diagonal entry of the global operator
        K(:,ii) = K(:,ii) + K_ij
      end do
      !$omp end parallel do

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute matrices
        call fcb_calcMatrix(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                            i, j, dscale, K_ij, K_ji, D_ij)

        ! Apply matrices to the global operator
        K(:,ii) = K(:,ii)        + D_ij
        K(:,ij) = K(:,ij) + K_ij - D_ij
        K(:,ji) = K(:,ji) + K_ji - D_ij
        K(:,jj) = K(:,jj)        + D_ij
      end do
    end subroutine doOperatorMat79_3D

  end subroutine gfsys_buildDivOperatorScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildVecBlock(RcoeffMatrices, rafcstab, rx,&
      fcb_calcFlux, dscale, bclear, ry)

!<description>
    ! This subroutine assembles the divergence vector for block vectors.
    ! If the vector contains only one block, then the scalar
    ! counterpart of this routine is called with the scalar subvector.
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

    ! callback functions to compute local fluxes
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx,p_Dy
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim


    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildVecScalar(RcoeffMatrices, rafcstab,&
          rx%RvectorBlock(1), fcb_calcFlux, dscale, bclear,&
          ry%RvectorBlock(1))
      return
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecBlock')
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
        call doDivVectorMat79_1D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
            p_DcoeffX, p_Dx, dscale, p_Dy)
      case (NDIM2D)
        call doDivVectorMat79_2D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
            p_DcoeffX, p_DcoeffY, p_Dx, dscale, p_Dy)
      case (NDIM3D)
        call doDivVectorMat79_3D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%nblocks,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, p_Dy)
      end select

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doDivVectorMat79_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale, Dy)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: Dx_i,Dx_j,F_ij,F_ji
      integer :: iedge,ij,ji,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute the fluxes
        call fcb_calcFlux(Dx_i, Dx_j, C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)+F_ji
      end do
    end subroutine doDivVectorMat79_1D


    !**************************************************************
    ! Assemble divergence vector in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doDivVectorMat79_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale, Dy)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: Dx_i,Dx_j,F_ij,F_ji
      integer :: iedge,ij,ji,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute the fluxes
        call fcb_calcFlux(Dx_i, Dx_j, C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)+F_ji
      end do
    end subroutine doDivVectorMat79_2D


    !**************************************************************
    ! Assemble divergence vector in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doDivVectorMat79_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, Dy)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: Dx_i,Dx_j,F_ij,F_ji
      integer :: iedge,ij,ji,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute the fluxes
        call fcb_calcFlux(Dx_i, Dx_j, C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)+F_ji
      end do
    end subroutine doDivVectorMat79_3D
  end subroutine gfsys_buildVecBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildVecScalar(RcoeffMatrices, rafcstab, rx,&
      fcb_calcFlux, dscale, bclear, ry)

!<description>
    ! This subroutine assembles the divergence vector. Note that the
    ! vectors are required as scalar vectors which are stored in the
    ! interleave format.
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

    ! callback functions to compute local fluxes
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx,p_Dy
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecScalar')
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
        call doDivVectorMat79_1D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%NVAR,&
            p_DcoeffX, p_Dx, dscale, p_Dy)
      case (NDIM2D)
        call doDivVectorMat79_2D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%NVAR,&
            p_DcoeffX, p_DcoeffY, p_Dx, dscale, p_Dy)
      case (NDIM3D)
        call doDivVectorMat79_3D(p_IverticesAtEdge,&
            rafcstab%NEDGE, RcoeffMatrices(1)%NEQ, rx%NVAR,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, dscale, p_Dy)
      end select

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doDivVectorMat79_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale, Dy)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: iedge,ij,ji,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)+F_ji
      end do
    end subroutine doDivVectorMat79_1D


    !**************************************************************
    ! Assemble divergence vector in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doDivVectorMat79_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale, Dy)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: iedge,ij,ji,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)+F_ji
      end do
    end subroutine doDivVectorMat79_2D


    !**************************************************************
    ! Assemble divergence vector in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doDivVectorMat79_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, Dy)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij,F_ji
      integer :: iedge,ij,ji,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)+F_ji
      end do
    end subroutine doDivVectorMat79_3D

  end subroutine gfsys_buildVecScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildVecTVDBlock(RcoeffMatrices, rafcstab, rx,&
      fcb_calcFlux, fcb_calcCharacteristics, dscale, bclear, ry)

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
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim


    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildVecTVDScalar(RcoeffMatrices, rafcstab,&
          rx%RvectorBlock(1), fcb_calcFlux, fcb_calcCharacteristics,&
          dscale, bclear, ry%RvectorBlock(1))
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecTVDBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)


    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_Dpp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_Dpm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_Dqp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_Dqm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_Drp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_Drm)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecTVDBlock')
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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecTVDBlock')
      call sys_halt()
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,Dx_i,Dx_j
      integer :: iedge,ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's for X-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx_i, Dx_j, C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble high-order divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)+F_ji

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, XDir1D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for X-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, XDir1D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar = 1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)-F_ij
      end do
    end subroutine doLimitTVDMat79_1D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,Dx_i,Dx_j
      integer :: iedge,ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's for X-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx_i, Dx_j, C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble high-order divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)+F_ji

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, XDir2D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for X-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Clear P's and Q's for Y-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, XDir2D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij  =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij  =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij   =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar = 1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)-F_ij

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, YDir2D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for Y-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, YDir2D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar =1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)-F_ij
      end do

    end subroutine doLimitTVDMat79_2D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,Dx_i,Dx_j
      integer :: iedge,ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's for X-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx_i, Dx_j, C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble high-order divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)+F_ji

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, XDir3D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for X-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Clear P's and Q's for Y-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, XDir3D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij  =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij  =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij   =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar = 1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)-F_ij

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, YDir3D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for Y-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Clear P's and Q's for Z-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, YDir3D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar =1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)-F_ij

        ! Compute characteristic fluxes in Z-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, ZDir3D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffZ(ji)-DcoeffZ(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffZ(ij)+DcoeffZ(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for Z-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Get solution values at nodes
        Dx_i = Dx(i,:); Dx_j = Dx(j,:)

        ! Compute characteristic fluxes in Z-direction
        call fcb_calcCharacteristics(Dx_i, Dx_j, ZDir3D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffZ(ji)-DcoeffZ(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffZ(ji)-DcoeffZ(ij))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar =1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(i,:) = Dy(i,:)+F_ij
        Dy(j,:) = Dy(j,:)-F_ij
      end do

    end subroutine doLimitTVDMat79_3D

  end subroutine gfsys_buildVecTVDBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildVecTVDScalar(RcoeffMatrices, rafcstab, rx,&
      fcb_calcFlux, fcb_calcCharacteristics, dscale, bclear, ry)

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
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim


    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecTVDScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_Dpp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_Dpm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_Dqp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_Dqm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_Drp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_Drm)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecTVDScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecTVDScalar')
      call sys_halt()
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,Dx_i,Dx_j
      integer :: iedge,ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's for X-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble high-order divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)+F_ji

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), XDir1D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for X-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), XDir1D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij  =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij  =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij   =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar = 1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)-F_ij
      end do

    end subroutine doLimitTVDMat79_1D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,Dx_i,Dx_j
      integer :: iedge,ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble high-order divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)+F_ji

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), XDir2D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for X-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Clear P's and Q's for Y-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), XDir2D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij  =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij  =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij   =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar = 1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)-F_ij

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), YDir2D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for Y-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), YDir2D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar =1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)-F_ij
      end do

    end subroutine doLimitTVDMat79_2D


    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doLimitTVDMat79_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR*NVAR) :: R_ij,L_ij
      real(DP), dimension(NVAR) :: F_ij,F_ji,W_ij,Lbd_ij,ka_ij,ks_ij,Dx_i,Dx_j
      integer :: iedge,ij,ji,i,j,iloc,jloc,ivar


      ! Clear P's and Q's for X-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the fluxes
        call fcb_calcFlux(Dx(:,i), Dx(:,j), C_ij, C_ji,&
                          i, j, dscale, F_ij, F_ji)

        ! Assemble high-order divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)+F_ji

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), XDir3D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do
      end do


      ! Compute nodal correction factors for X-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Clear P's and Q's for Y-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute characteristic fluxes in X-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), XDir3D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij  =  0.5_DP*(DcoeffX(ji)-DcoeffX(ij))*Lbd_ij
        ks_ij  =  0.5_DP*(DcoeffX(ij)+DcoeffX(ji))*Lbd_ij
        F_ij   = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij   =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar = 1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)-F_ij

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), YDir3D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do
      end do


      ! Compute nodal correction factors for Y-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Clear P's and Q's for Z-direction
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute characteristic fluxes in Y-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), YDir3D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffY(ji)-DcoeffY(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffY(ij)+DcoeffY(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar =1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)-F_ij

        ! Compute characteristic fluxes in Z-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), ZDir3D, W_ij, Lbd_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffZ(ji)-DcoeffZ(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffZ(ij)+DcoeffZ(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij

        ! Update sums of downstream/upstream edge contributions
        do ivar = 1, NVAR

          ! Set node orientation
          if (ka_ij(ivar) > 0) then
            iloc = j; jloc = i
            F_ij(ivar) = -F_ij(ivar)
          else
            iloc = i; jloc = j
          end if

          ! Assemble P's and Q's
          if (F_ij(ivar) > 0) then
            Dpp(ivar,iloc) = Dpp(ivar,iloc)+F_ij(ivar)
            Dqm(ivar,iloc) = Dqm(ivar,iloc)-F_ij(ivar)
            Dqp(ivar,jloc) = Dqp(ivar,jloc)+F_ij(ivar)
          else
            Dpm(ivar,iloc) = Dpm(ivar,iloc)+F_ij(ivar)
            Dqp(ivar,iloc) = Dqp(ivar,iloc)-F_ij(ivar)
            Dqm(ivar,jloc) = Dqm(ivar,jloc)+F_ij(ivar)
          end if
        end do

      end do


      ! Compute nodal correction factors for Z-direction
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute characteristic fluxes in Z-direction
        call fcb_calcCharacteristics(Dx(:,i), Dx(:,j), ZDir3D, W_ij, Lbd_ij, R_ij)

        ! Compute antidiffusive fluxes
        ka_ij =  0.5_DP*(DcoeffZ(ji)-DcoeffZ(ij))*Lbd_ij
        ks_ij =  0.5_DP*(DcoeffZ(ij)+DcoeffZ(ji))*Lbd_ij
        F_ij  = -max(0.0_DP, min(abs(ka_ij)-ks_ij, 2.0_DP*abs(ka_ij)))*W_ij
        W_ij  =  abs(ka_ij)*W_ij

        ! Construct characteristic fluxes
        do ivar =1, NVAR

          ! Limit characteristic fluxes
          if (ka_ij(ivar) < 0) then
            if (F_ij(ivar) > 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,i)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,i)*F_ij(ivar)
            end if
          else
            if (F_ij(ivar) < 0) then
              W_ij(ivar) = W_ij(ivar)+Drp(ivar,j)*F_ij(ivar)
            else
              W_ij(ivar) = W_ij(ivar)+Drm(ivar,j)*F_ij(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        call DGEMV('n', NVAR, NVAR, dscale, R_ij,&
                   NVAR, W_ij, 1, 0.0_DP, F_ij, 1)

        ! Assemble high-resolution divergence vector
        Dy(:,i) = Dy(:,i)+F_ij
        Dy(:,j) = Dy(:,j)-F_ij
      end do

    end subroutine doLimitTVDMat79_3D

  end subroutine gfsys_buildVecTVDScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildVecFCTBlock(rlumpedMassMatrix,&
      rafcstab, rx, dscale, bclear, ioperationSpec, ry,&
      fcb_calcFluxTransformation, fcb_calcDiffTransformation)

!<description>
    ! This subroutine assembles the divergence vector for nonlinear
    ! FEM-FCT schemes.  If the vectors contain only one block, then
    ! the scalar counterpart of this routine is called with the scalar
    ! subvectors. Consider the documentation of subroutine
    ! 'gfsys_buildVecFCTScalar' for further details.
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

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation
    optional :: fcb_calcDiffTransformation
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry
    !</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1)) then
      call gfsys_buildVecFCTScalar(rlumpedMassMatrix,&
          rafcstab, rx%RvectorBlock(1), dscale, bclear,&
          ioperationSpec, ry%RvectorBlock(1),&
          fcb_calcFluxTransformation, fcb_calcDiffTransformation)
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBLock')
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

    if (iand(ioperationSpec, AFCSTAB_FCTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if

    if (iand(ioperationSpec, AFCSTAB_FCTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute sums of antidiffusive increments
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_Dpp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_Dpm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux)
      else
        call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
      end if

      ! Compute sums of antidiffusive increments
      if (present(fcb_calcFluxTransformation) .and.&
          present(fcb_calcDiffTransformation)) then
        call doADIncrementsTransformed(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            rafcstab%NVARtransformed, p_Dx,&
            p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
      else
        call doADIncrements(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
      end if

      ! Set specifiers
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute local bounds
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_Dqp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_Dqm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute local bounds
      if (present(fcb_calcDiffTransformation)) then
        call doBoundsTransformed(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            rafcstab%NVARtransformed, p_Dx, p_Dqp, p_Dqm)
      else
        call doBounds(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            p_Dx, p_Dqp, p_Dqm)
      end if

      ! Set specifiers
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)       .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_Dpp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_Dpm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_Dqp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_Dqm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_Drp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_Drm)

      ! Compute nodal correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call doLimitNodal(rafcstab%NEQ, rafcstab%NVARtransformed,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        call doLimitNodalConstrained(rafcstab%NEQ, rafcstab%NVARtransformed,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      end if

      ! Set specifier
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides nodal correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_Drp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_Drm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux0)

        if (present(fcb_calcFluxTransformation)) then
          call doLimitEdgewiseConstrainedTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, rafcstab%NVARtransformed, p_Dx,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          call doLimitEdgewiseConstrained(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      else

        if (present(fcb_calcFluxTransformation)) then
          call doLimitEdgewiseTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, rafcstab%NVARtransformed, p_Dx,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          call doLimitEdgewise(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      end if

      ! Set specifier
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides edgewise correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTBlock')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(ry, p_Dy)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Clear divergence vector?
      if (bclear) call lsysbl_clearVector(ry)

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then
        call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
        call doCorrectScaleByMass(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
      else
        call doCorrect(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            dscale, p_Dalpha, p_Dflux, p_Dy)
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation

    subroutine doADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dx, Dflux, Dalpha, Dpp, Dpm)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(out) :: Dpp,Dpm

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

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

            real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVARtransformed) :: F_ij,F_ji,Diff
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
        Diff = Dalpha(iedge) * Dflux(:,iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(i,:), Dx(j,:), Diff, F_ij, F_ji)

        ! Compute transformed solution difference
        call fcb_calcDiffTransformation(&
            Dx(i,:), Dx(j,:), Diff)

        ! Compute the sums of positive/negative antidiffusive increments
        Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
        Dpp(:,j) = Dpp(:,j)+max(0.0_DP, F_ji)
        Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
        Dpm(:,j) = Dpm(:,j)+min(0.0_DP, F_ji)
      end do
    end subroutine doADIncrementsTransformed

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! without transformation

    subroutine doBounds(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dx, Dqp, Dqm)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(out) :: Dqp,Dqm

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

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dqp,Dqm

      ! local variables
      real(DP), dimension(NVARtransformed) :: Diff
      integer :: iedge,i,j

      ! Clear Q`s
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute transformed solution difference
        call fcb_calcDiffTransformation(&
            Dx(i,:), Dx(j,:), Diff)

        ! Compute the distance to a local extremum
        ! of the predicted solution
        Dqp(:,i) = max(Dqp(:,i), Diff)
        Dqp(:,j) = max(Dqp(:,j),-Diff)
        Dqm(:,i) = min(Dqm(:,i), Diff)
        Dqm(:,j) = min(Dqm(:,j),-Diff)
      end do
    end subroutine doBoundsTransformed

    !**************************************************************
    ! Compute nodal correction factors without constraints

    subroutine doLimitNodal(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drp(:,ieq) = ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq)+SYS_EPSREAL)
      end do
      !$omp end parallel do

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drm(:,ieq) = ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq)-SYS_EPSREAL)
      end do
      !$omp end parallel do
    end subroutine doLimitNodal

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrained(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drp(:,ieq) = min(1.0_DP,&
            ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq)+SYS_EPSREAL))
      end do
      !$omp end parallel do

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drm(:,ieq) = min(1.0_DP,&
            ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq)-SYS_EPSREAL))
      end do
      !$omp end parallel do
    end subroutine doLimitNodalConstrained

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewise(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Drp, Drm, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVARtransformed) :: F_ij,F_ji,R_ij,R_ji
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(i,:), Dx(j,:), Dflux(:,iedge), F_ij, F_ji)

        ! Compute nodal correction factors
        R_ij = merge(Drp(:,i), Drm(:,i), F_ij .ge. 0.0_DP)
        R_ji = merge(Drp(:,j), Drm(:,j), F_ji .ge. 0.0_DP)

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
      end do
    end subroutine doLimitEdgewiseTransformed

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrained(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux1, Dflux2, Drp, Drm, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVARtransformed) :: F1_ij,F1_ji,F2_ij,F2_ji,R_ij,R_ji
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(i,:), Dx(j,:), Dflux1(:,iedge), F1_ij, F1_ji)
        call fcb_calcFluxTransformation(&
            Dx(i,:), Dx(j,:), Dflux2(:,iedge), F2_ij, F2_ji)

        ! Compute nodal correction factors
        where (F1_ij*F2_ij .le. 0.0_DP)
          R_ij = 0.0_DP
        elsewhere
          R_ij = min(1.0_DP,&
              F1_ij/F2_ij*merge(Drp(:,i), Drm(:,i), F1_ij .ge. 0.0_DP))
        end where

        where (F1_ji*F2_ji .le. 0.0_DP)
          R_ji = 0.0_DP
        elsewhere
          R_ji = min(1.0_DP,&
              F1_ji/F2_ji*merge(Drp(:,j), Drm(:,j), F1_ji .ge. 0.0_DP))
        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
      end do
    end subroutine doLimitEdgewiseConstrainedTransformed

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrect(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, dscale, Dalpha, Dflux, Dy)

      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

      real(DP), dimension(:), intent(in) :: Dalpha,ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

  end subroutine gfsys_buildVecFCTBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildVecFCTScalar(rlumpedMassMatrix,&
      rafcstab, rx, dscale, bclear, ioperationSpec, ry,&
      fcb_calcFluxTransformation, fcb_calcDiffTransformation)

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

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation
    optional :: fcb_calcDiffTransformation
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
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


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute sums of antidiffusive increments
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_Dpp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_Dpm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux)
      else
        call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
      end if

      ! Compute sums of antidiffusive increments
      if (present(fcb_calcFluxTransformation) .and.&
          present(fcb_calcDiffTransformation)) then
        if (rafcstab%bprelimiting) then
          call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux0)
          call doPreADIncrementsTransformed(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              rafcstab%NVARtransformed, p_Dx,&
              p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
        else
          call doADIncrementsTransformed(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              rafcstab%NVARtransformed, p_Dx,&
              p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        end if
      else
        if (rafcstab%bprelimiting) then
          call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux0)
          call doPreADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dx, p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
        else
          call doADIncrements(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        end if
      end if

      ! Set specifiers
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute local bounds
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_Dqp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_Dqm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute local bounds
      if (present(fcb_calcDiffTransformation)) then
        call doBoundsTransformed(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            rafcstab%NVARtransformed, p_Dx, p_Dqp, p_Dqm)
      else
        call doBounds(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            p_Dx, p_Dqp, p_Dqm)
      end if

      ! Set specifiers
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)       .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_Dpp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_Dpm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_Dqp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_Dqm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_Drp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_Drm)

      ! Compute nodal correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call doLimitNodal(rafcstab%NEQ, rafcstab%NVARtransformed,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        call doLimitNodalConstrained(rafcstab%NEQ, rafcstab%NVARtransformed,&
            dscale, p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      end if

      ! Set specifier
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides nodal correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_Drp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_Drm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux0)

        if (present(fcb_calcFluxTransformation)) then
          call doLimitEdgewiseConstrainedTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, rafcstab%NVARtransformed, p_Dx,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          call doLimitEdgewiseConstrained(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      else

        if (present(fcb_calcFluxTransformation)) then
          call doLimitEdgewiseTransformed(&
              p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, rafcstab%NVARtransformed, p_Dx,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          call doLimitEdgewise(p_IverticesAtEdge,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      end if

      ! Set specifier
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides edgewise correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(ry, p_Dy)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Clear divergence vector?
      if (bclear) call lsyssc_clearVector(ry)

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then
        call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
        call doCorrectScaleByMass(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
      else
        call doCorrect(p_IverticesAtEdge,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            dscale, p_Dalpha, p_Dflux, p_Dy)
      end if
    end if
    
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dx, Dflux, Dalpha, Dpp, Dpm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(out) :: Dpp,Dpm

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

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: Diff
      real(DP), dimension(NVARtransformed) :: F_ij,F_ji
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
        Diff = Dalpha(iedge)*Dflux(:,iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(:,i), Dx(:,j), Diff, F_ij, F_ji)
        
        ! Compute the sums of positive/negative antidiffusive increments
        Dpp(:,i) = Dpp(:,i)+max(0.0_DP, F_ij)
        Dpp(:,j) = Dpp(:,j)+max(0.0_DP, F_ji)
        Dpm(:,i) = Dpm(:,i)+min(0.0_DP, F_ij)
        Dpm(:,j) = Dpm(:,j)+min(0.0_DP, F_ji)
      end do
    end subroutine doADIncrementsTransformed

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and with prelimiting

    subroutine doPreADIncrements(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dx, Dflux, Dflux0, Dalpha, Dpp, Dpm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(NVAR,NEQ), intent(out) :: Dpp,Dpm

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

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: Diff
      real(DP), dimension(NVARtransformed) :: F_ij,F_ji,G_ij,G_ji
      real(DP) :: alpha_ij,alpha_ji
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
        Diff = Dalpha(iedge)*Dflux(:,iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(:,i), Dx(:,j), Diff, F_ij, F_ji)
        
        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(:,i), Dx(:,j), Dflux0(:,iedge), G_ij, G_ji)
        
        ! MinMod prelimiting
        alpha_ij = minval(mprim_minmod3(F_ij, G_ij, F_ij))
        alpha_ji = minval(mprim_minmod3(F_ji, G_ji, F_ji))

        ! Synchronisation of correction factors
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
    end subroutine doPreADIncrementsTransformed
    
    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! without transformation

    subroutine doBounds(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dx, Dqp, Dqm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(out) :: Dqp,Dqm

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

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dqp,Dqm

      ! local variables
      real(DP), dimension(NVARtransformed) :: Diff
      integer :: iedge,i,j

      ! Clear Q`s
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute transformed solution difference
        call fcb_calcDiffTransformation(&
            Dx(:,i), Dx(:,j), Diff)

        ! Compute the distance to a local extremum
        ! of the predicted solution
        Dqp(:,i) = max(Dqp(:,i), Diff)
        Dqp(:,j) = max(Dqp(:,j),-Diff)
        Dqm(:,i) = min(Dqm(:,i), Diff)
        Dqm(:,j) = min(Dqm(:,j),-Diff)
      end do
    end subroutine doBoundsTransformed

    !**************************************************************
    ! Compute nodal correction factors without constraints

    subroutine doLimitNodal(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drp(:,ieq) = ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq)+SYS_EPSREAL)
      end do
      !$omp end parallel do

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drm(:,ieq) = ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq)-SYS_EPSREAL)
      end do
      !$omp end parallel do
    end subroutine doLimitNodal

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrained(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      integer :: ieq

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drp(:,ieq) = min(1.0_DP,&
            ML(ieq)*Dqp(:,ieq)/(dscale*Dpp(:,ieq)+SYS_EPSREAL))
      end do
      !$omp end parallel do

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drm(:,ieq) = min(1.0_DP,&
            ML(ieq)*Dqm(:,ieq)/(dscale*Dpm(:,ieq)-SYS_EPSREAL))
      end do
      !$omp end parallel do

    end subroutine doLimitNodalConstrained

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewise(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux, Drp, Drm, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVARtransformed) :: F_ij,F_ji,R_ij,R_ji
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(:,i), Dx(:,j), Dflux(:,iedge), F_ij, F_ji)

        ! Compute nodal correction factors
        R_ij = merge(Drp(:,i), Drm(:,i), F_ij .ge. 0.0_DP)
        R_ji = merge(Drp(:,j), Drm(:,j), F_ji .ge. 0.0_DP)

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
      end do
    end subroutine doLimitEdgewiseTransformed

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrained(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, Dflux1, Dflux2, Drp, Drm, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVARtransformed) :: F1_ij,F1_ji,F2_ij,F2_ji,R_ij,R_ji
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute transformed fluxes
        call fcb_calcFluxTransformation(&
            Dx(:,i), Dx(:,j), Dflux1(:,iedge), F1_ij, F1_ji)
        call fcb_calcFluxTransformation(&
            Dx(:,i), Dx(:,j), Dflux2(:,iedge), F2_ij, F2_ji)

        ! Compute nodal correction factors
        where (F1_ij*F2_ij .le. 0.0_DP)
          R_ij = 0.0_DP
        elsewhere
          R_ij = min(1.0_DP,&
              F1_ij/F2_ij*merge(Drp(:,i), Drm(:,i), F1_ij .ge. 0.0_DP))
        end where

        where (F1_ji*F2_ji .le. 0.0_DP)
          R_ji = 0.0_DP
        elsewhere
          R_ji = min(1.0_DP,&
              F1_ji/F2_ji*merge(Drp(:,j), Drm(:,j), F1_ji .ge. 0.0_DP))
        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
      end do
    end subroutine doLimitEdgewiseConstrainedTransformed

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrect(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, dscale, Dalpha, Dflux, Dy)

      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

      real(DP), dimension(:), intent(in) :: Dalpha,ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

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

  end subroutine gfsys_buildVecFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildFluxFCTBlock(rlumpedMassMatrix,&
      RcoeffMatrices, rafcstab, rx1, rx2, fcb_calcFluxFCT,&
      theta, tstep, dscale, binit, rconsistentMassMatrix)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for FEM-FCT schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vectors
    type(t_vectorBlock), intent(in) :: rx1, rx2

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

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rconsistentMassMatrix

    ! callback functions to compute local fluxes
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_MC,p_DcoeffX,p_DcoeffY,p_DcoeffZ
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dflux0,p_Dflux,p_Dalpha
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim

    ! Check if block vector contains only one block
    if ((rx1%nblocks .eq. 1) .and. (rx2%nblocks .eq. 1)) then
      call gfsys_buildFluxFCTScalar(rlumpedMassMatrix,&
          RcoeffMatrices, rafcstab, rx1%RvectorBlock(1),&
          rx2%RvectorBlock(1), fcb_calcFluxFCT,&
          theta, tstep, dscale, binit, rconsistentMassMatrix)
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Set pointers
    call lsysbl_getbase_double(rx1, p_Dx1)
    call lsysbl_getbase_double(rx2, p_Dx2)
    call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end select


    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

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

      ! The current flux is stored at position 2
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)

      ! The initial flux (including the iterative correction for
      ! the iterative FEM-FCT algorithm) is stored at position 3
      call lsyssc_getbase_double(rafcstab%RedgeVectors(3), p_Dflux0)

      ! Check if the amount of rejected antidiffusion should be
      ! included in the initial raw antidiffusive fluxes
      if (.not.binit .and.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_ITERATIVE) then

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
          call output_line('Stabilisation does not provide correction factors',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTBlock')
          call sys_halt()
        end if

        ! Set pointer
        call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)

        ! Subtract amount of rejected antidiffusion
        call doCombineFluxes(rafcstab%NVAR, rafcstab%NEDGE,&
            -1.0_DP, p_Dalpha, p_Dflux, p_Dflux0)
      end if

      ! Do we have to use the consistent mass matrix?
      if (present(rconsistentMassMatrix)) then

        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointer for consistent mass matrix
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            if (binit) then
              call doFluxesMat79ConsMass_1D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_Dx1, p_Dx2,&
                  -dscale/tstep, (1-theta)*dscale, p_Dflux0)
              call doFluxesMat79NoMass_1D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_DcoeffX, p_Dx2, dscale, p_Dflux)
            else
              call doFluxesMat79ConsMass_1D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_Dx1, p_Dx2,&
                  dscale/tstep, theta*dscale, p_Dflux)
              call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
            end if

          case (NDIM2D)
            if (binit) then
              call doFluxesMat79ConsMass_2D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_Dx1, p_Dx2,&
                  -dscale/tstep, (1-theta)*dscale, p_Dflux0)
              call doFluxesMat79NoMass_2D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_DcoeffX, p_DcoeffY, p_Dx2, dscale, p_Dflux)
            else
              call doFluxesMat79ConsMass_2D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_Dx1, p_Dx2,&
                  dscale/tstep, theta*dscale, p_Dflux)
              call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
            end if

          case (NDIM3D)
            if (binit) then
              call doFluxesMat79ConsMass_3D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx1, p_Dx2,&
                  -dscale/tstep, (1-theta)*dscale, p_Dflux0)
              call doFluxesMat79NoMass_3D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx2, dscale, p_Dflux)
            else
              call doFluxesMat79ConsMass_3D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx1, p_Dx2,&
                  dscale/tstep, theta*dscale, p_Dflux)
              call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
            end if
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
          call sys_halt()
        end select

      else

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call doFluxesMat79NoMass_1D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_Dx2, dscale, p_Dflux)

          case (NDIM2D)
            call doFluxesMat79NoMass_2D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_Dx2, dscale, p_Dflux)

          case (NDIM3D)
            call doFluxesMat79NoMass_3D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx2, dscale, p_Dflux)
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
          call sys_halt()
        end select

        ! Combine explicit and implicit fluxes
        if (binit) then
          call lalg_copyVector(p_Dflux, p_Dflux0)
        elseif (1-theta .gt. SYS_EPSREAL) then
          call lalg_vectorLinearComb(&
              p_Dflux0, p_Dflux, 1-theta, theta)
        elseif (theta .gt. SYS_EPSREAL) then
          call lalg_scaleVector(p_Dflux, theta)
        else
          call lalg_clearVector(p_Dflux)
        end if

      end if

      ! Do we have to store the initial fluxes separately?
      if (binit .and. (rafcstab%ctypeAFCstabilisation &
                       .eq. AFCSTAB_FEMFCT_IMPLICIT)) then
        ! The initial flux without contribution of the
        ! consistent mass matrix is stored at position 4
        call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
        call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux0)
        call lalg_copyVector(p_Dflux, p_Dflux0)
      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_FEMFCT_LINEARISED)

      !-------------------------------------------------------------------------
      ! Linearised FEM-FCT algorithm
      !-------------------------------------------------------------------------

      ! The linearised flux is stored at position 2
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)

      ! Do we have to use the consistent mass matrix?
      if (present(rconsistentMassMatrix)) then

        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointer for consistent mass matrix
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call doFluxesMat79ConsMass_1D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_MC, p_DcoeffX, p_Dx1, p_Dx2, dscale, dscale, p_Dflux)

          case (NDIM2D)
            call doFluxesMat79ConsMass_2D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_MC, p_DcoeffX, p_DcoeffY, p_Dx1, p_Dx2,&
                dscale, dscale, p_Dflux)

          case (NDIM3D)
            call doFluxesMat79ConsMass_3D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_MC, p_DcoeffX, p_DcoeffY, p_DcoeffZ,&
                p_Dx1, p_Dx2, dscale, dscale, p_Dflux)
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
          call sys_halt()
        end select

      else

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call doFluxesMat79NoMass_1D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_Dx2, dscale, p_Dflux)

          case (NDIM2D)
            call doFluxesMat79NoMass_2D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_Dx2, dscale, p_Dflux)

          case (NDIM3D)
            call doFluxesMat79NoMass_3D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx2, dscale, p_Dflux)
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
          call sys_halt()
        end select

      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)


    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Combine two fluxes: flux2 := flux2+dscale*alpha*flux2

    subroutine doCombineFluxes(NVAR, NEDGE, dscale, Dalpha, Dflux1, Dflux2)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR, NEDGE

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux2

      ! local variables
      integer :: iedge

      !$omp parallel do
      do iedge = 1, NEDGE
        Dflux2(:,iedge) = Dflux2(:,iedge) +&
            dscale * Dalpha(iedge) * Dflux1(:,iedge)
      end do
      !$omp end parallel do

    end subroutine doCombineFluxes

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with contribution
    ! of the consistent mass matrix in 1D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79ConsMass_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, MC, DcoeffX, Dx1, Dx2, dscale1, dscale2, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx1,Dx2
      real(DP), dimension(:), intent(in) :: MC,DcoeffX
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx1(i,:), Dx1(j,:), Dx2(i,:), Dx2(j,:), C_ij, C_ji,&
            i, j, dscale1*MC(ij), dscale2, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79ConsMass_1D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with contribution
    ! of the consistent mass matrix in 2D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79ConsMass_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, MC, DcoeffX, DcoeffY, Dx1, Dx2,&
        dscale1, dscale2, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx1,Dx2
      real(DP), dimension(:), intent(in) :: MC,DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx1(i,:), Dx1(j,:), Dx2(i,:), Dx2(j,:), C_ij, C_ji,&
            i, j, dscale1*MC(ij), dscale2, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79ConsMass_2D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with contribution
    ! of the consistent mass matrix in 3D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79ConsMass_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, MC, DcoeffX, DcoeffY, DcoeffZ, Dx1, Dx2,&
        dscale1, dscale2, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx1,Dx2
      real(DP), dimension(:), intent(in) :: MC,DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx1(i,:), Dx1(j,:), Dx2(i,:), Dx2(j,:), C_ij, C_ji,&
            i, j, dscale1*MC(ij), dscale2, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79ConsMass_3D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without contribution
    ! of the consistent mass matrix in 1D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79NoMass_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx(i,:), Dx(j,:), Dx(i,:), Dx(j,:), C_ij, C_ji,&
            i, j, 0.0_DP, dscale, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79NoMass_1D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without contribution
    ! of the consistent mass matrix in 2D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79NoMass_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx(i,:), Dx(j,:), Dx(i,:), Dx(j,:), C_ij, C_ji,&
            i, j, 0.0_DP, dscale, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79NoMass_2D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without contribution
    ! of the consistent mass matrix in 3D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79NoMass_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx(i,:), Dx(j,:), Dx(i,:), Dx(j,:), C_ij, C_ji,&
            i, j, 0.0_DP, dscale, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79NoMass_3D

  end subroutine gfsys_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildFluxFCTScalar(rlumpedMassMatrix,&
      RcoeffMatrices, rafcstab, rx1, rx2, fcb_calcFluxFCT,&
      theta, tstep, dscale, binit, rconsistentMassMatrix)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! FEM-FCT schemes. Note that the vectors are required as scalar
    ! vectors which are stored in the interleave format.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx1, rx2

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

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rconsistentMassMatrix

    ! callback functions to compute local fluxes
    include 'intf_gfsyscallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_MC,p_DcoeffX,p_DcoeffY,p_DcoeffZ
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dflux0,p_Dflux,p_Dalpha
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: ndim

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Set pointers
    call lsyssc_getbase_double(rx1, p_Dx1)
    call lsyssc_getbase_double(rx2, p_Dx2)
    call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

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
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end select


    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)

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

      ! The current flux is stored at position 2
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)

      ! The initial flux (including the iterative correction for
      ! the iterative FEM-FCT algorithm) is stored at position 3
      call lsyssc_getbase_double(rafcstab%RedgeVectors(3), p_Dflux0)

      ! Check if the amount of rejected antidiffusion should be
      ! included in the initial raw antidiffusive fluxes
      if (.not.binit .and.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_ITERATIVE) then

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
          call output_line('Stabilisation does not provide correction factors',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        ! Set pointer
        call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_Dalpha)

        ! Subtract amount of rejected antidiffusion
        call doCombineFluxes(rafcstab%NVAR, rafcstab%NEDGE,&
            -1.0_DP, p_Dalpha, p_Dflux, p_Dflux0)
      end if

      ! Do we have to use the consistent mass matrix?
      if (present(rconsistentMassMatrix)) then

        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointer for consistent mass matrix
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            if (binit) then
              call doFluxesMat79ConsMass_1D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_Dx1, p_Dx2,&
                  -dscale/tstep, (1-theta)*dscale, p_Dflux0)
              call doFluxesMat79NoMass_1D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_DcoeffX, p_Dx2, dscale, p_Dflux)
            else
              call doFluxesMat79ConsMass_1D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_Dx1, p_Dx2,&
                  dscale/tstep, theta*dscale, p_Dflux)
              call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
            end if

          case (NDIM2D)
            if (binit) then
              call doFluxesMat79ConsMass_2D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_Dx1, p_Dx2,&
                  -dscale/tstep, (1-theta)*dscale, p_Dflux0)
              call doFluxesMat79NoMass_2D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_DcoeffX, p_DcoeffY, p_Dx2, dscale, p_Dflux)
            else
              call doFluxesMat79ConsMass_2D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_Dx1, p_Dx2,&
                  dscale/tstep, theta*dscale, p_Dflux)
              call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
            end if

          case (NDIM3D)
            if (binit) then
              call doFluxesMat79ConsMass_3D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx1, p_Dx2,&
                  -dscale/tstep, (1-theta)*dscale, p_Dflux0)
              call doFluxesMat79NoMass_3D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx2, dscale, p_Dflux)
            else
              call doFluxesMat79ConsMass_3D(p_IverticesAtEdge,&
                  rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                  p_MC, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx1, p_Dx2,&
                  dscale/tstep, theta*dscale, p_Dflux)
              call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
            end if
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
          call sys_halt()
        end select

      else

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call doFluxesMat79NoMass_1D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_Dx2, dscale, p_Dflux)

          case (NDIM2D)
            call doFluxesMat79NoMass_2D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_Dx2, dscale, p_Dflux)

          case (NDIM3D)
            call doFluxesMat79NoMass_3D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx2, dscale, p_Dflux)
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
          call sys_halt()
        end select

        ! Combine explicit and implicit fluxes
        if (binit) then
          call lalg_copyVector(p_Dflux, p_Dflux0)
        elseif (1-theta .gt. SYS_EPSREAL) then
          call lalg_vectorLinearComb(&
              p_Dflux0, p_Dflux, 1-theta, theta)
        elseif (theta .gt. SYS_EPSREAL) then
          call lalg_scaleVector(p_Dflux, theta)
        else
          call lalg_clearVector(p_Dflux)
        end if

      end if

      ! Do we have to store the initial fluxes separately?
      if (binit .and. (rafcstab%ctypeAFCstabilisation &
                       .eq. AFCSTAB_FEMFCT_IMPLICIT)) then
        ! The initial flux without contribution of the
        ! consistent mass matrix is stored at position 4
        call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
        call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux0)
        call lalg_copyVector(p_Dflux, p_Dflux0)
      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_FEMFCT_LINEARISED)

      !-------------------------------------------------------------------------
      ! Linearised FEM-FCT algorithm
      !-------------------------------------------------------------------------

      ! Do we have to use the consistent mass matrix?
      if (present(rconsistentMassMatrix)) then
        
        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! The linearised flux is stored at position 2
        call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)

        ! Set pointer for consistent mass matrix
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)

        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call doFluxesMat79ConsMass_1D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_MC, p_DcoeffX, p_Dx1, p_Dx2, dscale, dscale, p_Dflux)

          case (NDIM2D)
            call doFluxesMat79ConsMass_2D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_MC, p_DcoeffX, p_DcoeffY, p_Dx1, p_Dx2,&
                dscale, dscale, p_Dflux)

          case (NDIM3D)
            call doFluxesMat79ConsMass_3D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_MC, p_DcoeffX, p_DcoeffY, p_DcoeffZ,&
                p_Dx1, p_Dx2, dscale, dscale, p_Dflux)
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
          call sys_halt()
        end select

      end if

      if (not(present(rconsistentMassMatrix)) .or.&
          rafcstab%bprelimiting) then

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        if (rafcstab%bprelimiting) then
          ! The flux for prelimiting is stored at position 4
          call lsyssc_getbase_double(rafcstab%RedgeVectors(4), p_Dflux)
        else
          ! The linearised flux is stored at position 2
          call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_Dflux)
        end if
          
        ! What kind of matrix are we?
        select case(RcoeffMatrices(1)%cmatrixFormat)
        case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
          !---------------------------------------------------------------------
          ! Matrix format 7 and 9
          !---------------------------------------------------------------------

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call doFluxesMat79NoMass_1D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_Dx2, dscale, p_Dflux)

          case (NDIM2D)
            call doFluxesMat79NoMass_2D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_Dx2, dscale, p_Dflux)

          case (NDIM3D)
            call doFluxesMat79NoMass_3D(p_IverticesAtEdge,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx2, dscale, p_Dflux)
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
          call sys_halt()
        end select
        
        ! Prelimiting is only necessary if the consistent mass matrix
        ! is built into the raw antidiffusive fluxes. However, if the
        ! switch for prelimiting was not set to .false. and no mass
        ! antidiffusion is built into the fluxes we can simply copy it
        if (not(present(rconsistentMassMatrix)) .and.&
            rafcstab%bprelimiting) then
          call lsyssc_copyVector(rafcstab%RedgeVectors(4),&
              rafcstab%RedgeVectors(2))
        end if
        
      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)


    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Combine two fluxes: flux2 := flux2+dscale*alpha*flux2

    subroutine doCombineFluxes(NVAR, NEDGE, dscale, Dalpha, Dflux1, Dflux2)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR, NEDGE

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux2

      ! local variables
      integer :: iedge

      !$omp parallel do
      do iedge = 1, NEDGE
        Dflux2(:,iedge) = Dflux2(:,iedge) +&
            dscale * Dalpha(iedge) * Dflux1(:,iedge)
      end do
      !$omp end parallel do

    end subroutine doCombineFluxes

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with contribution
    ! of the consistent mass matrix in 1D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79ConsMass_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, MC, DcoeffX, Dx1, Dx2, dscale1, dscale2, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx1,Dx2
      real(DP), dimension(:), intent(in) :: MC,DcoeffX
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx1(:,i), Dx1(:,j), Dx2(:,i), Dx2(:,j), C_ij, C_ji,&
            i, j, dscale1*MC(ij), dscale2, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79ConsMass_1D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with contribution
    ! of the consistent mass matrix in 2D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79ConsMass_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, MC, DcoeffX, DcoeffY, Dx1, Dx2,&
        dscale1, dscale2, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx1,Dx2
      real(DP), dimension(:), intent(in) :: MC,DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx1(:,i), Dx1(:,j), Dx2(:,i), Dx2(:,j), C_ij, C_ji,&
            i, j, dscale1*MC(ij), dscale2, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79ConsMass_2D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with contribution
    ! of the consistent mass matrix in 3D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79ConsMass_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, MC, DcoeffX, DcoeffY, DcoeffZ,&
        Dx1, Dx2, dscale1, dscale2, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx1,Dx2
      real(DP), dimension(:), intent(in) :: MC,DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx1(:,i), Dx1(:,j), Dx2(:,i), Dx2(:,j), C_ij, C_ji,&
            i, j, dscale1*MC(ij), dscale2, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79ConsMass_3D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without contribution
    ! of the consistent mass matrix in 1D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79NoMass_1D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, Dx, dscale, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx(:,i), Dx(:,j), Dx(:,i), Dx(:,j), C_ij, C_ji,&
            i, j, 0.0_DP, dscale, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79NoMass_1D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without contribution
    ! of the consistent mass matrix in 2D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79NoMass_2D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, Dx, dscale, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx(:,i), Dx(:,j), Dx(:,i), Dx(:,j), C_ij, C_ji,&
            i, j, 0.0_DP, dscale, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79NoMass_2D

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without contribution
    ! of the consistent mass matrix in 3D.
    ! All matrices are stored in matrix format 7 and 9

    subroutine doFluxesMat79NoMass_3D(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, DcoeffX, DcoeffY, DcoeffZ, Dx, dscale, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,ij,ji,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute the raw antidiffusives fluxes
        call fcb_calcFluxFCT(&
            Dx(:,i), Dx(:,j), Dx(:,i), Dx(:,j), C_ij, C_ji,&
            i, j, 0.0_DP, dscale, Dflux(:,iedge))
      end do
    end subroutine doFluxesMat79NoMass_3D

  end subroutine gfsys_buildFluxFCTScalar

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<subroutine>

  subroutine gfsys_getbase_double(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the block matrix.
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
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_double')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),&
                                     rarray(iblock,jblock)%Da)
        else
          nullify(rarray(iblock,jblock)%Da)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull
  end subroutine gfsys_getbase_double

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_getbase_single(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the block matrix.
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
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_single')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          call lsyssc_getbase_single(rmatrix%RmatrixBlock(iblock,jblock),&
                                     rarray(iblock,jblock)%Fa)
        else
          nullify(rarray(iblock,jblock)%Fa)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull
  end subroutine gfsys_getbase_single

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_getbase_int(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the block matrix.
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
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsys_getbase_int')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          call lsyssc_getbase_int(rmatrix%RmatrixBlock(iblock,jblock),&
                                  rarray(iblock,jblock)%Ia)
        else
          nullify(rarray(iblock,jblock)%Da)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull
  end subroutine gfsys_getbase_int

end module groupfemsystem

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
!# to the book chapter
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
!# There are three types of routines. The gfsc_buildXXXOperator
!# routines can be used to assemble the discrete convection or
!# diffusion operators resulting from the standard Galerkin finite
!# element discretisation plus some discretely defined artificial
!# diffusion. For convective terms, this technique is termed discrete
!# upwinding, whereas for physical diffusion operators this approach
!# ensures that the discrete maximum principle holds. Consequently,
!# the term DMP is adopted.
!#
!# The second type of routines is given by
!# gfsc_buildConvVectorXXX. They can be used to update/initialise the
!# convectove term applying some sort of algebraic flux
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
!# 1.) gfsc_initStabilisation
!#     -> Initializes the stabilisation structure
!#
!# 2.) gfsc_isMatrixCompatible
!#     -> Checks wether a matrix and a stabilisation structure are compatible
!#
!# 3.) gfsc_isVectorCompatible
!#     -> Checks wether a vector and a stabilisation structure are compatible
!#
!# 4.) gfsc_buildConvectionOperator = gfsc_buildConvOperatorScalar /
!#                                    gfsc_buildConvOperatorBlock
!#     -> Assembles the convective part of the transport operator
!#
!# 5.) gfsc_buildDiffusionOperator = gfsc_buildDiffusionOperator
!#     -> Assembles the diffusive part of the transport operator
!#
!# 6.) gfsc_buildConvVectorFCT = gfsc_buildConvVecFCTScalar /
!#                               gfsc_buildConvVecFCTBlock
!#     -> Assembles the convective vector for AFC stabilisation of FCT type
!#
!# 7.) gfsc_buildConvVectorTVD = gfsc_buildConvVecTVDScalar /
!#                               gfsc_buildConvVecTVDBlock
!#     -> Assembles the convective vector for AFC stabilisation of TVD type
!#
!# 8.) gfsc_buildConvVectorGP = gfsc_buildConvVecGPScalar /
!#                              gfsc_buildConvVecGPBlock
!#     -> Assembles the convective vector for AFC stabilisation of general-purpose type
!#
!# 9.) gfsc_buildConvVectorSymm = gfsc_buildConvVecSymmScalar /
!#                                gfsc_buildConvVecSymmBlock
!#     -> Assembles the convective vector for stabilisation by means of
!#        symmetric flux limiting for diffusion operators
!#
!# 10.) gfsc_buildConvectionJacobian = gfsc_buildConvJacobianScalar /
!#                                     gfsc_buildConvJacobianBlock
!#     -> Assembles the Jacobian matrix for the convective part of
!#        the transport operator for a scalar convection equation
!#
!# 11.) gfsc_buildJacobianFCT = gfsc_buildJacLinearFCTScalar /
!#                              gfsc_buildJacLinearFCTBlock /
!#                              gfsc_buildJacobianFCTScalar /
!#                              gfsc_buildJacobianFCTBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of FCT type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 12.) gfsc_buildJacobianTVD = gfsc_buildJacLinearTVDScalar /
!#                              gfsc_buildJacLinearTVDBlock /
!#                              gfsc_buildJacobianTVDScalar /
!#                              gfsc_buildJacobianTVDBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of TVD type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 13.) gfsc_buildJacobianGP = gfsc_buildJacLinearGPScalar /
!#                             gfsc_buildJacLinearGPBlock /
!#                             gfsc_buildJacobianGPScalar /
!#                             gfsc_buildJacobianGPBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of general 
!#         purpose limiter. For the first two routines, the velocity is assumed
!#         to be linear which simplifies the evaluation of the Jacobian matrix
!#         significantly. For the second two routines, the velocity can be arbitrary.
!#
!# 14.) gfsc_buildJacobianSymm = gfsc_buildJacobianSymmScalar /
!#                               gfsc_buildJacobianSymmBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of symmetric
!#         flux limiting for diffusion operators
!#
!# 15.) gfsc_buildFluxFCT = gfsc_buildFluxFCTScalar /
!#                          gfsc_buildFluxFCTBlock
!#     -> Assembles the raw antidiffusive flux for the FEM-FCT stabilisation
!#
!# The following auxiliary routines are available:
!#
!# 1.) gfsc_hasOrientation
!#     -> Checks if the stabilisation technique requires an oriented structure
!#
!# </purpose>
!##############################################################################

module groupfemscalar

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
  
  public :: gfsc_initStabilisation
  public :: gfsc_isMatrixCompatible
  public :: gfsc_isVectorCompatible
  public :: gfsc_buildConvectionOperator
  public :: gfsc_buildDiffusionOperator
  public :: gfsc_buildConvVectorFCT
  public :: gfsc_buildConvVectorTVD
  public :: gfsc_buildConvVectorGP
  public :: gfsc_buildConvVectorSymm
  public :: gfsc_buildConvectionJacobian
  public :: gfsc_buildJacobianFCT
  public :: gfsc_buildJacobianTVD
  public :: gfsc_buildJacobianGP
  public :: gfsc_buildJacobianSymm
  public :: gfsc_buildFluxFCT
  
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface gfsc_buildConvectionOperator
    module procedure gfsc_buildConvOperatorScalar
    module procedure gfsc_buildConvOperatorBlock
  end interface

  interface gfsc_buildConvVectorFCT
    module procedure gfsc_buildConvVecFCTScalar
    module procedure gfsc_buildConvVecFCTBlock
  end interface

  interface gfsc_buildConvVectorTVD
    module procedure gfsc_buildConvVecTVDScalar
    module procedure gfsc_buildConvVecTVDBlock
  end interface

  interface gfsc_buildConvVectorGP
    module procedure gfsc_buildConvVecGPScalar
    module procedure gfsc_buildConvVecGPBlock
  end interface

  interface gfsc_buildConvVectorSymm
    module procedure gfsc_buildConvVecSymmScalar
    module procedure gfsc_buildConvVecSymmBlock
  end interface

  interface gfsc_buildConvectionJacobian
    module procedure gfsc_buildConvJacobianScalar
    module procedure gfsc_buildConvJacobianBlock
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

  interface gfsc_buildFluxFCT
    module procedure gfsc_buildFluxFCTScalar
    module procedure gfsc_buildFluxFCTBlock
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_initStabilisation(rmatrixTemplate, rafcstab)

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
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: i

    
    ! Set atomic data
    rafcstab%NVARtransformed = rmatrixTemplate%NVAR
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = int(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ))
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)
      
      ! Handle for IverticesAtEdge
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      
    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_IMPLICIT,&
          AFCSTAB_FEMFCT_ITERATIVE)
      
      ! Handle for IverticesAtEdge
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
      
      ! Handle for DcoefficientsAtEdge
      Isize = (/3, rafcstab%NEDGE/)
      if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_DcoefficientsAtEdge)
      call storage_new('gfsc_initStabilisation', 'DcoefficientsAtEdge',&
          Isize, ST_DOUBLE, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .false., ST_DOUBLE)
      end do

      ! Associated vectors
      rafcstab%p_rvectorPp => rafcstab%RnodalVectors(1)
      rafcstab%p_rvectorPm => rafcstab%RnodalVectors(2)
      rafcstab%p_rvectorQp => rafcstab%RnodalVectors(3)
      rafcstab%p_rvectorQm => rafcstab%RnodalVectors(4)
      rafcstab%p_rvectorRp => rafcstab%RnodalVectors(5)
      rafcstab%p_rvectorRm => rafcstab%RnodalVectors(6)

      ! We need 1 nodal vector for the predictor
      allocate(rafcstab%RnodalBlockVectors(1))
      call lsysbl_createVectorBlock(&
          rafcstab%RnodalBlockVectors(1),&
          rafcstab%NEQ, 1, .false., ST_DOUBLE)
      
      ! Associate vector
      rafcstab%p_rvectorPredictor => rafcstab%RnodalBlockVectors(1)
      
      if (rafcstab%bprelimiting .or.&
          rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        ! We need 4 edgewise vectors for the correction factors and
        ! for the raw antidiffusive fluxes and the constraints
        allocate(rafcstab%RedgeVectors(4))
        rafcstab%p_rvectorPrelimit  => rafcstab%RedgeVectors(4)
        
        do i = 1, 4
          call lsyssc_createVector(rafcstab%RedgeVectors(i),&
              rafcstab%NEDGE, .false., ST_DOUBLE)
        end do
      else
        ! We need 3 edgewise vectors for the correction factors
        ! and for the raw antidiffusive fluxes
        allocate(rafcstab%RedgeVectors(3))
        do i = 1, 3
          call lsyssc_createVector(rafcstab%RedgeVectors(i),&
              rafcstab%NEDGE, .false., ST_DOUBLE)
        end do
      end if

      ! Associate vectors
      rafcstab%p_rvectorAlpha => rafcstab%RedgeVectors(1)
      rafcstab%p_rvectorFlux0 => rafcstab%RedgeVectors(2)
      rafcstab%p_rvectorFlux  => rafcstab%RedgeVectors(3)

      
    case (AFCSTAB_FEMTVD,&
          AFCSTAB_FEMGP)
      
      ! Handle for IverticesAtEdge
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
      
      ! Handle for DcoefficientsAtEdge
      Isize = (/3, rafcstab%NEDGE/)
      if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_DcoefficientsAtEdge)
      call storage_new('gfsc_initStabilisation', 'DcoefficientsAtEdge',&
          Isize, ST_DOUBLE, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .false., ST_DOUBLE)
      end do
      
      ! Associated vectors
      rafcstab%p_rvectorPp => rafcstab%RnodalVectors(1)
      rafcstab%p_rvectorPm => rafcstab%RnodalVectors(2)
      rafcstab%p_rvectorQp => rafcstab%RnodalVectors(3)
      rafcstab%p_rvectorQm => rafcstab%RnodalVectors(4)
      rafcstab%p_rvectorRp => rafcstab%RnodalVectors(5)
      rafcstab%p_rvectorRm => rafcstab%RnodalVectors(6)
      
      ! We need 3 edgewise vectors for the correction factors
      ! and for the raw antidiffusive fluxes
      allocate(rafcstab%RedgeVectors(3))
      do i = 1, 3
        call lsyssc_createVector(rafcstab%RedgeVectors(i),&
            rafcstab%NEDGE, .false., ST_DOUBLE)
      end do
      
      ! Associated vectors
      rafcstab%p_rvectorAlpha => rafcstab%RedgeVectors(1)
      rafcstab%p_rvectorFlux0 => rafcstab%RedgeVectors(2)
      rafcstab%p_rvectorFlux  => rafcstab%RedgeVectors(3)
      

    case (AFCSTAB_FEMFCT_LINEARISED)
      
      ! Handle for IverticesAtEdge
      Isize = (/4, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
      
      ! Handle for DcoefficientsAtEdge
      Isize = (/3, rafcstab%NEDGE/)
      if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_DcoefficientsAtEdge)
      call storage_new('gfsc_initStabilisation', 'DcoefficientsAtEdge',&
          Isize, ST_DOUBLE, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .false., ST_DOUBLE)
      end do

      ! Associated vectors
      rafcstab%p_rvectorPp => rafcstab%RnodalVectors(1)
      rafcstab%p_rvectorPm => rafcstab%RnodalVectors(2)
      rafcstab%p_rvectorQp => rafcstab%RnodalVectors(3)
      rafcstab%p_rvectorQm => rafcstab%RnodalVectors(4)
      rafcstab%p_rvectorRp => rafcstab%RnodalVectors(5)
      rafcstab%p_rvectorRm => rafcstab%RnodalVectors(6)
      
      ! We need 2-3 edgewise vectors at most
      if (rafcstab%bprelimiting) then
        allocate(rafcstab%RedgeVectors(3))
      else
        allocate(rafcstab%RedgeVectors(2))
      end if

      ! We need 2 edgewise vectors for the edgewise correction
      ! factors and the raw antidiffusive fluxes
      do i = 1, 2
        call lsyssc_createVector(rafcstab%RedgeVectors(i),&
            rafcstab%NEDGE, .false., ST_DOUBLE)
      end do

      ! Associated vectors
      rafcstab%p_rvectorAlpha => rafcstab%RedgeVectors(1)
      rafcstab%p_rvectorFlux => rafcstab%RedgeVectors(2)

      ! If the raw antidiffusive fluxes should be prelimited
      ! then 1 additional edgewie vector is required
      if (rafcstab%bprelimiting) then
        call lsyssc_createVector(rafcstab%RedgeVectors(3),&
            rafcstab%NEDGE, .false., ST_DOUBLE)
        rafcstab%p_rvectorPrelimit  => rafcstab%RedgeVectors(3)
      end if

      ! We need 1 nodal block-vector for the 
      ! approximation to the time derivative
      allocate(rafcstab%RnodalBlockVectors(1))
      call lsysbl_createVectorBlock(&
          rafcstab%RnodalBlockVectors(1),&
          rafcstab%NEQ, 1, .false., ST_DOUBLE)

      ! Associate vector
      rafcstab%p_rvectorPredictor  => rafcstab%RnodalBlockVectors(1)


    case (AFCSTAB_SYMMETRIC)

      ! Handle for IverticesAtEdge
      Isize = (/2, rafcstab%NEDGE/)
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      call storage_new('gfsc_initStabilisation', 'IverticesAtEdge',&
          Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)

      ! Handle for DcoefficientsAtEdge
      Isize = (/2, rafcstab%NEDGE/)
      if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_DcoefficientsAtEdge)
      call storage_new('gfsc_initStabilisation', 'DcoefficientsAtEdge',&
          Isize, ST_DOUBLE, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)

      ! We need 6 nodal vectors for P`s, Q`s and R`s
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
            rafcstab%NEQ, .false., ST_DOUBLE)
      end do
      
      ! Associated vectors
      rafcstab%p_rvectorPp => rafcstab%RnodalVectors(1)
      rafcstab%p_rvectorPm => rafcstab%RnodalVectors(2)
      rafcstab%p_rvectorQp => rafcstab%RnodalVectors(3)
      rafcstab%p_rvectorQm => rafcstab%RnodalVectors(4)
      rafcstab%p_rvectorRp => rafcstab%RnodalVectors(5)
      rafcstab%p_rvectorRm => rafcstab%RnodalVectors(6)

      ! We need 1 edgewise vector for the fluxes
      allocate(rafcstab%RedgeVectors(1))
      call lsyssc_createVector(rafcstab%RedgeVectors(1),&
          rafcstab%NEDGE, .false., ST_DOUBLE)

      ! Associate vector
      rafcstab%p_rvectorFlux  => rafcstab%RedgeVectors(1)

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_initStabilisation')
      call sys_halt()
    end select

    ! Set specifier
    rafcstab%iSpec = AFCSTAB_INITIALISED

  end subroutine gfsc_initStabilisation

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_isMatrixCompatible(rafcstab, rmatrix, bcompatible)

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
        rafcstab%NEDGE .ne. int(0.5*(rmatrix%NA-rmatrix%NEQ),I32)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_isMatrixCompatible')
        call sys_halt()
      end if
    end if
  end subroutine gfsc_isMatrixCompatible

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_isVectorCompatible(rafcstab, rvector, bcompatible)

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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_isVectorCompatible')
        call sys_halt()
      end if
    end if
  end subroutine gfsc_isVectorCompatible

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvOperatorBlock(RcoeffMatrices, rafcstab, rx,&
      fcb_calcMatrix, bbuildStabilisation, bclear, rconvMatrix,&
      bisConservative)
    
!<description>
    ! This subroutine assembles the discrete transport operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !   $$ \nabla\cdot({\bf v}u) $$
    !
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! The array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! The solution vector
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorBlock), intent(in) :: rx
    
    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for (non-)conservative matrix assembly
    ! TRUE  : assemble conservative convection operator (default)
    ! FALSE : assemble non-conservative convection operator
    logical, intent(in), optional :: bisConservative

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The transport operator
    type(t_matrixScalar), intent(inout) :: rconvMatrix
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvOperatorScalar(RcoeffMatrices, rafcstab,&
          rx%RvectorBlock(1), fcb_calcMatrix, bbuildStabilisation, bclear,&
          rconvMatrix, bisConservative)

    end if

  end subroutine gfsc_buildConvOperatorBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvOperatorScalar(RcoeffMatrices, rafcstab,&
      rx, fcb_calcMatrix, bbuildStabilisation, bclear, rconvMatrix,&
      bisConservative)

!<description>
    ! This subroutine assembles the discrete transport operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !     $$ \nabla\cdot({\bf v}u) $$
    !
    ! This routine can be used to apply the following discretisations:
    !
    ! (1) the standard Galerkin finite element method
    !     which will be referred to as high-order approximation
    !
    ! (2) the discrete upwinding operator which results from the conservative
    !     elimination of negative off-diagonal entries from the Galerkin
    !     operator. This technique is for instance described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    !     conservation laws, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 155-206.
    !
    ! (3) In addition to (2), discrete upwinding is performed and auxiliary
    !     data required for flux limiting is generated in the optional 
    !     stabilisation structure rafcstab. For symmetric flux limiters,
    !     the artificil diffusio coefficient d_ij and the entries of the
    !     low-order operator l_ij=k_ij+d_ij are stored as well as the node
    !     numbers i and j and the matrix positions ij/ji of the edge (i,j).
    !     For upwind-biased flux limiters the same data are stored but the
    !     following orientation convention is applied. Node i is located 
    !     upwind and corresponds to the column number whose entry has been
    !     eliminated.
!</description>

!<input>
    ! The array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! The solution vector
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorScalar), intent(in) :: rx

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for (non-)conservative matrix assembly
    ! TRUE  : assemble conservative convection operator (default)
    ! FALSE : assemble non-conservative convection operator
    logical, intent(in), optional :: bisConservative

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The transport operator
    type(t_matrixScalar), intent(inout) :: rconvMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_CoeffX,p_CoeffY,p_CoeffZ,p_ConvOp,p_Dx
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: ndim
    logical :: bconservative


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    ! Let us check if the edge-based data structure has been generated
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(RcoeffMatrices(1), rafcstab)
    end if

    ! Check if conservative of non-conservative convection operator is required
    bconservative = .true.
    if (present(bisConservative)) bconservative = bisConservative
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rconvMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsyssc_getbase_double(rconvMatrix, p_ConvOp)
    call lsyssc_getbase_double(rx, p_Dx)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_CoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      call sys_halt()
    end select
    

    ! What kind of matrix are we?
    select case(rconvMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 7 and 9
      !-------------------------------------------------------------------------
      
      ! Set diagonal pointer
      if (rconvMatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then
        call lsyssc_getbase_Kld(RcoeffMatrices(1), p_Kdiagonal)
      else
        call lsyssc_getbase_Kdiagonal(RcoeffMatrices(1), p_Kdiagonal)
      end if

      ! Do we have to perform stabilisation?
      if (bbuildStabilisation) then
        
        ! Set additional pointers
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
        
        ! Do we need edge orientation?
        if (gfsc_hasOrientation(rafcstab)) then
          
          ! Adopt orientation convention IJ, such that L_ij < L_ji
          ! and generate edge structure for the flux limiter
          
          if (bconservative) then
            
            ! Conservative formulation of convection operator
            
            select case(ndim)
            case (NDIM1D)
              call doOperatorOAFCMat79Cons1D(p_Kdiagonal,&
                  rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_Dx, p_ConvOp,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doOperatorOAFCMat79Cons2D(p_Kdiagonal,&
                  rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_Dx, p_ConvOp,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doOperatorOAFCMat79Cons3D(p_Kdiagonal,&
                  rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_ConvOp,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else

            ! Non-conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doOperatorOAFCMat79Nonc1D(p_Kdiagonal,&
                  rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_Dx, p_ConvOp,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doOperatorOAFCMat79Nonc2D(p_Kdiagonal,&
                  rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_Dx, p_ConvOp,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doOperatorOAFCMat79Nonc3D(p_Kdiagonal,&
                  rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_ConvOp,&
                  p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION)

        else   ! bhasOrientation == no

          ! Adopt no orientation convention and generate edge structure
          
          if (bconservative) then

            ! Conservative formulation of convection operator
            
            select case(ndim)
            case (NDIM1D)
              call doOperatorAFCMat79Cons1D(p_Kdiagonal,&
                  p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_Dx, p_ConvOp, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doOperatorAFCMat79Cons2D(p_Kdiagonal,&
                  p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_Dx, p_ConvOp, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doOperatorAFCMat79Cons3D(p_Kdiagonal,&
                  p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_ConvOp, p_DcoefficientsAtEdge)
            end select

          else

            ! Non-conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doOperatorAFCMat79Nonc1D(p_Kdiagonal,&
                  p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_Dx, p_ConvOp, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doOperatorAFCMat79Nonc2D(p_Kdiagonal,&
                  p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_Dx, p_ConvOp, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doOperatorAFCMat79Nonc3D(p_Kdiagonal,&
                  p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                  p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_ConvOp, p_DcoefficientsAtEdge)
            end select
            
          end if
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)
          rafcstab%iSpec = iand(rafcstab%iSpec, not(AFCSTAB_HAS_EDGEORIENTATION))

        end if   ! bhasOrientation
        
      else   ! bbuildStabilisation == no
        
        ! Apply standard discretisation without stabilisation
        
        if (bconservative) then
          
          ! Conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doOperatorMat79Cons1D(p_Kdiagonal,&
                p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                p_CoeffX, p_Dx, p_ConvOp)
          case (NDIM2D)
            call doOperatorMat79Cons2D(p_Kdiagonal,&
                p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                p_CoeffX, p_CoeffY, p_Dx, p_ConvOp)
          case (NDIM3D)
            call doOperatorMat79Cons3D(p_Kdiagonal,&
                p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_ConvOp)
          end select
          
        else
          
          ! Non-conservative formulation of convection operator
          
          select case(ndim)
          case (NDIM1D)
            call doOperatorMat79Nonc1D(p_Kdiagonal,&
                p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                p_CoeffX, p_Dx, p_ConvOp)
          case (NDIM2D)
            call doOperatorMat79Nonc2D(p_Kdiagonal,&
                p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                p_CoeffX, p_CoeffY, p_Dx, p_ConvOp)
          case (NDIM3D)
            call doOperatorMat79Nonc3D(p_Kdiagonal,&
                p_IverticesAtEdge, rafcstab%NEDGE, rconvMatrix%NEQ,&
                p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_ConvOp)
          end select
          
        end if
        
        ! Set state of stabilisation
        rafcstab%iSpec = iand(rafcstab%iSpec, not(AFCSTAB_HAS_EDGEVALUES))
        rafcstab%iSpec = iand(rafcstab%iSpec, not(AFCSTAB_HAS_EDGEORIENTATION))

      end if
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routine follow
        
    !**************************************************************
    ! Assemble convection operator operator L in 1D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation
    
    subroutine doOperatorMat79Cons1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, Dx, L)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,d_ij
      integer :: iedge,ii,jj,ij,ji,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ii)

        ! Compute coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion (if any)
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
      end do

    end subroutine doOperatorMat79Cons1D

    !**************************************************************
    ! Assemble convection operator operator L in 1D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation
    
    subroutine doOperatorMat79Nonc1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, Dx, L)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,d_ij
      integer :: iedge,ii,jj,ij,ji,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion (if any)
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
      end do

    end subroutine doOperatorMat79Nonc1D
    
    !**************************************************************
    ! Assemble convection operator operator L in 2D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation
    
    subroutine doOperatorMat79Cons2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, Dx, L)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,d_ij
      integer :: iedge,ii,jj,ij,ji,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ii); C_ij(2) = CoeffY(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion (if any)
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
      end do

    end subroutine doOperatorMat79Cons2D

    !**************************************************************
    ! Assemble convection operator operator L in 2D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation
    
    subroutine doOperatorMat79Nonc2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, Dx, L)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,d_ij
      integer :: iedge,ii,jj,ij,ji,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion (if any)
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
      end do

    end subroutine doOperatorMat79Nonc2D

    
    !**************************************************************
    ! Assemble convection operator operator L in 3D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation
    
    subroutine doOperatorMat79Cons3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, CoeffZ, Dx, L)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP):: k_ij,k_ji,d_ij
      integer :: iedge,ii,jj,ij,ji,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ii); C_ij(2) = CoeffY(ii); C_ij(3) = CoeffZ(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion (if any)
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
      end do

    end subroutine doOperatorMat79Cons3D

    !**************************************************************
    ! Assemble convection operator operator L in 3D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation
    
    subroutine doOperatorMat79Nonc3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, CoeffZ, Dx, L)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP):: k_ij,k_ji,d_ij
      integer :: iedge,ii,jj,ij,ji,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

        ! Apply artificial diffusion (if any)
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
      end do

    end subroutine doOperatorMat79Nonc3D
    
    !**************************************************************
    ! Assemble convection operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation
    
    subroutine doOperatorAFCMat79Cons1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, Dx, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficient
        C_ij(1) = CoeffX(ii)
        
        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
      end do
      
    end subroutine doOperatorAFCMat79Cons1D

    !**************************************************************
    ! Assemble convection operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation
    
    subroutine doOperatorAFCMat79Nonc1D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, Dx, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
      end do
      
    end subroutine doOperatorAFCMat79Nonc1D

    !**************************************************************
    ! Assemble convection operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation

    subroutine doOperatorAFCMat79Cons2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, Dx, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ii); C_ij(2) = CoeffY(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
      end do
      
    end subroutine doOperatorAFCMat79Cons2D

    !**************************************************************
    ! Assemble convection operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation

    subroutine doOperatorAFCMat79Nonc2D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, Dx, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ
      
      real(DP), dimension(:), intent(inout) :: L

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
      end do
      
    end subroutine doOperatorAFCMat79Nonc2D
    
    !**************************************************************
    ! Assemble convection operator L and AFC data w/o edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation

    subroutine doOperatorAFCMat79Cons3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, CoeffZ, Dx, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      

      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ii); C_ij(2) = CoeffY(ii); C_ij(3) = CoeffZ(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
      end do
      
    end subroutine doOperatorAFCMat79Cons3D

    !**************************************************************
    ! Assemble convection operator L and AFC data w/o edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation

    subroutine doOperatorAFCMat79Nonc3D(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, CoeffX, CoeffY, CoeffZ, Dx, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
            
      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
       
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
 
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
      end do
      
    end subroutine doOperatorAFCMat79Nonc3D
    
    !**************************************************************
    ! Assemble convection operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation

    subroutine doOperatorOAFCMat79Cons1D(Kdiagonal, NEDGE, NEQ,&
        CoeffX, Dx, L, IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = CoeffX(ii)
        
        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
        
        ! AFC with edge orientation
        if (k_ij < k_ji) then
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
        else
          IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
        end if
      end do

    end subroutine doOperatorOAFCMat79Cons1D

    !**************************************************************
    ! Assemble convection operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation

    subroutine doOperatorOAFCMat79Nonc1D(Kdiagonal, NEDGE, NEQ,&
        CoeffX, Dx, L, IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
        
        ! AFC with edge orientation
        if (k_ij < k_ji) then
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
        else
          IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
        end if
      end do

    end subroutine doOperatorOAFCMat79Nonc1D
    
    !**************************************************************
    ! Assemble convection operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation

    subroutine doOperatorOAFCMat79Cons2D(Kdiagonal, NEDGE, NEQ,&
        CoeffX, CoeffY, Dx, L, IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j

      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ij(1) = CoeffX(ii); C_ij(2) = CoeffY(ii)
        
        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ji, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
        
        ! AFC with edge orientation
        if (k_ij < k_ji) then
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
        else
          IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
        end if
      end do

    end subroutine doOperatorOAFCMat79Cons2D

    !**************************************************************
    ! Assemble convection operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation

    subroutine doOperatorOAFCMat79Nonc2D(Kdiagonal, NEDGE, NEQ,&
        CoeffX, CoeffY, Dx, L, IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
        
        ! AFC with edge orientation
        if (k_ij < k_ji) then
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
        else
          IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
        end if
      end do

    end subroutine doOperatorOAFCMat79Nonc2D
    
    !**************************************************************
    ! Assemble convection operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7 and 9
    ! Conservative formulation

    subroutine doOperatorOAFCMat79Cons3D(Kdiagonal, NEDGE, NEQ,&
        CoeffX, CoeffY, CoeffZ, Dx, L, IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE, NEQ

      real(DP), dimension(:), intent(inout) :: L
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Loop over all rows
      !$omp parallel do private(ii,C_ij,k_ij,k_ji,d_ij)
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ii); C_ij(2) = CoeffY(ii); C_ij(3) = CoeffZ(ii)
        
        ! Compute convection coefficients for diagonal
        call fcb_calcMatrix(Dx(i), Dx(i),&
            C_ij, C_ij, i, i, k_ij, k_ji, d_ij)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ij
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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - d_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - d_ij
          
        ! AFC with edge orientation
        if (k_ij < k_ji) then
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
        else
          IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
        end if
      end do

    end subroutine doOperatorOAFCMat79Cons3D

    !**************************************************************
    ! Assemble convection operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7 and 9
    ! Non-conservative formulation

    subroutine doOperatorOAFCMat79Nonc3D(Kdiagonal, NEDGE, NEQ,&
        CoeffX, CoeffY, CoeffZ, Dx, L, IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE, NEQ

      real(DP), dimension(:), intent(inout) :: L
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      

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
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute convection coefficients
        call fcb_calcMatrix(Dx(i), Dx(j),&
            C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
        
        ! Apply artificial diffusion
        k_ij = k_ij + d_ij
        k_ji = k_ji + d_ij
        
        ! Assemble the global operator
        L(ii) = L(ii) - k_ij
        L(ij) = L(ij) + k_ij 
        L(ji) = L(ji) + k_ji
        L(jj) = L(jj) - k_ji
        
        ! AFC with edge orientation
        if (k_ij < k_ji) then
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
        else
          IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
        end if
      end do
      
    end subroutine doOperatorOAFCMat79Nonc3D
    
  end subroutine gfsc_buildConvOperatorScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildDiffusionOperator(rcoeffMatrix, rafcstab,&
      bbuildStabilisation, bclear, rdiffMatrix)

!<description>
    ! This subroutine assembles the diffusive part of the discrete
    ! transport operator which results from the discretisation of the
    ! scalar convection-diffusion-reaction equation.  The matrix
    ! rmatrixS holds the unmodified diffusion operator which is
    ! possible modified by this routine. If the optional argument
    ! rmatrixDest is given, then the content of rmatrixS is first
    ! copied/added to the destination matrix prior to performing
    ! further modifications.  If the parameter bclear is TRUE the
    ! destination matrix is cleared, that is, the content of rmatrixS
    ! is copied. Otherwise, the its content is combined linearly with
    ! that of the destination matrix.  If the parameter bStabilse is
    ! TRUE, then symmetric stabilisation is applied so that the
    ! resulting diffusion operator is guaranted to ensure the discrete
    ! maximum principle.  If the optional parameter rafcstab is given,
    ! then the auxiliary data structures for of the stabilisation
    ! structure are updated.
!</description>

!<input>
    ! (anisotropic) diffusion operator
    type(t_matrixScalar), intent(in) :: rcoeffMatrix
    
    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab

    ! diffusion operator
    type(t_matrixScalar), intent(inout) :: rdiffMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_S,p_DiffOp
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffusionOperator')
      call sys_halt()
    end if
    
    ! Check if stabilisation provides edge-based structure
    ! Let us check if the edge-based data structure has been generated
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(rcoeffMatrix, rafcstab)
    end if
    
    ! Should matrix be cleared?
    if (bclear) call lsyssc_clearMatrix(rdiffMatrix)

    ! Set pointers
     call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call lsyssc_getbase_double(rcoeffMatrix, p_S)
    call lsyssc_getbase_double(rdiffMatrix, p_DiffOp)      

    
    ! What kind of matrix are we?
    select case(rdiffMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 7 and 9
      !-------------------------------------------------------------------------
      
      ! Set diagonal pointer
      if (rdiffMatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then
        call lsyssc_getbase_Kld(rcoeffMatrix, p_Kdiagonal)
      else
        call lsyssc_getbase_Kdiagonal(rcoeffMatrix, p_Kdiagonal)
      end if
      
      ! Do we have to perform stabilisation?
      if (bbuildStabilisation) then
        
        ! Set additional pointers
        call afcstab_getbase_DcoeffsAtEdge(rafcstab,&
            p_DcoefficientsAtEdge)
        
        call doOperatorAFCMat79(p_Kdiagonal, p_IverticesAtEdge,&
            rafcstab%NEDGE, rdiffMatrix%NEQ, p_S, p_DiffOp,&
            p_DcoefficientsAtEdge)
        
        ! Set state of stabilisation
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)

      else
        
        call doOperatorMat79(p_Kdiagonal, p_IverticesAtEdge,&
            rafcstab%NEDGE, rdiffMatrix%NEQ, p_S, p_DiffOp)

        ! Set state of stabilisation
        rafcstab%iSpec = iand(rafcstab%iSpec, not(AFCSTAB_HAS_EDGEVALUES))
        
      end if
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffusionOperator')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble low-order diffusion operator S.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorMat79(Kdiagonal, IverticesAtEdge, NEDGE, NEQ, S, L)

      real(DP), dimension(:), intent(in) :: S
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L

      ! local variables
      real(DP) :: d_ij
      integer :: iedge,ii,ij,ji,jj,i,j
      

      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)
        
        ! Artificial diffusion coefficient
        d_ij = max(0.0_DP, -S(ij)) 
        
        ! Assemble the global operator
        L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
        L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij
      end do
    end subroutine doOperatorMat79

    !**************************************************************
    ! Assemble low-order diffusion operator S and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorAFCMat79(Kdiagonal, IverticesAtEdge,&
        NEDGE, NEQ, S, L, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(in) :: S
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE,NEQ

      real(DP), dimension(:), intent(inout) :: L
      
      real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
      
      ! local variables
      real(DP) :: d_ij,s_ij
      integer :: ii,ij,ji,jj,iedge,i,j
      

      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)
        ii = Kdiagonal(i)
        jj = Kdiagonal(j)
        
        ! Artificial diffusion coefficient
        d_ij = max(0.0_DP, -S(ij))
        s_ij = max(0.0_DP,  S(ij))
        
        ! Assemble the global operator
        L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
        L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij
        
        ! AFC w/o edge orientation
        DcoefficientsAtEdge(:,iedge) = (/d_ij, s_ij/)
      end do

    end subroutine doOperatorAFCMat79

  end subroutine gfsc_buildDiffusionOperator
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecFCTBlock(rlumpedMassMatrix,&
      rafcstab, rx, dscale, bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-FCT type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
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
      
      call gfsc_buildConvVecFCTScalar(rlumpedMassMatrix,&
          rafcstab, rx%RvectorBlock(1), dscale, bclear,&
          ioperationSpec, ry%RvectorBlock(1))
      
    end if
    
  end subroutine gfsc_buildConvVecFCTBlock

  ! *****************************************************************************
  
!<subroutine>
  
  subroutine gfsc_buildConvVecFCTScalar(rlumpedMassMatrix,&
      rafcstab, rx, dscale, bclear, ioperationSpec, ry)

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
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Check if stabilisation is prepeared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
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
    ! 4) Apply the limited antiddifusive fluxes to the convective term
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
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
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
        call lsyssc_getbase_double(rafcstab%p_rvectorPrelimit, p_Dflux)
      else
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      end if

      ! Compute sums of antidiffusive increments
      if (rafcstab%bprelimiting) then
        call lsyssc_getbase_double(rafcstab%p_rvectorPrelimit, p_Dflux0)
        call doPreADIncrements(p_IverticesAtEdge,&
            rafcstab%NEDGE, p_Dflux, p_Dflux0, p_Dalpha, p_Dpp, p_Dpm)
      else
        call doADIncrements(p_IverticesAtEdge,&
            rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if
      
      ! Set pointers
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
      call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute bounds
      call doBounds(p_IverticesAtEdge,&
          rafcstab%NEDGE, p_Dx, p_Dqp, p_Dqm)

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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
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
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        call doLimitNodal(rafcstab%NEQ, dscale,&
            p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        call doLimitNodalConstrained(rafcstab%NEQ, dscale,&
            p_ML, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
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
        call output_line('Stabilisation does not provide nodal correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
      call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_FEMFCT_IMPLICIT) then
        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorPrelimit, p_Dflux0)
        call doLimitEdgewiseConstrained(p_IverticesAtEdge,&
            rafcstab%NEDGE, p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
      else
        call doLimitEdgewise(p_IverticesAtEdge,&
            rafcstab%NEDGE, p_Dflux, p_Drp, p_Drm, p_Dalpha)
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
        call output_line('Stabilisation does not provide edgewise correction factors',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Set pointers
      call lsyssc_getbase_double(ry, p_Dy)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

      ! Clear convective vector?
      if (bclear) call lsyssc_clearVector(ry)

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then
        call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
        call doCorrectScaleByMass(p_IverticesAtEdge,&
            rafcstab%NEDGE, dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
      else
        call doCorrect(p_IverticesAtEdge,&
            rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Dy)
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble sums of antidiffusive increments
    ! for the given antidiffusive fluxes
    
    subroutine doADIncrements(IverticesAtEdge,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)
      
      real(DP), dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(:), intent(out) :: Dpp,Dpm
      
      ! local variables
      real(DP) :: f_ij
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
        f_ij = Dalpha(iedge)*Dflux(iedge)
        
        ! Compute the sums of antidiffusive increments
        Dpp(i) = Dpp(i)+max(0.0_DP, f_ij)
        Dpp(j) = Dpp(j)+max(0.0_DP,-f_ij)
        Dpm(i) = Dpm(i)+min(0.0_DP, f_ij)
        Dpm(j) = Dpm(j)+min(0.0_DP,-f_ij)
      end do
    end subroutine doADIncrements

    !**************************************************************
    ! Assemble sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and with prelimiting

    subroutine doPreADIncrements(IverticesAtEdge,&
        NEDGE, Dflux, Dflux0, Dalpha, Dpp, Dpm)
      
      real(DP), dimension(:), intent(in) :: Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Dalpha
      real(DP), dimension(:), intent(out) :: Dpp,Dpm
      
      ! local variables
      real(DP) :: f_ij,alpha_ij
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
        f_ij = Dalpha(iedge)*Dflux(iedge)
        
        ! MinMod prelimiting
        alpha_ij = mprim_minmod3(f_ij, Dflux0(iedge), f_ij)

        ! Synchronisation of correction factors
        Dalpha(iedge) = Dalpha(iedge) * alpha_ij

        ! Update the raw antidiffusive Dflux
        F_ij = alpha_ij * F_ij

        ! Compute the sums of antidiffusive increments
        Dpp(i) = Dpp(i)+max(0.0_DP, f_ij)
        Dpp(j) = Dpp(j)+max(0.0_DP,-f_ij)
        Dpm(i) = Dpm(i)+min(0.0_DP, f_ij)
        Dpm(j) = Dpm(j)+min(0.0_DP,-f_ij)
      end do
    end subroutine doPreADIncrements

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    
    subroutine doBounds(IverticesAtEdge, NEDGE, Dx, Dqp, Dqm)
      
      real(DP), dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(out) :: Dqp,Dqm
      
      ! local variables
      real(DP) :: Dx_ij
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
        Dx_ij = Dx(j)-Dx(i)
        
        ! Compute the distance to a local extremum
        ! of the predicted solution
        Dqp(i) = max(Dqp(i), Dx_ij)
        Dqp(j) = max(Dqp(j),-Dx_ij)
        Dqm(i) = min(Dqm(i), Dx_ij)
        Dqm(j) = min(Dqm(j),-Dx_ij)
      end do
    end subroutine doBounds

    !**************************************************************
    ! Compute nodal correction factors without constraints
    
    subroutine doLimitNodal(NEQ, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)
      
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ
      
      real(DP), dimension(:), intent(inout) :: Drp,Drm
      
      ! local variables
      integer :: ieq

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drp(ieq) = ML(ieq)*Dqp(ieq)/(dscale*Dpp(ieq)+SYS_EPSREAL)
      end do
      !$omp end parallel do

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drm(ieq) = ML(ieq)*Dqm(ieq)/(dscale*Dpm(ieq)-SYS_EPSREAL)
      end do
      !$omp end parallel do
    end subroutine doLimitNodal

    !**************************************************************
    ! Compute nodal correction factors with constraints
    
    subroutine doLimitNodalConstrained(NEQ, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)
      
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ
      
      real(DP), dimension(:), intent(inout) :: Drp,Drm
      
      ! local variables
      integer :: ieq

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drp(ieq) = min(1.0_DP, ML(ieq)*Dqp(ieq)/(dscale*Dpp(ieq)+SYS_EPSREAL))
      end do
      !$omp end parallel do

      ! Loop over all vertices
      !$omp parallel do
      do ieq = 1, NEQ
        Drm(ieq) = min(1.0_DP, ML(ieq)*Dqm(ieq)/(dscale*Dpm(ieq)-SYS_EPSREAL))
      end do
      !$omp end parallel do
    end subroutine doLimitNodalConstrained    

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    
    subroutine doLimitEdgewise(IverticesAtEdge,&
        NEDGE, Dflux, Drp, Drm, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dalpha
      
      ! local variables
      real(DP) :: f_ij,r_ij
      integer :: iedge,i,j
      
      ! Loop over all edges
      !$omp parallel do private(i,j,f_ij,r_ij)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)
 
        ! Compute nodal correction factors
        if (f_ij .ge. 0.0_DP) then
          r_ij = min(Drp(i),Drm(j))
        else
          r_ij = min(Drp(j),Drm(i))
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do
    end subroutine doLimitEdgewise

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes
    
    subroutine doLimitEdgewiseConstrained(IverticesAtEdge,&
        NEDGE, Dflux1, Dflux2, Drp, Drm, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP) :: f1_ij,f2_ij,r_ij
      integer :: iedge,i,j
      
      ! Loop over all edges
      !$omp parallel do private(i,j,f_ij,r_ij)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
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
    end subroutine doLimitEdgewiseConstrained

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    
    subroutine doCorrect(IverticesAtEdge,&
        NEDGE, dscale, Dalpha, Dflux, Dy)
      
      real(DP), dimension(:), intent(in) :: Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dy
      
      ! local variables
      real(DP) :: f_ij
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Correct antidiffusive flux
        f_ij = dscale * Dalpha(iedge) * Dflux(iedge)
        
        ! Apply limited antidiffusive fluxes
        Dy(i) = Dy(i) + f_ij
        Dy(j) = Dy(j) - f_ij
      end do
    end subroutine doCorrect

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix
    
    subroutine doCorrectScaleByMass(IverticesAtEdge,&
        NEDGE, dscale, ML, Dalpha, Dflux, Dy)
      
      real(DP), dimension(:), intent(in) :: ML,Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dy
      
      ! local variables
      real(DP) :: f_ij
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Correct antidiffusive flux
        f_ij = dscale * Dalpha(iedge) * Dflux(iedge)
        
        ! Apply limited antidiffusive fluxes
        Dy(i) = Dy(i) + f_ij/ML(i)
        Dy(j) = Dy(j) - f_ij/ML(j)
      end do
    end subroutine doCorrectScaleByMass

  end subroutine gfsc_buildConvVecFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecTVDBlock(rx, tstep, ry, rafcstab)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-TVD type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! time step size
    real(DP), intent(in) :: tstep
!</input>

!<inputoutput>
    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry

    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks   .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDBlock')
      call sys_halt()

    else

      call gfsc_buildConvVecTVDScalar(rx%RvectorBlock(1),&
          tstep, ry%RvectorBlock(1), rafcstab)
      
    end if

  end subroutine gfsc_buildConvVecTVDBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecTVDScalar(rx, tstep, ry, rafcstab)

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

    ! time step size
    real(DP), intent(in) :: tstep
!</input>

!<inputoutput>
    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry

    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dx,p_Dy,p_Dflux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)  .eq.0) .or. &
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION).eq.0) .or. &
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)     .eq.0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDScalar')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
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
    call doLimit_TVD(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
        p_Dx, tstep, rafcstab%NEDGE, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
        p_Drp, p_Drm, p_Dflux, p_Dy)

    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)
    
  contains

    ! Here, the working routine follows
    
    !**************************************************************
    ! The FEM-TVD limiting procedure
    
    subroutine doLimit_TVD(IverticesAtEdge, DcoefficientsAtEdge,&
        Dx, tstep, NEDGE, Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dy)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: tstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm,Dflux,Dy

      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff
      integer :: iedge,ij,i,j
      
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(Dpp)
      call lalg_clearVectorDble(Dpm)
      call lalg_clearVectorDble(Dqp)
      call lalg_clearVectorDble(Dqm)

      ! Assemble P`s and Q`s
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        l_ij = DcoefficientsAtEdge(2,iedge)
        l_ji = DcoefficientsAtEdge(3,iedge)
        
        ! Determine solution difference
        diff = tstep*(Dx(i)-Dx(j))
        
        ! Prelimit the antidiffusive flux F`_IJ=MIN(-P_IJ,L_JI)(DX_I-DX_J)
        f_ij = min(d_ij,l_ji)*diff; Dflux(iedge) = f_ij
        
        ! Assemble P`s accordingly
        Dpp(i) = Dpp(i)+max(0.0_DP, f_ij)
        Dpm(i) = Dpm(i)+min(0.0_DP, f_ij)
        
        ! Assemble Q`s
        Dqp(i) = Dqp(i)+max(0.0_DP,-f_ij)
        Dqp(j) = Dqp(j)+max(0.0_DP, f_ij)
        Dqm(i) = Dqm(i)+min(0.0_DP,-f_ij)
        Dqm(j) = Dqm(j)+min(0.0_DP, f_ij)
      end do
      
      ! Apply the nodal limiter
      Drp = afcstab_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)

      ! Apply limiter
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed raw antidiffusive flux
        f_ij = Dflux(iedge)
        
        ! Apply correction factor and store limite flux
        f_ij = merge(Drp(i), Drm(i), f_ij > 0)*f_ij
        
        ! Update the vector
        Dy(i) = Dy(i)+f_ij
        Dy(j) = Dy(j)-f_ij
      end do
    end subroutine doLimit_TVD
    
  end subroutine gfsc_buildConvVecTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecGPBlock(rconsistentMassMatrix, rx, rx0,&
      theta, tstep, ry, rafcstab)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-GP type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep
!</input>

!<inputoutput>
    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry

    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
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

      call gfsc_buildConvVecGPScalar(rconsistentMassMatrix,&
          rx%RvectorBlock(1), rx0%RvectorBlock(1),&
          theta, tstep, ry%RvectorBlock(1), rafcstab)
      
    end if
  end subroutine gfsc_buildConvVecGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecGPScalar(rconsistentMassMatrix, rx, rx0,&
      theta, tstep, ry, rafcstab)

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
    type(t_matrixScalar), intent(in) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! initial solution vector
    type(t_vectorScalar), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep
!</input>

!<inputoutput>
    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry

    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_MC,p_Dx,p_Dx0,p_Dy,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or. &
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or. &
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecGPScalar')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(rx0, p_Dx0)
    call lsyssc_getbase_double(ry, p_Dy)

    ! Perform flux limiting by the general purpose limiter
    call doLimit_GP(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
        p_MC, p_Dx, p_Dx0, theta, tstep, rafcstab%NEDGE,&
        p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm,&
        p_Dflux, p_Dflux0, p_Dy)
    
    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)
    
  contains

    ! Here, the working routine follows
    
    !**************************************************************
    ! The FEM-GP limiting procedure
    
    subroutine doLimit_GP(IverticesAtEdge, DcoefficientsAtEdge,&
        MC, Dx, Dx0, theta, tstep, NEDGE, Dpp, Dpm, Dqp, Dqm,&
        Drp, Drm, Dflux, Dflux0, Dy)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0
      real(DP), intent(in) :: theta,tstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(:), intent(inout) :: Dflux,Dflux0,Dy
      
      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,m_ij,p_ij,pf_ij,df_ij,q_ij,q_ji
      real(DP) :: diff,diff0,diff1
      integer :: iedge,ij,i,j
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(Dpp)
      call lalg_clearVectorDble(Dpm)
      call lalg_clearVectorDble(Dqp)
      call lalg_clearVectorDble(Dqm)
      
      ! Assemble P`s and Q`s
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        l_ij = DcoefficientsAtEdge(2,iedge)
        l_ji = DcoefficientsAtEdge(3,iedge)
        m_ij = MC(ij)
        
        ! Compute: diff1 = dt*theta*(Dx_i-Dx_j) + dt*(1-theta)*(Dx0_i-Dx0_j)
        diff1 = Dx(i)-Dx(j); diff0 = Dx0(i)-Dx0(j)
        diff  = tstep*(theta*diff1+(1.0_DP-theta)*diff0)
        
        ! Compute antidiffusive flux f_ij=min(0,p_ij)*(Dx_j-Dx_i)
        if (abs(diff) < SYS_EPSREAL) then
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
        
        q_ij = m_ij/tstep+l_ij
        q_ji = m_ij/tstep+l_ji

        ! Assemble Q`s
        Dqp(i) = Dqp(i)+q_ij*max(0.0_DP,-diff)
        Dqm(i) = Dqm(i)+q_ij*min(0.0_DP,-diff)
        Dqp(j) = Dqp(j)+q_ji*max(0.0_DP, diff)
        Dqm(j) = Dqm(j)+q_ji*min(0.0_DP, diff)
      end do

      ! Apply nodal limiter
      Drp = afcstab_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)

      ! Apply limiter
      do iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
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
    end subroutine doLimit_GP
    
  end subroutine gfsc_buildConvVecGPScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecSymmBlock(rx, dscale, ry, rafcstab)

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

    ! scaling parameter
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry

    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or. ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecSymmBlock')
      call sys_halt()

    else

      call gfsc_buildConvVecSymmScalar(rx%RvectorBlock(1),&
          dscale, ry%RvectorBlock(1), rafcstab)

    end if
  end subroutine gfsc_buildConvVecSymmBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecSymmScalar(rx, dscale, ry, rafcstab)

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

    ! scaling parameter
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry

    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dx,p_Dy,p_Dflux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if ((rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_SYMMETRIC) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0)    .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecSymmScalar')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
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
    call doLimit_Symmetric(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
        p_Dx, dscale, rafcstab%NEDGE, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
        p_Drp, p_Drm, p_Dflux, p_Dy)

    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_NODELIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)
    
  contains
    
    ! Here, the working routine follows
    
    !**************************************************************
    ! Perform symmetric flux limiting
    
    subroutine doLimit_Symmetric(IverticesAtEdge, DcoefficientsAtEdge,&
        Dx, dscale, NEDGE, Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dy)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm,Dflux,Dy

      ! local variables
      real(DP) :: d_ij,f_ij,s_ij,diff
      integer :: iedge,ij,i,j
      
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(Dpp)
      call lalg_clearVectorDble(Dpm)
      call lalg_clearVectorDble(Dqp)
      call lalg_clearVectorDble(Dqm)
      
      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
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
      
      ! Apply the nodal limiter
      Drp = afcstab_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
      
      ! Apply limiter
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
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
    end subroutine doLimit_Symmetric

  end subroutine gfsc_buildConvVecSymmScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianBlock(RcoeffMatrices, rx,&
      fcb_calcMatrix, hstep, bbuildStabilisation, bclear, rjacobianMatrix)

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

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then
      
      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvJacobianScalar(RcoeffMatrices,&
          rx%RvectorBlock(1), fcb_calcMatrix, hstep,&
          bbuildStabilisation, bclear, rjacobianMatrix)
      
    end if

  end subroutine gfsc_buildConvJacobianBlock

   !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianScalar(RcoeffMatrices, rx,&
      fcb_calcMatrix, hstep, bbuildStabilisation, bclear, rjacobianMatrix)

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

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    real(DP), dimension(:), pointer :: p_CoeffX,p_CoeffY,p_CoeffZ,p_Jac,p_Dx
    integer :: h_Ksep,ndim
    
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)
    
    ! Set pointers
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      
    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_CoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianScalar')
      call sys_halt()
    end select
    
    
    ! What kind of matrix are we?
    select case(rjacobianMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bbuildStabilisation) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Jac)
        end select

      else   ! bbuildStabilisation

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Jac)
        end select

      end if   ! bbuildStabilisation

      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bbuildStabilisation) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Jac)
        end select
      
      else   ! bbuildStabilisation

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobianMatrix%NEQ, p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Jac)
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

    subroutine doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, CoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients k_ij and k_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij ,k_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

    subroutine doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, CoeffX, CoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients k_ij and k_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij ,k_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

    subroutine doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, CoeffX, CoeffY, CoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
          C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients k_ij and k_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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
        CoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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
        CoeffX, CoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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
        CoeffX, CoeffY, CoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
          C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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

    subroutine doUpwindMat7_1D(Kld, Kcol, Ksep, NEQ, CoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

    subroutine doUpwindMat7_2D(Kld, Kcol, Ksep, NEQ, CoeffX, CoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

    subroutine doUpwindMat7_3D(Kld, Kcol, Ksep, NEQ, CoeffX, CoeffY, CoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
          C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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
        CoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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
        CoeffX, CoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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
        CoeffX, CoeffY, CoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx
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
          C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
          C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
          C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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

          ! Compute perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
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
      bclear, rafcstab, rjacobianMatrix,rconsistentMassMatrix)

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
    type(t_matrixScalar), intent(in), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearFCTBlock')
      call sys_halt()

    else

      call gfsc_buildJacLinearFCTScalar(rx%RvectorBlock(1),&
          theta, tstep, hstep, bclear, rafcstab,&
          rjacobianMatrix, rconsistentMassMatrix)

    end if
  end subroutine gfsc_buildJacLinearFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearFCTScalar(rx, theta, tstep, hstep,&
      bclear, rafcstab, rjacobianMatrix, rconsistentMassMatrix)

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
    type(t_matrixScalar), intent(in), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dflux,p_Dflux0
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_Dx
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearFCTScalar')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT_IMPLICIT)
      
      ! What kind of matrix are we?
      select case(rjacobianMatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        !-------------------------------------------------------------------------
        ! Matrix format 7
        !-------------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
        
        if (present(rconsistentMassMatrix)) then
          call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
          call doJacobian_implFCTconsMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_MC, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        end if
        
      case(LSYSSC_MATRIX9)
        !-------------------------------------------------------------------------
        ! Matrix format 9
        !-------------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
        call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
        
        if (present(rconsistentMassMatrix)) then
          call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
          call doJacobian_implFCTconsMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_MC, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Dx,&
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

    subroutine doJacobian_implFCTnoMass(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
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

    subroutine doJacobian_implFCTconsMass(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, MC, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
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
      bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

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
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearTVDBlock')
      call sys_halt()

    else

      call gfsc_buildJacLinearTVDScalar(rx%RvectorBlock(1), tstep,&
          hstep, bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacLinearTVDBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearTVDScalar(rx, tstep, hstep,&
      bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

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
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm
    real(DP), dimension(:), pointer :: p_Jac,p_Dx,p_Dflux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    
    integer :: h_Ksep
    logical :: bisExtended


    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)          .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearTVDScalar')
      call sys_halt()
    end if
      
    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_generateOffdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if


    ! What kind of matrix format are we?
    select case(rjacobianMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      
      ! Create diagonal separator and increase it by one
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_TVD(&
          p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
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
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      
      call doJacobianMat79_TVD(&
          p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
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
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
        Dx, Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
              IverticesAtEdge, DcoefficientsAtEdge, Dx,&
              Dpp, Dpm, Dqp, Dqm, tstep, hstep, iedge, iloc, k,&
              Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              IverticesAtEdge, DcoefficientsAtEdge, Dx,&
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
                IverticesAtEdge, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IverticesAtEdge, Kdiagonal, Dflux,&
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
    subroutine updateJacobianMat79_TVD(IverticesAtEdge,&
        DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        iedge, iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)

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
    subroutine assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal,&
        Dflux, Kloc, Drploc, Drmloc, Dfluxloc, hstep, iedge, iloc, k, l,&
        bisExtended, Ksep, Jac)

      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
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

  subroutine gfsc_buildJacLinearGPBlock(rconsistentMassMatrix, rx,&
      rx0, theta, tstep, hstep, bclear, rafcstab, rjacobianMatrix,&
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
    type(t_matrixScalar), intent(in) :: rconsistentMassMatrix

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
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1 .or.&
        rx0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearGPBlock')
      call sys_halt()

    else

      call gfsc_buildJacLinearGPScalar(&
          rconsistentMassMatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), theta, tstep, hstep,&
          bclear, rafcstab, rjacobianMatrix, bextendedSparsity)
      
    end if
  end subroutine gfsc_buildJacLinearGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearGPScalar(rconsistentMassMatrix, rx,&
      rx0, theta, tstep, hstep, bclear, rafcstab, rjacobianMatrix,&
      bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rconsistentMassMatrix

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
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_Dx,p_Dx0,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep
    logical :: bisExtended

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)          .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearGPScalar')
      call sys_halt()
    end if
    
    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_generateOffdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
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
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(rx0, p_Dx0)

    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if


    ! What kind of matrix format are we?
    select case(rjacobianMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_GP(&
          p_IsuperdiagEdgesIdx, p_IverticesAtEdge, p_IsubdiagEdgesIdx,&
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
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      
      call doJacobianMat79_GP(&
          p_IsuperdiagEdgesIdx, p_IverticesAtEdge, p_IsubdiagEdgesIdx,&
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
    subroutine doJacobianMat79_GP(IsuperdiagEdgesIdx, IverticesAtEdge,&
        IsubdiagEdgesIdx, IsubdiagEdges, DcoefficientsAtEdge, Kld,&
        Kcol, Kdiagonal, MC, Dx, Dx0, Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp,&
        Drm, theta, tstep, hstep, NEQ, NEDGE, NNVEDGE, bisExtended,&
        bisMat7, Ksep, Jac)
    
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0,Dflux,Dflux0,Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
              IverticesAtEdge, DcoefficientsAtEdge, MC, Dx, Dx0,&
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
              IverticesAtEdge, DcoefficientsAtEdge, MC, Dx, Dx0,&
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
                IverticesAtEdge, Kdiagonal, Dflux, Dflux0, Drp, Drm,&
                Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IverticesAtEdge, Kdiagonal, Dflux, Dflux0, Drp, Drm,&
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
    subroutine updateJacobianMat79_GP(IverticesAtEdge,&
        DcoefficientsAtEdge, MC, Dx, Dx0, Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm,&
        theta, tstep, hstep, iedge, iloc, k, Dpploc, Dpmloc,&
        Dqploc, Dqmloc, Dfluxloc, Dfluxloc0, Kloc)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0,Dflux,Dflux0,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
      i  = IverticesAtEdge(1,iedge)
      j  = IverticesAtEdge(2,iedge)
      ij = IverticesAtEdge(3,iedge)

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
      diff = tstep*(theta*diff1+(1-theta)*diff0)

      ! Compute antidiffusive flux 
      if (abs(diff) < SYS_EPSREAL) then
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
          diff = tstep*(theta*diff1+(1-theta)*diff0)
          
          ! Compute antidiffusive flux
          if (abs(diff) < SYS_EPSREAL) then
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
          diff = tstep*(theta*diff1+(1-theta)*diff0)
          
          ! Compute antidiffusive flux
          if (abs(diff) < SYS_EPSREAL) then
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
    subroutine assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal,&
        Dflux, Dflux0, Drp, Drm, Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0,&
        hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(:), intent(in) :: Drp,Drm,Dflux,Dflux0
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
      i=IverticesAtEdge(1,iedge)
      j=IverticesAtEdge(2,iedge)
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
      fcb_calcMatrix, theta, tstep, hstep, bclear, rafcstab,&
      rjacobianMatrix, rconsistentMassMatrix)

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
    include 'intf_gfsccallback.inc'

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTBlock')
      call sys_halt()

    else

      call gfsc_buildJacobianFCTScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrix, theta,&
          tstep, hstep, bclear, rafcstab, rjacobianMatrix,&
          rconsistentMassMatrix)

    end if
  end subroutine gfsc_buildJacobianFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianFCTScalar(RcoeffMatrices, rx,&
      fcb_calcMatrix, theta, tstep, hstep, bclear, rafcstab,&
      rjacobianMatrix, rconsistentMassMatrix)

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
    include 'intf_gfsccallback.inc'
    
    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dflux,p_Dflux0,p_Dx,p_CoeffX,p_CoeffY,p_CoeffZ,p_MC,p_Jac
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer :: ndim
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTScalar')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_CoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTScalar')
      call sys_halt()
    end select
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT_IMPLICIT)
      
      ! What kind of matrix format are we?
      select case(rjacobianMatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        !-----------------------------------------------------------------------
        ! Matrix format 7
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
        
        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_1D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld,&
                p_CoeffX, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld,&
                p_CoeffX, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM2D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld,&
                p_CoeffX, p_CoeffY, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE,  p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld,&
                p_CoeffX, p_CoeffY, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM3D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld,&
                p_CoeffX, p_CoeffY, p_CoeffZ, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld,&
                p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
        end select
        
        
      case(LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
        
        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_1D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_CoeffX, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_CoeffX, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
            
        case (NDIM2D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_CoeffX, p_CoeffY, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_CoeffX, p_CoeffY, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM3D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_CoeffX, p_CoeffY, p_CoeffZ, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Dflux, p_Dflux0,&
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
    subroutine doJacobian_implFCTnoMass_1D(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, CoeffX, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
    subroutine doJacobian_implFCTconsMass_1D(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, CoeffX, MC, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j) 
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
    subroutine doJacobian_implFCTnoMass_2D(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, CoeffX, CoeffY, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
    subroutine doJacobian_implFCTconsMass_2D(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, CoeffX, CoeffY, MC, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j) 
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
    subroutine doJacobian_implFCTnoMass_3D(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, CoeffX, CoeffY, CoeffZ, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep


        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
    subroutine doJacobian_implFCTconsMass_3D(IverticesAtEdge,&
        DcoefficientsAtEdge, Kdiagonal, CoeffX, CoeffY, CoeffZ, MC, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine matrix indices
        ij = IverticesAtEdge(3,iedge)
        ji = IverticesAtEdge(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = CoeffX(ij); C_ji(1) = CoeffX(ji)
        C_ij(2) = CoeffY(ij); C_ji(2) = CoeffY(ji)
        C_ij(3) = CoeffZ(ij); C_ji(3) = CoeffZ(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j) 
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
      fcb_calcMatrix, tstep, hstep, bclear, rafcstab,&
      rjacobianMatrix, bextendedSparsity)

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
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then
      
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDBlock')
      call sys_halt()

    else

      call gfsc_buildJacobianTVDScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrix, tstep,&
          hstep, bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianTVDBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianTVDScalar(RcoeffMatrices, rx,&
      fcb_calcMatrix, tstep, hstep, bclear, rafcstab,&
      rjacobianMatrix, bextendedSparsity)

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
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Dflux
    real(DP), dimension(:), pointer :: p_CoeffX,p_CoeffY,p_CoeffZ,p_Jac,p_Dx
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisExtended
    

    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)          .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_CoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDScalar')
      call sys_halt()
    end select

    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_generateOffdiagEdges(rafcstab)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if

    
    ! What kind of matrix format are we?
    select case(rjacobianMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_TVD_1D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_CoeffX, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_CoeffX, p_CoeffY, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Dflux, p_Dpp, p_Dpm,&
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
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix,   p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_TVD_1D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_CoeffX, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_CoeffX, p_CoeffY, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_CoeffX, p_CoeffY, p_CoeffZ, p_Dx, p_Dflux, p_Dpp, p_Dpm,&
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
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
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
    ! Assemble the Jacobian matrix for FEM-TVD in 1D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD_1D(IsuperdiagEdgesIdx,&
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
        CoeffX, Dx, Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = CoeffX(ij)
          c_ji = CoeffX(ji)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = CoeffX(ij)
          c_ji = CoeffX(ji)

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
                IverticesAtEdge, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IverticesAtEdge, Kdiagonal, Dflux,&
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
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, CoeffX, CoeffY, Dx, Dflux,&
        Dpp, Dpm, Dqp, Dqm, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji)/)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji)/)

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
                IverticesAtEdge, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IverticesAtEdge, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
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
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, CoeffX, CoeffY, CoeffZ, Dx,&
        Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij),CoeffZ(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji),CoeffZ(ji)/)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij),CoeffZ(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji),CoeffZ(ji)/)

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
                IverticesAtEdge, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IverticesAtEdge, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
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
        
        ! Compute perturbed coefficients k_ij and k_ji
        call fcb_calcMatrix(Dx(i)+dsign*hstep_ik, Dx(j)+dsign*hstep_jk,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
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
    subroutine assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal,&
        Dflux, Kloc, Drploc, Drmloc, Dfluxloc, hstep, iedge, iloc, k, l,&
        bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge,Kloc
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
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
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

  subroutine gfsc_buildJacobianGPBlock(RcoeffMatrices,&
      rconsistentMassMatrix, rx, rx0, fcb_calcMatrix, theta, tstep,&
      hstep, bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

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
    type(t_matrixScalar), intent(in) :: rconsistentMassMatrix

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
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1 .or.&
        rx0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPBlock')
      call sys_halt()
      
    else
      
      call gfsc_buildJacobianGPScalar(&
          RcoeffMatrices, rconsistentMassMatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), fcb_calcMatrix, theta, tstep, hstep,&
          bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianGPScalar(RcoeffMatrices,&
      rconsistentMassMatrix, rx, rx0, fcb_calcMatrix, theta, tstep,&
      hstep, bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

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
    type(t_matrixScalar), intent(in) :: rconsistentMassMatrix

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
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm,p_Dflux,p_Dflux0
    real(DP), dimension(:), pointer :: p_CoeffX,p_CoeffY,p_CoeffZ,p_MC,p_Jac,p_Dx,p_Dx0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisExtended
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)          .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
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
    call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(rx0, p_Dx0)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_CoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_CoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_CoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPScalar')
      call sys_halt()
    end select

    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_generateOffdiagEdges(rafcstab)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if

    
    ! What kind of matrix format are we?
    select case(rjacobianMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_GP_1D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_CoeffX, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_GP_2D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_CoeffX, p_CoeffY, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_GP_3D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_CoeffX, p_CoeffY, p_CoeffZ, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
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
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
            
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_GP_1D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_CoeffX, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_GP_2D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_CoeffX, p_CoeffY, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_GP_3D(&
            p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_CoeffX, p_CoeffY, p_CoeffZ, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
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
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
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
    ! Assemble the Jacobian matrix for FEM-GP in 1D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP_1D(IsuperdiagEdgesIdx,&
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, CoeffX, MC, Dx, Dx0,&
        Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp, Drm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,MC,Dx,Dx0,Dflux,Dflux0
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = CoeffX(ij)
          c_ji = CoeffX(ji)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = CoeffX(ij)
          c_ji = CoeffX(ji)
   
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
                IverticesAtEdge, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IverticesAtEdge, Kdiagonal, Dflux,&
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
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, CoeffX, CoeffY, MC, Dx, Dx0,&
        Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp, Drm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,MC,Dx,Dx0,Dflux,Dflux0
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji)/)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji)/)
   
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
                IverticesAtEdge, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IverticesAtEdge, Kdiagonal, Dflux,&
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
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, CoeffX, CoeffY, CoeffZ, MC, Dx,&
        Dx0, Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp, Drm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: CoeffX,CoeffY,CoeffZ,MC,Dx,Dx0,Dflux,Dflux0
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij),CoeffZ(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji),CoeffZ(ji)/)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/CoeffX(ij),CoeffY(ij),CoeffZ(ij)/)
          c_ji = (/CoeffX(ji),CoeffY(ji),CoeffZ(ji)/)
   
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
                IverticesAtEdge, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IverticesAtEdge, Kdiagonal, Dflux,&
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
      diff = tstep*(theta*diff1+(1-theta)*diff0)

      ! Compute antidiffusive flux 
      if (abs(diff) < SYS_EPSREAL) then
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
        
        ! Compute perturbed velocity
        call fcb_calcMatrix(Dx(i)+dsign*hstep_ik, Dx(j)+dsign*hstep_jk,&
            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
          diff = tstep*(theta*diff1+(1-theta)*diff0)

          ! Compute antidiffusive flux
          if (abs(diff) < SYS_EPSREAL) then
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
          if (abs(diff) < SYS_EPSREAL) then
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
    subroutine assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal,&
        Dflux, Dflux0, Drp, Drm, Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0,&
        hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(:), intent(in) :: Dflux,Dflux0,Drp,Drm
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge,Kloc
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
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
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
      rafcstab, rjacobianMatrix, bextendedSparsity)

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

    ! scaling parameter
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
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (rx%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianSymmBlock')
      call sys_halt()

    else

      call gfsc_buildJacobianSymmScalar(rx%RvectorBlock(1), dscale,&
          hstep, bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianSymmBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianSymmScalar(rx, dscale, hstep, bclear,&
      rafcstab, rjacobianMatrix, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete diffusion operator for a scalar convection equation.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling parameter
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
    type(t_matrixScalar), intent(inout) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dflux,p_Dx,p_Jac
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep
    logical :: bisExtended

    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADINCREMENTS)  .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_BOUNDS)        .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianSymmScalar')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)
    
    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_generateOffdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
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
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if


    ! What kind of matrix format are we?
    select case(rjacobianMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_Symm(&
          p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
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
      call lsyssc_getbase_Kld(rjacobianMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rjacobianMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_Symm(&
          p_IsuperdiagEdgesIdx, p_IverticesAtEdge,&
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
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
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
    ! Assemble the Jacobian matrix for symmetric flux limiting
    subroutine doJacobianMat79_Symm(IsuperdiagEdgesIdx,&
        IverticesAtEdge, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, Dx, Dflux, Dpp, Dpm,&
        Dqp, Dqm, Drp, Drm, dscale, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: dscale,hstep
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
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
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

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
                IverticesAtEdge, Kld, Kcol, Dflux, Drp, Drm,&
                Kloc, Drploc, Drmloc, Dfluxloc, dscale, hstep,&
                iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
          
          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_Symm(&
                IverticesAtEdge, Kld, Kcol, Dflux, Drp, Drm,&
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
    subroutine assembleJacobianMat79_Symm(IverticesAtEdge, Kdiagonal,&
        Kcol, Dflux, Drp, Drm, Kloc, Drploc, Drmloc, Dfluxloc, dscale,&
        hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), dimension(:), intent(in) :: Drp,Drm,Dflux
      real(DP), intent(in) :: dscale,hstep
      integer, dimension(:,:), intent(in)  :: IverticesAtEdge,Kloc
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
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)
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

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildFluxFCTBlock(rlumpedMassMatrix,&
      rafcstab, rx1, rx2, theta, tstep, dscale, binit,&
      rconsistentMassMatrix)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for FEM-FCT schemes.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

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
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if ((rx1%nblocks .ne. 1) .or. (rx2%nblocks .ne. 1)) then
      
      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTBlock')
      call sys_halt()
      
    else

      call gfsc_buildFluxFCTScalar(rlumpedMassMatrix,&
          rafcstab, rx1%RvectorBlock(1), rx2%RvectorBlock(1),&
          theta, tstep, dscale, binit, rconsistentMassMatrix)

    end if
    
  end subroutine gfsc_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildFluxFCTScalar(rlumpedMassMatrix,&
      rafcstab, rx1, rx2, theta, tstep, dscale, binit,&
      rconsistentMassMatrix)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for FEM-FCT schemes.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

    ! solution vectors
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
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_ML,p_MC,p_Dx1,p_Dx2
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux0,p_Dflux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
      call sys_halt()
    end if
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%iSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
      call sys_halt()
    end if
    
    ! Set pointers
    call lsyssc_getbase_double(rx1, p_Dx1)
    call lsyssc_getbase_double(rx2, p_Dx2)
    call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)


    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)
    
      if (rafcstab%bprelimiting) then
        print *, "Prelimiting not implemented"
        stop
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
      
      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
      
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
        call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

        ! Subtract amount of rejected antidiffusion
        call doCombineFluxes(rafcstab%NEDGE,&
            -1.0_DP, p_Dalpha, p_Dflux, p_Dflux0)
      end if

      ! Do we have to use the consistent mass matrix?
      if (present(rconsistentMassMatrix)) then

        !-----------------------------------------------------------------------
        ! Include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------
  
        ! Set pointer for consistent mass matrix
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)

        ! Do we have to compute the initial fluxes?
        if (binit) then
          call doFluxesConsMass(p_IverticesAtEdge,&
              p_DcoefficientsAtEdge, rafcstab%NEDGE, p_MC, p_Dx1, p_Dx2,&
              -dscale/tstep, (1-theta)*dscale, p_Dflux0)
          call doFluxesNoMass(p_IverticesAtEdge,&
              p_DcoefficientsAtEdge,rafcstab%NEDGE, p_Dx2, dscale, p_Dflux)
          
        else
          call doFluxesConsMass(p_IverticesAtEdge,&
              p_DcoefficientsAtEdge, rafcstab%NEDGE, p_MC, p_Dx1, p_Dx2,&
              dscale/tstep, theta*dscale, p_Dflux)
          call lalg_vectorLinearComb(p_Dflux0, p_Dflux, 1.0_DP, 1.0_DP)
        end if
        
      else

        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------
        
        call doFluxesNoMass(p_IverticesAtEdge,&
            p_DcoefficientsAtEdge, rafcstab%NEDGE, p_Dx2, dscale, p_Dflux)
        
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
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rafcstab%p_rvectorPrelimit, p_Dflux0)
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
        
        ! Set pointers
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
        
        ! There are only initial fluxes
        call doFluxesConsMass(p_IverticesAtEdge,&
            p_DcoefficientsAtEdge, rafcstab%NEDGE, p_MC, p_Dx1, p_Dx2,&
            dscale, dscale, p_Dflux)
      end if

      if (.not.(present(rconsistentMassMatrix)) .or.&
          rafcstab%bprelimiting) then
        
        !-----------------------------------------------------------------------
        ! Do not include contribution of the consistent mass matrix
        !-----------------------------------------------------------------------

        ! Set pointer
        if (rafcstab%bprelimiting) then
          call lsyssc_getbase_double(rafcstab%p_rvectorPrelimit, p_Dflux)
        else
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        end if
        
        ! There are only initial fluxes
        call doFluxesNoMass(p_IverticesAtEdge,&
            p_DcoefficientsAtEdge,rafcstab%NEDGE, p_Dx2, dscale, p_Dflux)
        
        ! Prelimiting is only necessary if the consistent mass matrix
        ! is built into the raw antidiffusive fluxes. However, if the
        ! switch for prelimiting was not set to .false. and no mass
        ! antidiffusion is built into the fluxes we can simply copy it
        if (.not.(present(rconsistentMassMatrix)) .and.&
            rafcstab%bprelimiting) then
          call lsyssc_copyVector(rafcstab%p_rvectorPrelimit,&
              rafcstab%p_rvectorFlux)
        end if
      end if
      
      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_ADFLUXES)
      
      
    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains
    
    ! Here, the working routines follow

    !**************************************************************
    ! Combine two fluxes: flux2 := flux2+dscale*alpha*flux2

    subroutine doCombineFluxes(NEDGE, dscale, Dalpha, Dflux1, Dflux2)
      
      real(DP), dimension(:), intent(in) :: Dflux1,Dalpha
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Dflux2

      ! local variables
      integer :: iedge

      !$omp parallel do
      do iedge = 1, NEDGE
        Dflux2(iedge) = Dflux2(iedge) +&
            dscale * Dalpha(iedge) * Dflux1(iedge)
      end do
      !$omp end parallel do
      
    end subroutine doCombineFluxes

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with
    ! contribution of the consistent mass matrix.

    subroutine doFluxesConsMass(IverticesAtEdge, DcoefficientsAtEdge,&
        NEDGE, MC, Dx1, Dx2, dscale1, dscale2, Dflux)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx1,Dx2
      real(DP), intent(in) :: dscale1,dscale2
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(out) :: Dflux


      ! local variables
      integer :: iedge,ij,i,j
      
      if (dscale1 .eq. 0.0_DP) then
        
        if (dscale2 .eq. 0.0_DP) then
          
          call lalg_clearVector(Dflux)
          
        else
          
          call doFluxesNoMass(IverticesAtEdge,&
              DcoefficientsAtEdge, NEDGE, Dx2, dscale2, Dflux)
          
        end if
        
      else
        
        if (dscale2 .eq. 0.0_DP) then
          
          !$omp parallel do private(i,j)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IverticesAtEdge(1,iedge)
            j  = IverticesAtEdge(2,iedge)
            ij = IverticesAtEdge(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = dscale1*MC(ij)*(Dx1(i)-Dx1(j))
          end do
          !$omp end parallel do

        else

          !$omp parallel do private(i,j)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IverticesAtEdge(1,iedge)
            j  = IverticesAtEdge(2,iedge)
            ij = IverticesAtEdge(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = dscale1*MC(ij)*(Dx1(i)-Dx1(j))+&
                dscale2*DcoefficientsAtEdge(1,iedge)*(Dx2(i)-Dx2(j))
          end do
          !$omp end parallel do
          
        end if

      end if
    end subroutine doFluxesConsMass

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without
    ! contribution of the consistent mass matrix.

    subroutine doFluxesNoMass(IverticesAtEdge, DcoefficientsAtEdge,&
        NEDGE, Dx, dscale, Dflux)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      if (dscale .eq. 0.0_DP) then
        
        call lalg_clearVector(Dflux)

      else

        !$omp parallel do private(i,j)
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          
          ! Compute raw antidiffusive flux
          Dflux(iedge) = dscale*DcoefficientsAtEdge(1,iedge)*(Dx(i)-Dx(j))
        end do
        !$omp end parallel do
        
      end if
    end subroutine doFluxesNoMass
    
  end subroutine gfsc_buildFluxFCTScalar

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<function>

  pure function gfsc_hasOrientation(rafcstab) result(bhasOrientation)

!<description>
    ! This function returns .TRUE. if the given stabilisation structure
    ! requires an oriented edge data structure. Otherwise it returns .FALSE.
!</description>

!<input>
    ! stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<result>
    ! =.TRUE. if edge orientation is required
    ! =.FALSE. otherwise
    logical :: bhasOrientation
!</result>
!</function>

    select case(rafcstab%ctypeAFCstabilisation)
    case(AFCSTAB_FEMTVD,&
         AFCSTAB_FEMGP)
      bhasOrientation = .true.
      
    case DEFAULT
      bhasOrientation = .false.
    end select
  end function gfsc_hasOrientation

end module groupfemscalar

!##############################################################################
!# ****************************************************************************
!# <name> groupfemscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for applying the group-finite
!# element formulation to scalar problems, i.e. conservation laws.
!# The technique was proposed by C.A.J. Fletcher in:
!#
!#     C.A.J. Fletcher, "The group finite element formulation"
!#     Computer Methods in Applied Mechanics and Engineering (ISSN 0045-7825),
!#     vol. 37, April 1983, p. 225-244.
!#
!# The group finite element formulation uses the same basis functions for the
!# unknown solution and the fluxes. This allows for an efficient matrix assemble,
!# whereby the constant coefficient matrices can be assembled once and for all
!# at the beginning of the simulation and each time the grid is modified.
!#
!# Moreover, this module allows to modifying discrete operators by means of
!# the algebraic flux correction (AFC) methodology proposed by Kuzmin, Moeller
!# and Turek in a series of publications. As a starting point for scalar
!# conservation laws, the reader is referred to the book chapter
!#
!#     D. Kuzmin and M. Moeller, "Algebraic flux correction I. Scalar conservation
!#     laws", In: D. Kuzmin et al. (eds), Flux-Corrected Transport: Principles, 
!#     Algorithms, and Applications, Springer, 2005, 155-206.
!#
!# A more detailed description of the algorithms is given in the comments of the
!# subroutine implementing the corresponding discretisation schemes. All methods
!# are based on the stabilisation structure t_afcstab which is defined in the 
!# underlying module "afcstabilisation". The initialisation as a scalar stabilisation
!# structure is done by the routine gfsc_initStabilisation. 
!# 
!# There are three types of routines. The gfsc_buildXXXOperator routines can be
!# used to assemble the discrete convection or diffusion operators resulting
!# from the standard Galerkin finite element discretisation plus some discretely
!# defined artificial diffusion. For convective terms, this technique is termed
!# "discrete upwinding", whereas for physical diffusion operators this approach
!# ensures that the discrete maximum principle holds. Consequently, the term
!# DMP is adopted.
!#
!# The second type of routines is given by gfsc_buildResidualXXX. They can be
!# used to update/initialise the residual term applying some sort of algebraic
!# flux correction. Importantly, the family of AFC schemes gives rise to nonlinear
!# algebraic equations that need to be solved iteratively. Thus, it is usefull to
!# build the compensating antidiffusion into the residual vector rather than the
!# right hand side. However, it is still possible to give a negative scaling factor.
!#
!# The third type of routines is used to assemble the Jacobian matrix for Newton's
!# method. Here, the exact Jacobian matrix is approximated by means of second-order
!# divided differences whereby the "perturbation parameter" is specified by the
!# user. You should be aware of the fact, that in general the employed  flux limiters
!# are not differentiable globally so that the construction of the Jacobian matrix
!# is somehow "delicate". Even though the routines will produce some matrix without
!# warnings, this matrix may be singular and/or ill-conditioned.
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
!# 6.) gfsc_buildResidualFCT = gfsc_buildResScalarFCT /
!#                             gfsc_buildResBlockFCT
!#     -> Assembles the residual vector for AFC stabilisation of FCT type
!#
!# 7.) gfsc_buildResidualTVD = gfsc_buildResScalarTVD /
!#                             gfsc_buildResBlockTVD
!#     -> Assembles the residual vector for AFC stabilisation of TVD type
!#
!# 8.) gfsc_buildResidualGP = gfsc_buildResScalarGP /
!#                            gfsc_buildResBlockGP
!#     -> Assembles the residual vector for AFC stabilisation of general-purpose type
!#
!# 9.) gfsc_buildResidualSymm = gfsc_buildResScalarSymm /
!#                              gfsc_buildResBlockSymm
!#     -> Assembles the residual vector for stabilisation by means of
!#        symmetric flux limiting for diffusion operators
!#
!# 10.) gfsc_buildConvectionJacobian = gfsc_buildConvJacobianScalar /
!#                                     gfsc_buildConvJacobianBlock
!#     -> Assembles the Jacobian matrix for the convective part of
!#        the transport operator for a scalar convection equation
!#
!# 11.) gfsc_buildJacobianFCT = gfsc_buildJacLinearScalarFCT /
!#                              gfsc_buildJacLinearBlockFCT /
!#                              gfsc_buildJacobianScalarFCT /
!#                              gfsc_buildJacobianBlockFCT
!#      -> Assembles the Jacobian matrix for the stabilisation part of FCT type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 12.) gfsc_buildJacobianTVD = gfsc_buildJacLinearScalarTVD /
!#                              gfsc_buildJacLinearBlockTVD /
!#                              gfsc_buildJacobianScalarTVD /
!#                              gfsc_buildJacobianBlockTVD
!#      -> Assembles the Jacobian matrix for the stabilisation part of TVD type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 13.) gfsc_buildJacobianGPD = gfsc_buildJacLinearScalarGP /
!#                              gfsc_buildJacLinearBlockGP /
!#                              gfsc_buildJacobianScalarGP /
!#                              gfsc_buildJacobianBlockGP
!#      -> Assembles the Jacobian matrix for the stabilisation part of general purpose limiter; 
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 14.) gfsc_buildJacobianSymm = gfsc_buildJacobianScalarSymm /
!#                               gfsc_buildJacobianBlockSymm
!#      -> Assembles the Jacobian matrix for the stabilisation part of symmetric
!#         flux limiting for diffusion operators
!#
!# 
!# The following auxiliary routines are available:
!#
!# 1.) gfsc_hasOrientation
!#     -> Checks if the stabilisation techniques requires an oriented structure
!#
!# </purpose>
!##############################################################################

module groupfemscalar

  use afcstabilisation
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use storage

  implicit none

  private
  public :: gfsc_initStabilisation
  public :: gfsc_isMatrixCompatible
  public :: gfsc_isVectorCompatible
  public :: gfsc_buildConvectionOperator
  public :: gfsc_buildDiffusionOperator
  public :: gfsc_buildResidualFCT
  public :: gfsc_buildResidualTVD
  public :: gfsc_buildResidualGP
  public :: gfsc_buildResidualSymm
  public :: gfsc_buildConvectionJacobian
  public :: gfsc_buildJacobianFCT
  public :: gfsc_buildJacobianTVD
  public :: gfsc_buildJacobianGP
  public :: gfsc_buildJacobianSymm

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface gfsc_buildConvectionOperator
    module procedure gfsc_buildConvOperatorScalar
    module procedure gfsc_buildConvOperatorBlock
  end interface

  interface gfsc_buildResidualFCT
    module procedure gfsc_buildResScalarFCT
    module procedure gfsc_buildResBlockFCT
  end interface

  interface gfsc_buildResidualTVD
    module procedure gfsc_buildResScalarTVD
    module procedure gfsc_buildResBlockTVD
  end interface

  interface gfsc_buildResidualGP
    module procedure gfsc_buildResScalarGP
    module procedure gfsc_buildResBlockGP
  end interface

  interface gfsc_buildResidualSymm
    module procedure gfsc_buildResScalarSymm
    module procedure gfsc_buildResBlockSymm
  end interface

  interface gfsc_buildConvectionJacobian
    module procedure gfsc_buildConvJacobianScalar
    module procedure gfsc_buildConvJacobianBlock
  end interface

  interface gfsc_buildJacobianFCT
    module procedure gfsc_buildJacLinearScalarFCT
    module procedure gfsc_buildJacLinearBlockFCT
    module procedure gfsc_buildJacobianScalarFCT
    module procedure gfsc_buildJacobianBlockFCT
  end interface

  interface gfsc_buildJacobianTVD
    module procedure gfsc_buildJacLinearScalarTVD
    module procedure gfsc_buildJacLinearBlockTVD
    module procedure gfsc_buildJacobianScalarTVD
    module procedure gfsc_buildJacobianBlockTVD
  end interface

  interface gfsc_buildJacobianGP
    module procedure gfsc_buildJacLinearScalarGP
    module procedure gfsc_buildJacLinearBlockGP
    module procedure gfsc_buildJacobianScalarGP
    module procedure gfsc_buildJacobianBlockGP
  end interface

  interface gfsc_buildJacobianSymm
    module procedure gfsc_buildJacobianScalarSymm
    module procedure gfsc_buildJacobianBlockSymm
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
    type(t_matrixScalar), intent(IN) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: i

    
    ! Set atomic data
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = int(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ))

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)
      ! do nothing

      
    case (AFCSTAB_FEMFCT,&
          AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMGP,&
          AFCSTAB_FEMTVD)
      
      ! Handle for IsuperdiagonalEdgesIdx
      if (rafcstab%h_IsuperdiagonalEdgesIdx .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IsuperdiagonalEdgesIdx)
      call storage_new('gfsc_initStabilisation', 'IsuperdiagonalEdgesIdx',&
          rafcstab%NEQ+1, ST_INT, rafcstab%h_IsuperdiagonalEdgesIdx, ST_NEWBLOCK_NOINIT)
      
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

      ! We need 6 nodal vectors for P's, Q's and R's
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
                                 rafcstab%NEQ, .false., ST_DOUBLE)
      end do

      ! We need 2 edgewise vectors for the explicit and implicit fluxes
      allocate(rafcstab%RedgeVectors(2))
      do i = 1, 2
        call lsyssc_createVector(rafcstab%RedgeVectors(i),&
                                 rafcstab%NEDGE, .false., ST_DOUBLE)
      end do


    case (AFCSTAB_SYMMETRIC)

      ! Handle for IsuperdiagonalEdgesIdx
      if (rafcstab%h_IsuperdiagonalEdgesIdx .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IsuperdiagonalEdgesIdx)
      call storage_new('gfsc_initStabilisation', 'IsuperdiagonalEdgesIdx',&
          rafcstab%NEQ+1, ST_INT, rafcstab%h_IsuperdiagonalEdgesIdx, ST_NEWBLOCK_NOINIT)
      
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

      ! We need 6 nodal vectors for P's, Q's and R's
      allocate(rafcstab%RnodalVectors(6))
      do i = 1, 6
        call lsyssc_createVector(rafcstab%RnodalVectors(i),&
                                 rafcstab%NEQ, .false., ST_DOUBLE)
      end do
      
      ! We need 1 edgewise vector for the fluxes
      allocate(rafcstab%RedgeVectors(1))
      call lsyssc_createVector(rafcstab%RedgeVectors(1),&
                               rafcstab%NEDGE, .false., ST_DOUBLE)


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
    type(t_matrixScalar), intent(IN) :: rmatrix

    ! The stabilisation structure
    type(t_afcstab), intent(IN)      :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether matrix and stabilisation are compatible or
    ! not.  If not given, an error will inform the user if the
    ! matrix/operator are not compatible and the program will halt.
    logical, intent(OUT), optional :: bcompatible
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
    type(t_vectorScalar), intent(IN) :: rvector

    ! Teh stabilisation structure
    type(t_afcstab), intent(IN)      :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether matrix and stabilisation are compatible or
    ! not. If not given, an error will inform the user if the
    ! matrix/operator are not compatible and the program will halt.
    logical, intent(OUT), optional :: bcompatible
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

  subroutine gfsc_buildConvOperatorBlock(RcoeffMatrices, ru, fcb_calcConvection,&
                                         bStabilise, bclear, rconvMatrix, rafcstab,&
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
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! The solution vector
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorBlock), intent(IN) :: ru
    
    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN) :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for (non-)conservative matrix assembly
    ! TRUE  : assemble conservative convection operator (default)
    ! FALSE : assemble non-conservative convection operator
    logical, intent(IN), optional :: bisConservative

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! The transport operator
    type(t_matrixScalar), intent(INOUT) :: rconvMatrix
    
    ! OPTIONAL: the stabilisation structure
    type(t_afcstab), intent(INOUT), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (ru%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvOperatorScalar(RcoeffMatrices, ru%RvectorBlock(1),&
                                        fcb_calcConvection, bStabilise, bclear,&
                                        rconvMatrix, rafcstab, bisConservative)

    end if

  end subroutine gfsc_buildConvOperatorBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvOperatorScalar(RcoeffMatrices, ru, fcb_calcConvection,&
                                          bStabilise, bclear, rconvMatrix, rafcstab,&
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
    !     D. Kuzmin and M. Moeller, "Algebraic flux correction I. Scalar
    !     conservation laws", In: D. Kuzmin et al. (eds), Flux-Corrected
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
    !     "upwind" and corresponds to the column number whose entry has been
    !     eliminated.
!</description>

!<input>
    ! The array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! The solution vector
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorScalar), intent(IN) :: ru

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN) :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for (non-)conservative matrix assembly
    ! TRUE  : assemble conservative convection operator (default)
    ! FALSE : assemble non-conservative convection operator
    logical, intent(IN), optional :: bisConservative

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! The transport operator
    type(t_matrixScalar), intent(INOUT) :: rconvMatrix
    
    ! OPTIONAL: the stabilisation structure
    type(t_afcstab), intent(INOUT), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Cx, p_Cy, p_Cz, p_ConvOp, p_u
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Ksep, p_Kdiagonal, p_IsuperdiagonalEdgesIdx
    integer :: h_Ksep, ndim
    logical :: bconservative


    ! Check if conservative of non-conservative convection operator is required
    bconservative = .true.
    if (present(bisConservative)) bconservative = bisConservative
        
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rconvMatrix)

    ! Set pointer
    call lsyssc_getbase_double(rconvMatrix, p_ConvOp)
    call lsyssc_getbase_double(ru, p_u)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      call sys_halt()
    end select
    

    ! What kind of matrix are we?
    select case(rconvMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rconvMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rconvMatrix, p_Kcol)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rconvMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rconvMatrix%NEQ+1)

      ! Do we have a stabilisation structure?
      if (present(rafcstab)) then

        ! Check if stabilisation has been initialised
        if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
          call output_line('Stabilisation has not been initialised',&
                           OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
          call sys_halt()
        end if
        
        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        ! Do we need edge orientation?
        if (gfsc_hasOrientation(rafcstab)) then
          
          ! Adopt orientation convention IJ, such that L_ij < L_ji
          ! and generate edge structure for the flux limiter

          if (bconservative) then

            ! Conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doUpwindOAFCMat7Cons1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindOAFCMat7Cons2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindOAFCMat7Cons3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else

            ! Non-conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doUpwindOAFCMat7Nonc1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindOAFCMat7Nonc2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindOAFCMat7Nonc3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION)

        else   ! bhasOrientation == no

          ! Adopt no orientation convention and generate edge structure
          
          if (bconservative) then

            ! Conservative formulation of convection operator
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindAFCMat7Cons1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                         p_Cx, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindAFCMat7Cons2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                         p_Cx, p_Cy, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindAFCMat7Cons3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                         p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else

            ! Non-conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doUpwindAFCMat7Nonc1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                         p_Cx, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindAFCMat7Nonc2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                         p_Cx, p_Cy, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindAFCMat7Nonc3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                         p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

        end if   ! bhasOrientation
        
      elseif (bStabilise) then   ! present(rafcstab) == no
        
        ! Perform discrete upwinding but do not generate the edge data structure

        if (bconservative) then

          ! Conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doUpwindMat7Cons1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindMat7Cons2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindMat7Cons3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        else

          ! Non-conservative formulation of convection operator
          
          select case(ndim)
          case (NDIM1D)
            call doUpwindMat7Nonc1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindMat7Nonc2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindMat7Nonc3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        end if

      else   ! present(rafcstab) == no and bStabilise == no
        
        ! Apply standard Galerkin discretisation
        
        if (bconservative) then

          ! Conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doGalerkinMat7Cons1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinMat7Cons2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinMat7Cons3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        else

          ! Non-conservative formulation of convection operator
          
          select case(ndim)
          case (NDIM1D)
            call doGalerkinMat7Nonc1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinMat7Nonc2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinMat7Nonc3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        end if

      end if   ! present(rafcstab)


      ! Release diagonal separator
      call storage_free(h_Ksep)


    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rconvMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rconvMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rconvMatrix, p_Kdiagonal)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rconvMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rconvMatrix%NEQ+1)

      ! Do we have a stabilisation structure?
      if (present(rafcstab)) then
        
        ! Check if stabilisation has been initialised
        if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
          call output_line('Stabilisation has not been initialised',&
                           OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
          call sys_halt()
        end if
        
        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)

        ! Do we need edge orientation?
        if (gfsc_hasOrientation(rafcstab)) then

          ! Adopt orientation convention IJ, such that L_ij < L_ji 
          ! and generate edge structure for the flux limiter

          if (bconservative) then

            ! Conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doUpwindOAFCMat9Cons1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindOAFCMat9Cons2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindOAFCMat9Cons3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else

            ! Non-conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doUpwindOAFCMat9Nonc1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindOAFCMat9Nonc2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindOAFCMat9Nonc3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION)

        else   ! bhasOrientation == no

          ! Adopt no orientation convention and generate edge structure
          
          if (bconservative) then
            
            ! Conservative formulation of convection operator
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindAFCMat9Cons1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                         rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindAFCMat9Cons2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                         rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindAFCMat9Cons3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                         rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else

            ! Non-conservative formulation of convection operator

            select case(ndim)
            case (NDIM1D)
              call doUpwindAFCMat9Nonc1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                         rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindAFCMat9Nonc2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                         rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindAFCMat9Nonc3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                         rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                         p_IsuperdiagonalEdgesIdx,&
                                         p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

        end if   ! bhasOrientation
          
      elseif (bStabilise) then   ! present(rafcstab) == no
        
        ! Perform discrete upwinding but do not generate the edge data structure
        
        if (bconservative) then
          
          ! Conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doUpwindMat9Cons1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindMat9Cons2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindMat9Cons3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        else

          ! Non-conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doUpwindMat9Nonc1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindMat9Nonc2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindMat9Nonc3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        end if
                  
      else   ! present(rafcstab) == no and bStabilise == no
        
        ! Apply standard Galerkin discretisation
        
        if (bconservative) then
          
          ! Conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doGalerkinMat9Cons1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinMat9Cons2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinMat9Cons3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        else

          ! Non-conservative formulation of convection operator

          select case(ndim)
          case (NDIM1D)
            call doGalerkinMat9Nonc1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinMat9Nonc2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinMat9Nonc3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        end if
        
      end if   ! present(rafcstab)
      
      
      ! Release diagonal separator
      call storage_free(h_Ksep)
       
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routine follow
        
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation
    
    subroutine doGalerkinMat7Cons1D(Kld, Kcol, Ksep, NEQ, Cx, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: k_ii,k_ij,k_ji
      integer :: ii,ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Compute coefficients
        C_ii(1) = Cx(ii)

        ! Compute coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii) + k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
        end do
      end do
    end subroutine doGalerkinMat7Cons1D

    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation
    
    subroutine doGalerkinMat7Nonc1D(Kld, Kcol, Ksep, NEQ, Cx, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
      integer :: ii,jj,ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j); jj = Kld(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ii) = K(ii) - k_ij
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
          K(jj) = K(jj) - k_ji
        end do
      end do
    end subroutine doGalerkinMat7Nonc1D    
    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation
    
    subroutine doGalerkinMat7Cons2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: k_ii,k_ij,k_ji
      integer :: ii,ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
        end do
      end do
    end subroutine doGalerkinMat7Cons2D

    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation
    
    subroutine doGalerkinMat7Nonc2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
      integer :: ii,jj,ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j); jj = Kld(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ii) = K(ii) - k_ij
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
          K(jj) = K(jj) - k_ji
        end do
      end do
    end subroutine doGalerkinMat7Nonc2D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation
    
    subroutine doGalerkinMat7Cons3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP):: k_ii,k_ij,k_ji
      integer :: ii,ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii) + k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
        end do
      end do
    end subroutine doGalerkinMat7Cons3D

    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation
    
    subroutine doGalerkinMat7Nonc3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP):: k_ij,k_ji
      integer :: ii,jj,ij,ji,i,j
      
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j); jj = Kld(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ii) = K(ii) - k_ij
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
          K(jj) = K(jj) - k_ji
        end do
      end do
    end subroutine doGalerkinMat7Nonc3D

    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation
    
    subroutine doGalerkinMat9Cons1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: k_ii,k_ij,k_ji
      integer :: ii,ij,ji,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficient
        C_ii(1) = Cx(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
        end do
      end do
    end subroutine doGalerkinMat9Cons1D
    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation
    
    subroutine doGalerkinMat9Nonc1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
      integer :: ii,jj,ij,ji,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
                
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; jj = Kdiagonal(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ii) = K(ii) - k_ij
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
          K(jj) = K(jj) - k_ji
        end do
      end do
    end subroutine doGalerkinMat9Nonc1D
    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation
    
    subroutine doGalerkinMat9Cons2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: k_ii,k_ij,k_ji
      integer :: ii,ij,ji,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii) + k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
        end do
      end do
    end subroutine doGalerkinMat9Cons2D

    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation
    
    subroutine doGalerkinMat9Nonc2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
      integer :: ii,jj,ij,ji,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; jj = Kdiagonal(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ii) = K(ii) - k_ij
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
          K(jj) = K(jj) - k_ji
        end do
      end do
    end subroutine doGalerkinMat9Nonc2D
    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation
    
    subroutine doGalerkinMat9Cons3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: k_ii,k_ij,k_ji
      integer :: ii,ij,ji,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        K(ii) = K(ii) + k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
        end do
      end do
    end subroutine doGalerkinMat9Cons3D

    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation
    
    subroutine doGalerkinMat9Nonc3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
      integer :: ii,jj,ij,ji,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; jj = Kdiagonal(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Assemble the global operator
          K(ii) = K(ii) - k_ij
          K(ij) = K(ij) + k_ij
          K(ji) = K(ji) + k_ji
          K(jj) = K(jj) - k_ji
        end do
      end do
    end subroutine doGalerkinMat9Nonc3D

    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation
    
    subroutine doUpwindMat7Cons1D(Kld, Kcol, Ksep, NEQ, Cx, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficient
        C_ij(1) = Cx(ij)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
        end do
      end do
    end subroutine doUpwindMat7Cons1D
    
    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation
    
    subroutine doUpwindMat7Nonc1D(Kld, Kcol, Ksep, NEQ, Cx, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
        end do
      end do
    end subroutine doUpwindMat7Nonc1D

    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 7
    ! Conservative treatement
    
    subroutine doUpwindMat7Cons2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
        end do
      end do
    end subroutine doUpwindMat7Cons2D

    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative treatement
    
    subroutine doUpwindMat7Nonc2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
        end do
      end do
    end subroutine doUpwindMat7Nonc2D
    
    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation
    
    subroutine doUpwindMat7Cons3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
        end do
      end do
    end subroutine doUpwindMat7Cons3D

    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation
    
    subroutine doUpwindMat7Nonc3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
        end do
      end do
    end subroutine doUpwindMat7Nonc3D

    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation
    
    subroutine doUpwindMat9Cons1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficient
        C_ij(1) = Cx(ij)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
        end do
      end do
    end subroutine doUpwindMat9Cons1D

    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation
    
    subroutine doUpwindMat9Nonc1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
        end do
      end do
    end subroutine doUpwindMat9Nonc1D
    
    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation
    
    subroutine doUpwindMat9Cons2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP):: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
        end do
      end do
    end subroutine doUpwindMat9Cons2D

    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation
    
    subroutine doUpwindMat9Nonc2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP):: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
        end do
      end do
    end subroutine doUpwindMat9Nonc2D
    
    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation
    
    subroutine doUpwindMat9Cons3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
        end do
      end do
    end subroutine doUpwindMat9Cons3D

    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation
    
    subroutine doUpwindMat9Nonc3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,i,j
      
            
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
        end do
      end do
    end subroutine doUpwindMat9Nonc3D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation
    
    subroutine doUpwindAFCMat7Cons1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,& 
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficient
        C_ii(1) = Cx(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat7Cons1D

    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation
    
    subroutine doUpwindAFCMat7Nonc1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,& 
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat7Nonc1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation

    subroutine doUpwindAFCMat7Cons2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat7Cons2D

    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation

    subroutine doUpwindAFCMat7Nonc2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat7Nonc2D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation

    subroutine doUpwindAFCMat7Cons3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
            
      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat7Cons3D

    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation

    subroutine doUpwindAFCMat7Nonc3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep
            
      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat7Nonc3D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindAFCMat9Cons1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficient
        C_ii(1) = Cx(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat9Cons1D

    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation

    subroutine doUpwindAFCMat9Nonc1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat9Nonc1D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindAFCMat9Cons2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      

      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
                
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat9Cons2D

    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation

    subroutine doUpwindAFCMat9Nonc2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      

      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
                
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat9Nonc2D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge 
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindAFCMat9Cons3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat9Cons3D

    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge 
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation

    subroutine doUpwindAFCMat9Nonc3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                     IsuperdiagonalEdgesIdx,&
                                     IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindAFCMat9Nonc3D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation

    subroutine doUpwindOAFCMat7Cons1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat7Cons1D

    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation

    subroutine doUpwindOAFCMat7Nonc1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          
          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat7Nonc1D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation

    subroutine doUpwindOAFCMat7Cons2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)
        
        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat7Cons2D

    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation

    subroutine doUpwindOAFCMat7Nonc2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat7Nonc2D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7
    ! Conservative formulation

    subroutine doUpwindOAFCMat7Cons3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat7Cons3D

    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7
    ! Non-conservative formulation

    subroutine doUpwindOAFCMat7Nonc3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     =(/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) =(/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat7Nonc3D
        
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindOAFCMat9Cons1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat9Cons1D

    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindOAFCMat9Nonc1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat9Nonc1D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindOAFCMat9Cons2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii)
        
        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii

        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat9Cons2D

    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation

    subroutine doUpwindOAFCMat9Nonc2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat9Nonc2D
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9
    ! Conservative formulation

    subroutine doUpwindOAFCMat9Cons3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ
      
      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ii,C_ij,C_ji
      real(DP) :: d_ij,k_ii,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficients
        C_ii(1) = Cx(ii); C_ii(2) = Cy(ii); C_ii(3) = Cz(ii)

        ! Compute convection coefficients for diagonal
        call fcb_calcConvection(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii)
        
        ! Update the diagonal coefficient
        L(ii) = L(ii) + k_ii
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - d_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat9Cons3D

    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9
    ! Non-conservative formulation

    subroutine doUpwindOAFCMat9Nonc3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
                                      IsuperdiagonalEdgesIdx,&
                                      IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ
      
      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,k_ij,k_ji
      integer :: ii,ij,ji,jj,iedge,i,j
      
      
      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

          ! Compute convection coefficients
          call fcb_calcConvection(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij + d_ij
          k_ji = k_ji + d_ij
          
          ! Assemble the global operator
          L(ii) = L(ii) - k_ij
          L(ij) = L(ij) + k_ij 
          L(ji) = L(ji) + k_ji
          L(jj) = L(jj) - k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC with edge orientation
          if (k_ij < k_ji) then
            IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)
          else
            IverticesAtEdge(:,iedge)     = (/j, i, ji, ij/)
            DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ji, k_ij/)
          end if
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindOAFCMat9Nonc3D

  end subroutine gfsc_buildConvOperatorScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildDiffusionOperator(rcoeffMatrix, bStabilise, bclear,&
                                         rdiffMatrix, rafcstab)

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
    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN) :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear
!</input>

!<inputoutput>
    ! (anisotropic) diffusion operator
    type(t_matrixScalar), intent(INOUT) :: rcoeffMatrix

    ! diffusion operator
    type(t_matrixScalar), intent(INOUT) :: rdiffMatrix

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(INOUT), optional :: rafcstab    
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_S,p_DiffOp
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal,p_IsuperdiagonalEdgesIdx
    integer :: h_Ksep

    
    ! Should matrix be cleared?
    if (bclear) call lsyssc_clearMatrix(rdiffMatrix)

    ! Set pointers
    call lsyssc_getbase_double(rcoeffMatrix, p_S)
    call lsyssc_getbase_double(rdiffMatrix, p_DiffOp)      

    
    ! What kind of matrix are we?
    select case(rdiffMatrix%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rcoeffMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rcoeffMatrix, p_Kcol)

      ! Do we have a stabilisation structure?
      if (present(rafcstab)) then

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rdiffMatrix%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rdiffMatrix%NEQ+1)

        ! Check if stabilisation has been initialised
        if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
          call output_line('Stabilisation has not been initialised',&
                            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffusionOperator')
          call sys_halt()
        end if
        
        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    
        call doLoworderMat7_AFC(p_Kld, p_Kcol, p_Ksep, rdiffMatrix%NEQ,&
                                p_S, p_DiffOp, p_IsuperdiagonalEdgesIdx,&
                                p_IverticesAtEdge, p_DcoefficientsAtEdge)

        ! Release diagonal separator
        call storage_free(h_Ksep)

        ! Set state of stabilisation
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

      elseif (bStabilise) then   ! present(rafcstab)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rdiffMatrix%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rdiffMatrix%NEQ+1)

        call doLoworderMat7(p_Kld, p_Kcol, p_Ksep, rdiffMatrix%NEQ, p_S, p_DiffOp)

        ! Release diagonal separator
        call storage_free(h_Ksep)

      else   ! present(rafcstab) and bStabilise
        
        call lsyssc_duplicateMatrix(rcoeffMatrix, rdiffMatrix,&
                                    LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)

      end if   ! present(rafcstab)
 
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rcoeffMatrix, p_Kld)
      call lsyssc_getbase_Kcol(rcoeffMatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rcoeffMatrix, p_Kdiagonal)

      ! Do we have a stabilisation structure?
      if (present(rafcstab)) then

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rdiffMatrix%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rdiffMatrix%NEQ+1)

        ! Check if stabilisation has been initialised
        if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
          call output_line('Stabilisation has not been initialised',&
                            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffusionOperator')
          call sys_halt()
        end if
        
        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        call doLoworderMat9_AFC(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rdiffMatrix%NEQ,&
                                p_S, p_DiffOp, p_IsuperdiagonalEdgesIdx,&
                                p_IverticesAtEdge, p_DcoefficientsAtEdge)

        ! Release diagonal separator
        call storage_free(h_Ksep)

        ! Set state of stabilisation
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
        
      elseif (bStabilise) then   ! present(rafcstab)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rdiffMatrix%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rdiffMatrix%NEQ+1)

        call doLoworderMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rdiffMatrix%NEQ, p_S, p_DiffOp)

        ! Release diagonal separator
        call storage_free(h_Ksep)

      else   ! present(rafcstab) and bStabilise
        
        call lsyssc_duplicateMatrix(rcoeffMatrix, rdiffMatrix,&
                                    LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPYOVERWRITE)

      end if   ! present(rafcstab)


    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildDiffusionOperator')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble low-order diffusion operator S.
    ! All matrices are stored in matrix format 7
    
    subroutine doLoworderMat7(Kld, Kcol, Ksep, NEQ, S, L)

      real(DP), dimension(:), intent(IN) :: S
      integer, dimension(:), intent(IN) :: Kld, Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP) :: d_ij
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Artificial diffusion coefficient
          d_ij = max(0.0_DP, -S(ij)) 
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doLoworderMat7

    
    !**************************************************************
    ! Assemble low-order diffusion operator S.
    ! All matrices are stored in matrix format 9
    
    subroutine doLoworderMat9(Kld, Kcol, Kdiagonal, Ksep, NEQ, S, L)

      real(DP), dimension(:), intent(IN) :: S
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ
      
      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP) :: d_ij
      integer :: ii,ij,ji,jj,i,j


      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Artificial diffusion coefficient
          d_ij = max(0.0_DP, -S(ij)) 
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doLoworderMat9


    !**************************************************************
    ! Assemble low-order diffusion operator S and AFC data.
    ! All matrices are stored in matrix format 7
    
    subroutine doLoworderMat7_AFC(Kld, Kcol, Ksep, NEQ, S, L,&
                                  IsuperdiagonalEdgesIdx,&
                                  IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: S
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx
      
      ! local variables
      real(DP) :: d_ij,s_ij
      integer :: ii,ij,ji,jj,iedge,i,j
      

      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kld(i)

        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Artificial diffusion coefficient
          d_ij = max(0.0_DP, -S(ij))
          s_ij = max(0.0_DP,  S(ij))
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij

          ! Increase edge counter
          iedge = iedge+1

          ! AFC
          IverticesAtEdge(:,iedge)     = (/i, j/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, s_ij/)
        end do
      end do

      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doLoworderMat7_AFC


    !**************************************************************
    ! Assemble low-order diffusion operator S and AFC data.
    ! All matrices are stored in matrix format 9
    
    subroutine doLoworderMat9_AFC(Kld, Kcol, Kdiagonal, Ksep, NEQ, S, L,&
                                  IsuperdiagonalEdgesIdx,&
                                  IverticesAtEdge, DcoefficientsAtEdge)

      real(DP), dimension(:), intent(IN) :: S
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: L
      integer, dimension(:), intent(INOUT) :: Ksep

      real(DP), dimension(:,:), intent(OUT) :: DcoefficientsAtEdge
      integer, dimension(:,:), intent(OUT) :: IverticesAtEdge
      integer, dimension(:), intent(OUT) :: IsuperdiagonalEdgesIdx

      ! local variable
      real(DP) :: d_ij,s_ij
      integer :: ii,ij,ji,jj,iedge,i,j
      

      ! Initialize edge counter
      iedge = 0

      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Set initial edge number for row I
        IsuperdiagonalEdgesIdx(i) = iedge+1

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Artificial diffusion coefficient
          d_ij = max(0.0_DP, -S(ij))
          s_ij = max(0.0_DP,  S(ij))
          
          ! Assemble the global operator
          L(ii) = L(ii)-d_ij; L(ij) = L(ij)+d_ij
          L(ji) = L(ji)+d_ij; L(jj) = L(jj)-d_ij

          ! Increase edge counter
          iedge = iedge+1

          ! AFC
          IverticesAtEdge(:,iedge)     = (/i, j/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, s_ij/)
        end do
      end do

      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doLoworderMat9_AFC

  end subroutine gfsc_buildDiffusionOperator

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResBlockFCT(rlumpedMassMatrix, ru, theta, tstep, binit,&
                                   rres, rafcstab, rconsistentMassMatrix)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-FCT type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(IN) :: rlumpedMassMatrix

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep
    
    ! Switch for residual
    ! TRUE  : build the initial residual
    ! FALSE : build an intermediate residual
    logical, intent(IN) :: binit

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(IN), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (ru%nblocks .ne. 1 .or. rres%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockFCT')
      call sys_halt()

    else

      call gfsc_buildResScalarFCT(rlumpedMassMatrix, ru%RvectorBlock(1),&
                                  theta, tstep, binit, rres%RvectorBlock(1),&
                                  rafcstab, rconsistentMassMatrix)

    end if

  end subroutine gfsc_buildResBlockFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResScalarFCT(rlumpedMassMatrix, ru, theta, tstep, binit,&
                                    rres, rafcstab, rconsistentMassMatrix)

!<description>

    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-FCT type. The idea of flux corrected
    ! transport can be traced back to the early SHASTA algorithm by
    ! Boris and Bock in the early 1970s. Zalesak suggested a fully
    ! multi-dimensional generalisation of this approach and paved the
    ! way for a large family of FCT algorithms.
    !
    ! This subroutine provides different algorithms:
    !
    ! 1. Semi-explicit FEM-FCT algorithm
    !
    !    This is the classical algorithm which makes use of Zalesak's
    !    flux limiter and recomputes and auxiliary positivity-
    !    preserving solution in each iteration step. 
    !    The details of this method can be found in:
    !
    !    D. Kuzmin and M. Moeller, "Algebraic flux correction I. Scalar
    !    conservation laws", Ergebnisberichte Angew. Math. 249,
    !    University of Dortmund, 2004.
    !
    ! 2. Semi-implicit FEM-FCT algorithm
    !
    !    This is the FCT algorithm that should be used by default. It
    !    is quite efficient since the nodal correction factors are
    !    only computed in the first iteration and used to limit the
    !    antidiffusive flux from the first iteration. This explicit
    !    predictor is used in all subsequent iterations to constrain
    !    the actual target flux.
    !    The details of this method can be found in:
    !
    !    D. Kuzmin and D. Kourounis, "A semi-implicit FEM-FCT
    !    algorithm for efficient treatment of time-dependent
    !    problems", Ergebnisberichte Angew. Math. 302, University of
    !    Dortmund, 2005.
    !
    ! 3. Linearized FEM-FCT algorithm
    !
    !    A new trend in the development of FCT algorithms is to
    !    linearise the raw antidiffusive fluxes about an intermediate
    !    solution computed by a positivity-preserving low-order
    !    scheme. By virtue of this linearisation, the costly
    !    evaluation of correction factors needs to be performed just
    !    once per time step. Furthermore, no questionable
    !    `prelimiting' of antidiffusive fluxes is required, which
    !    eliminates the danger of artificial steepening.
    !    The details of this method can be found in:
    !
    !    D. Kuzmin, "Explicit and implicit FEM-FCT algorithms with
    !    flux linearization", Ergebnisberichte Angew. Math. 358,
    !    University of Dortmund, 2008.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(IN) :: rlumpedMassMatrix

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep
    
    ! Switch for residual
    ! TRUE  : build the initial residual
    ! FALSE : build an intermediate residual
    logical, intent(IN) :: binit

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(IN), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_fluxImpl,p_fluxExpl
    real(DP), dimension(:), pointer :: p_u,p_res,p_uPredictor
    real(DP), dimension(:), pointer :: p_MC,p_ML
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
      call sys_halt()
    end if
 
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_fluxImpl)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_fluxExpl)
    call lsyssc_getbase_double(rlumpedMassMatrix, p_ML)
    call lsyssc_getbase_double(rres, p_res)
    call lsyssc_getbase_double(ru, p_u)


    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT)

      ! Should we build up the initial residual?
      if (binit) then
        
        ! Do we have a fully implicit time discretisation?
        if (theta < 1.0_DP) then

          ! Compute the low-order predictor
          !
          !   $ \tilde u=u^n+(1-\theta)\Delta t M_L^{-1}Lu^n $
          ! 
          ! whereby the residual of the 0-th iteration is assumed to be
          !
          !  $ r^{(0)}=\Delta t Lu^n $

          call lsyssc_invertedDiagMatVec(rlumpedMassMatrix, rres, 1.0_DP-theta,&
                                         rafcstab%RnodalVectors(5))
          call lsyssc_vectorLinearComb(ru, rafcstab%RnodalVectors(5), 1.0_DP, 1.0_DP)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_uPredictor)
          
          ! Initialise the flux limiter
          if (present(rconsistentMassMatrix)) then
            call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
            call doInit_implFCTconsMass(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_ML,&
                p_u, p_uPredictor, theta, tstep, rafcstab%NEDGE,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_fluxExpl)
          else
            call doInit_implFCTnoMass(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_ML,&
                p_u, p_uPredictor, theta, tstep, rafcstab%NEDGE,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_fluxExpl)
          end if
          
        else   ! theta < 1
          
          ! The low-order predictor is simply given by
          !
          !   $ \tilde u=u^n $
          
          ! Initialise the flux limiter with u in lieu of uLow
          if (present(rconsistentMassMatrix)) then
            call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
            call doInit_implFCTconsMass(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_ML,&
                p_u, p_u, theta, tstep, rafcstab%NEDGE,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_fluxExpl)
          else
            call doInit_implFCTnoMass(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_ML,&
                p_u, p_u, theta, tstep, rafcstab%NEDGE,&
                p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_fluxExpl)
          end if
          
        end if   ! theta < 1
        
        ! Set specifier
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)        

      end if   ! binit
      

      ! Check if correction factors and fluxes are available
      if (iand(rafcstab%iSpec, AFCSTAB_LIMITER) .eq. 0 .or.&
          iand(rafcstab%iSpec, AFCSTAB_FLUXES)  .eq. 0) then
        call output_line('Stabilisation does not provide precomputed fluxes &
                         &and/or nodal correction factors',&
                         OU_CLASS_ERROR,OU_MODE_STD, 'gfsc_buildResScalarFCT')
        call sys_halt()
      end if
      
      ! Apply the limited antidiffusion
      if (present(rconsistentMassMatrix)) then
        call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
        call doLimit_implFCTconsMass(&
            p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_u,&
            p_fluxImpl, p_fluxExpl, theta, tstep, rafcstab%NEDGE, p_res)
      else
        call doLimit_implFCTnoMass(&
            p_IverticesAtEdge, p_DcoefficientsAtEdge, p_u,&
            p_fluxImpl, p_fluxExpl, theta, tstep, rafcstab%NEDGE, p_res)
      end if


    case (AFCSTAB_FEMFCT_CLASSICAL)
      
      ! Should we build up the initial residual?
      if (binit) then
        
        ! Initialise the flux limiter
        if (present(rconsistentMassMatrix)) then
          call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
          call doInit_explFCTconsMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC,&
              p_u, theta, tstep, rafcstab%NEDGE, p_fluxExpl)
        else
          call doInit_explFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_u, theta, tstep, rafcstab%NEDGE, p_fluxExpl)
        end if

        ! Set specifier
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES) 

      end if   ! binit

      ! Check if correction factors and fluxes are available
      if (iand(rafcstab%iSpec, AFCSTAB_LIMITER) .eq. 0 .or.&
          iand(rafcstab%iSpec, AFCSTAB_FLUXES)  .eq. 0) then
        call output_line('Stabilisation does not provide precomputed fluxes &
                         &and/or nodal correction factors',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
        call sys_halt()
      end if

      ! Do we have a fully implicit time discretisation?
      if (theta < 1.0_DP) then
        
        ! Compute the low-order predictor
        !
        !   $ \tilde u=u^n+(1-\theta)\Delta t M_L^{-1}Lu^n $
        ! 
        ! whereby the residual of the 0-th iteration is assumed to be
        !
        !   $ r^{(0)}=\Delta tLu^n $
        
        call lsyssc_invertedDiagMatVec(rlumpedMassMatrix, rres, 1.0_DP-theta,&
                                       rafcstab%RnodalVectors(5))
        call lsyssc_vectorLinearComb(ru, rafcstab%RnodalVectors(5), 1.0_DP, 1.0_DP)
        call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_uPredictor)
        
        ! Apply the flux limiter
        if (present(rconsistentMassMatrix)) then
          call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
          call doLimit_explFCTconsMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_ML,&
              p_u, p_uPredictor, p_fluxExpl, theta, tstep, rafcstab%NEDGE,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_res)
        else
          call doLimit_explFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_ML,&
              p_u, p_uPredictor, p_fluxExpl, theta, tstep, rafcstab%NEDGE,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_res)
        end if
        
      else   ! theta < 1
        
        ! The low-order predictor is simply given by
        !
        !   $ \tilde u=u^n $
        
        ! Apply the flux limiter with u in leu of ulow
        if (present(rconsistentMassMatrix)) then
          call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
          call doLimit_explFCTconsMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_MC, p_ML,&
              p_u, p_u, p_fluxExpl, theta, tstep, rafcstab%NEDGE,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_res)
        else
          call doLimit_explFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_ML,&
              p_u, p_u, p_fluxExpl, theta, tstep, rafcstab%NEDGE,&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_fluxImpl, p_res)
        end if

      end if   ! theta < 1


    case (AFCSTAB_FEMFCT_LINEARIZED)
      call output_line('The linearized FEM-FCT algorithm is not implemented yet!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
      call sys_halt()
      

    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
      call sys_halt()
    end select

  contains
      
    ! Here, the working routine follow

    !**************************************************************
    ! Initialisation of the semi-implicit FEM-FCT procedure,
    ! whereby no mass antidiffusion is built into the residual
    
    subroutine doInit_implFCTnoMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                    ML, u, ulow, theta, tstep, NEDGE,&
                                    pp, pm, qp, qm, rp, rm, fluxImpl, fluxExpl)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: ML
      real(DP), dimension(:), intent(IN) :: u,ulow
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT):: pp,pm,qp,qm,rp,rm,fluxImpl,fluxExpl

      ! local variables
      real(DP) :: d_ij,f_ij,diff
      integer :: iedge,ij,i,j
      

      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)
      
      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        
        ! Determine fluxes
        diff = u(i)-u(j); f_ij = tstep*d_ij*diff
        fluxImpl(iedge) = f_ij
          
        ! Sum of positive/negative fluxes
        pp(i) = pp(i)+max(0.0_DP,f_ij); pp(j) = pp(j)+max(0.0_DP,-f_ij)
        pm(i) = pm(i)+min(0.0_DP,f_ij); pm(j) = pm(j)+min(0.0_DP,-f_ij)
        
        ! Upper/lower bounds
        diff = ulow(j)-ulow(i)
        qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
        qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
      end do
      
      ! Adopt the explicit part (if required)
      if (theta < 1.0_DP) then
        call lalg_copyVectorDble(fluxImpl, fluxExpl)
        call lalg_scaleVectorDble(fluxExpl, 1.0_DP-theta)
      else
        call lalg_clearVectorDble(fluxExpl)
      end if
      
      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit( pp, rp, 0.0_DP)
      rm = ML*qm; rm = afcstab_limit( pm, rm, 0.0_DP)
      
      ! Limiting procedure
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        if (fluxImpl(iedge) > 0.0_DP) then
          fluxImpl(iedge) = min(rp(i), rm(j))*fluxImpl(iedge)
        else
          fluxImpl(iedge) = min(rm(i), rp(j))*fluxImpl(iedge)
        end if
      end do

    end subroutine doInit_implFCTnoMass

    !**************************************************************
    ! Initialisation of the semi-implicit FEM-FCT procedure,
    ! whereby consistent mass antidiffusion is built into the residual
    
    subroutine doInit_implFCTconsMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                         MC, ML, u, ulow, theta, tstep, NEDGE,&
                                         pp, pm, qp, qm, rp, rm, fluxImpl, fluxExpl)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,ML
      real(DP), dimension(:), intent(IN) :: u,ulow
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT):: pp,pm,qp,qm,rp,rm,fluxImpl,fluxExpl

      ! local variables
      real(DP) :: d_ij,f_ij,m_ij,diff
      integer :: iedge,ij,i,j
      

      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)
      
      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
        
        ! Determine fluxes
        diff = u(i)-u(j); f_ij = tstep*d_ij*diff
        fluxImpl(iedge) = f_ij; fluxExpl(iedge) = -m_ij*diff
        
        ! Sum of positive/negative fluxes
        pp(i) = pp(i)+max(0.0_DP,f_ij); pp(j) = pp(j)+max(0.0_DP,-f_ij)
        pm(i) = pm(i)+min(0.0_DP,f_ij); pm(j) = pm(j)+min(0.0_DP,-f_ij)
        
        ! Upper/lower bounds
        diff = ulow(j)-ulow(i)
        qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
        qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
      end do
      
      ! Adopt the explicit part (if required)
      if (theta < 1.0_DP) then
        call lalg_vectorLinearCombDble(fluxImpl, fluxExpl, 1.0_DP-theta, 1.0_DP)
      end if
      
      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit( pp, rp, 0.0_DP)
      rm = ML*qm; rm = afcstab_limit( pm, rm, 0.0_DP)
      
      ! Limiting procedure
      do iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        if (fluxImpl(iedge) > 0.0_DP) then
          fluxImpl(iedge) = min(rp(i), rm(j))*fluxImpl(iedge)
        else
          fluxImpl(iedge) = min(rm(i), rp(j))*fluxImpl(iedge)
        end if
      end do

    end subroutine doInit_implFCTconsMass
    

    !**************************************************************
    ! The semi-implicit FEM-FCT limiting procedure,
    ! whereby no mass antidiffusion is built into the residual
    
    subroutine doLimit_implFCTnoMass(IverticesAtEdge, DcoefficientsAtEdge, u,&
                                     fluxImpl, fluxExpl, theta, tstep, NEDGE, res)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,fluxImpl,fluxExpl
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: res
      
      ! local variables
      real(DP) :: d_ij,f_ij
      integer :: iedge,ij,i,j
      
      
      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        
        ! Determine fluxes
        f_ij = (theta*tstep*d_ij)*(u(i)-u(j))+fluxExpl(iedge)
        
        if (f_ij > 0.0_DP) then
          f_ij = min(f_ij,max(fluxImpl(iedge),0.0_DP))
        else
          f_ij = max(f_ij,min(fluxImpl(iedge),0.0_DP))
        end if
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
      
    end subroutine doLimit_implFCTnoMass


    !**************************************************************
    ! The semi-implicit FEM-FCT limiting procedure,
    ! whereby consistent mass antidiffusion is built into the residual
    
    subroutine doLimit_implFCTconsMass(IverticesAtEdge, DcoefficientsAtEdge, MC, u,&
                                          fluxImpl, fluxExpl, theta, tstep, NEDGE, res)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u,fluxImpl,fluxExpl
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: res
      
      ! local variables
      real(DP) :: d_ij,f_ij,m_ij
      integer :: iedge,ij,i,j
      
      
      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
        
        ! Determine fluxes
        f_ij = (m_ij+theta*tstep*d_ij)*(u(i)-u(j))+fluxExpl(iedge)
        
        if (f_ij > 0.0_DP) then
          f_ij = min(f_ij,max(fluxImpl(iedge),0.0_DP))
        else
          f_ij = max(f_ij,min(fluxImpl(iedge),0.0_DP))
        end if
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
      
    end subroutine doLimit_implFCTconsMass

    
    !**************************************************************
    ! Initialisation of the semi-explicit FEM-FCT procedure,
    ! whereby no mass antidiffusion is built into the residual
    
    subroutine doInit_explFCTnoMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                    u, theta, tstep, NEDGE, fluxExpl)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: fluxExpl
      
      ! local variables
      real(DP) :: d_ij,diff
      integer :: iedge,ij,i,j
      
      
      ! Should we use semi-implicit scheme?
      if (theta < 1.0_DP) then
        
        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          
          ! Determine solution difference
          diff = u(i)-u(j)
          
          ! Determine explicit antidiffusive flux
          fluxExpl(iedge) = (1-theta)*tstep*d_ij*diff
        end do
        
      else
        
        ! Initialise explicit fluxes by zero
        call lalg_clearVectorDble(fluxExpl)
        
      end if
      
    end subroutine doInit_explFCTnoMass


    !**************************************************************
    ! Initialisation of the semi-explicit FEM-FCT procedure,
    ! whereby no mass antidiffusion is built into the residual
    
    subroutine doInit_explFCTconsMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                         MC, u, theta, tstep, NEDGE, fluxExpl)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: fluxExpl
      
      ! local variables
      real(DP) :: d_ij,m_ij,diff
      integer :: iedge,ij,i,j
      
      
      ! Should we use semi-implicit scheme?
      if (theta < 1.0_DP) then
        
        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
          
          ! Determine solution difference
          diff = u(i)-u(j)
          
          ! Determine explicit antidiffusive flux
          fluxExpl(iedge) = -m_ij*diff+(1-theta)*tstep*d_ij*diff
        end do
        
      else
        
        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          m_ij = MC(ij)
          
          ! Determine solution difference
          diff = u(i)-u(j)
          
          ! Determine explicit antidiffusive flux
          fluxExpl(iedge) = -m_ij*diff
        end do
        
      end if
      
    end subroutine doInit_explFCTconsMass

    !**************************************************************
    ! The semi-explicit FEM-FCT limiting procedure,
    ! whereby no mass antidiffusion is built into the residual
    
    subroutine doLimit_explFCTnoMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                     ML, u, ulow, fluxExpl, theta, tstep, NEDGE,&
                                     pp, pm, qp, qm, rp, rm, fluxImpl, res)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: ML,u,ulow,fluxExpl
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT) :: pp,pm,qp,qm,rp,rm,fluxImpl,res

      ! local variables
      real(DP) :: diff,d_ij,f_ij
      integer :: iedge,ij,i,j
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)

      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge); diff=u(i)-u(j)
        
        ! Determine antidiffusive flux
        f_ij = fluxExpl(iedge)+theta*tstep*d_ij*diff
        
        ! Determine low-order solution difference
        diff = ulow(j)-ulow(i)
        
        ! Perform prelimiting
        if (f_ij*diff .ge. 0) f_ij = 0.0_DP         
        fluxImpl(iedge) = f_ij
        
        ! Sum of positive/negative fluxes
        pp(i) = pp(i)+max(0.0_DP,f_ij); pp(j) = pp(j)+max(0.0_DP,-f_ij)
        pm(i) = pm(i)+min(0.0_DP,f_ij); pm(j) = pm(j)+min(0.0_DP,-f_ij)
        
        ! Upper/lower bounds
        qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
        qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
      end do
        
      
      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit(pp, rp, 0.0_DP, 1.0_DP)
      rm = ML*qm; rm = afcstab_limit(pm, rm, 0.0_DP, 1.0_DP)
      
      ! Limiting procedure
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        if (fluxImpl(iedge) > 0.0_DP) then
          f_ij = min(rp(i), rm(j))*fluxImpl(iedge)
        else
          f_ij = min(rm(i), rp(j))*fluxImpl(iedge)
        end if

        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do

    end subroutine doLimit_explFCTnoMass


    !**************************************************************
    ! The semi-explicit FEM-FCT limiting procedure,
    ! whereby consistent mass antidiffusion is built into the residual
    
    subroutine doLimit_explFCTconsMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                          MC, ML, u, ulow, fluxExpl, theta, tstep, NEDGE,&
                                          pp, pm, qp, qm, rp, rm, fluxImpl, res)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,ML,u,ulow,fluxExpl
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT) :: pp,pm,qp,qm,rp,rm,fluxImpl,res

      ! local variables
      real(DP) :: diff,d_ij,f_ij,m_ij
      integer :: iedge,ij,i,j
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)

      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IverticesAtEdge(1,iedge)
        j  = IverticesAtEdge(2,iedge)
        ij = IverticesAtEdge(3,iedge)
        
        ! Determine coefficients and solution difference
        d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij); diff=u(i)-u(j)
        
        ! Determine antidiffusive flux
        f_ij = fluxExpl(iedge)+m_ij*diff+theta*tstep*d_ij*diff
        
        ! Determine low-order solution difference
        diff = ulow(j)-ulow(i)
        
        ! Perform prelimiting
        if (f_ij*diff .ge. 0) f_ij = 0.0_DP         
        fluxImpl(iedge) = f_ij
        
        ! Sum of positive/negative fluxes
        pp(i) = pp(i)+max(0.0_DP,f_ij); pp(j) = pp(j)+max(0.0_DP,-f_ij)
        pm(i) = pm(i)+min(0.0_DP,f_ij); pm(j) = pm(j)+min(0.0_DP,-f_ij)
        
        ! Upper/lower bounds
        qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
        qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
      end do
      
      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit(pp, rp, 0.0_DP, 1.0_DP)
      rm = ML*qm; rm = afcstab_limit(pm, rm, 0.0_DP, 1.0_DP)
      
      ! Limiting procedure
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        if (fluxImpl(iedge) > 0.0_DP) then
          f_ij = min(rp(i), rm(j))*fluxImpl(iedge)
        else
          f_ij = min(rm(i), rp(j))*fluxImpl(iedge)
        end if

        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
    end subroutine doLimit_explFCTconsMass

  end subroutine gfsc_buildResScalarFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResBlockTVD(ru, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-TVD type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! time step size
    real(DP), intent(IN) :: tstep
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (ru%nblocks   .ne. 1 .or.&
        rres%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockTVD')
      call sys_halt()

    else

      call gfsc_buildResScalarTVD(ru%RvectorBlock(1), tstep,&
                                  rres%RvectorBlock(1), rafcstab)
      
    end if

  end subroutine gfsc_buildResBlockTVD
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResScalarTVD(ru, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-TVD type.
    !
    ! A detailed description of the FEM-TVD limiter in general is given in:
    !
    !     D. Kuzmin and S. Turek, "Multidimensional FEM-TVD paradigm
    !     for convection-dominated flows" In:  Proceedings of the 
    !     IV European Congress on Computational Methods in Applied Sciences
    !     and Engineering (ECCOMAS 2004). Vol. II, ISBN 951-39-1869-6.
    !
    ! The method actually implemented in this routine is described in:
    !
    !     D. Kuzmin, "Algebraic flux correction for finite element
    !     discretizations of coupled systems" In: E. Onate,
    !     M. Papadrakakis and B. Schrefler (eds.) Computational
    !     Methods for Coupled Problems in Science and Engineering II,
    !     CIMNE, Barcelona, 2007, 653-656.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! time step size
    real(DP), intent(IN) :: tstep
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_u,p_res,p_flux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)  .eq.0 .or. &
        iand(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION).eq.0 .or. &
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)     .eq.0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarTVD')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(rres, p_res)

    ! Perform flux limiting of TVD-type
    call doLimit_TVD(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_u, tstep,&
                     rafcstab%NEDGE, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)
    
  contains

    ! Here, the working routine follows
    
    !**************************************************************
    ! The FEM-TVD limiting procedure
    
    subroutine doLimit_TVD(IverticesAtEdge, DcoefficientsAtEdge, u, tstep,&
                           NEDGE, pp, pm, qp, qm, rp, rm, flux, res)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u
      real(DP), intent(IN) :: tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT) :: pp,pm,qp,qm,rp,rm,flux,res

      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff
      integer :: iedge,ij,i,j
      
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)

      ! Assemble P's and Q's
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        l_ij = DcoefficientsAtEdge(2,iedge)
        l_ji = DcoefficientsAtEdge(3,iedge)
        
        ! Determine solution difference
        diff = tstep*(u(i)-u(j))
        
        ! Prelimit the antidiffusive flux F'_IJ=MIN(-P_IJ,L_JI)(U_I-U_J)
        f_ij = min(d_ij,l_ji)*diff; flux(iedge) = f_ij
        
        ! Assemble P's accordingly
        pp(i) = pp(i)+max(0.0_DP,f_ij); pm(i) = pm(i)+min(0.0_DP,f_ij)
        
        ! Assemble Q's
        qp(i) = qp(i)+max(0.0_DP,-f_ij); qp(j) = qp(j)+max(0.0_DP, f_ij)
        qm(i) = qm(i)+min(0.0_DP,-f_ij); qm(j) = qm(j)+min(0.0_DP, f_ij)
      end do
      
      ! Apply the nodal limiter
      rp = afcstab_limit(pp, qp, 0.0_DP, 1.0_DP)
      rm = afcstab_limit(pm, qm, 0.0_DP, 1.0_DP)

      ! Apply limiter
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed raw antidiffusive flux
        f_ij = flux(iedge)
        
        ! Apply correction factor and store limite flux
        f_ij = merge(rp(i), rm(i), f_ij > 0)*f_ij
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
    end subroutine doLimit_TVD
    
  end subroutine gfsc_buildResScalarTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResBlockGP(rconsistentMassMatrix, ru, ru0, theta, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-GP type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (ru%nblocks   .ne. 1 .or.&
        ru0%nblocks  .ne. 1 .or.&
        rres%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockGP')
      call sys_halt()

    else

      call gfsc_buildResScalarGP(rconsistentMassMatrix, ru%RvectorBlock(1),&
                                 ru0%RvectorBlock(1), theta, tstep,&
                                 rres%RvectorBlock(1), rafcstab)
      
    end if
  end subroutine gfsc_buildResBlockGP

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResScalarGP(rconsistentMassMatrix, ru, ru0, theta, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation using the general purpose limiter.
    !
    ! A detailed description of the FEM-GP limiter in general is given in:
    !
    !     D. Kuzmin, "On the design of general-purpose flux 
    !     limiters for implicit FEM with a consistent mass matrix.
    !     I. Scalar convection."
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
    type(t_matrixScalar), intent(IN) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorScalar), intent(IN) :: ru0

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_MC,p_u,p_u0,p_res,p_fluxImpl,p_fluxExpl
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)  .eq.0 .or. &
        iand(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION).eq.0 .or. &
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)     .eq.0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarGP')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_fluxImpl)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_fluxExpl)
    call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(ru0, p_u0)
    call lsyssc_getbase_double(rres, p_res)

    ! Perform flux limiting by the general purpose limiter
    call doLimit_GP(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
                    p_MC, p_u, p_u0, theta, tstep, rafcstab%NEDGE,&
                    p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
                    p_fluxImpl, p_fluxExpl, p_res)
    
    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)
    
  contains

    ! Here, the working routine follows
    
    !**************************************************************
    ! The FEM-GP limiting procedure
    
    subroutine doLimit_GP(IverticesAtEdge, DcoefficientsAtEdge, MC, u, u0,&
                          theta, tstep, NEDGE, pp, pm, qp, qm, rp, rm,&
                          fluxImpl, fluxExpl, res)


      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u,u0
      real(DP), intent(IN) :: theta,tstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(:), intent(INOUT) :: fluxImpl,fluxExpl,res
      
      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,m_ij,p_ij,pf_ij,df_ij,q_ij,q_ji
      real(DP) :: diff,diff0,diff1
      integer :: iedge,ij,i,j
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)
      
      ! Assemble P' and Q's
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
        
        ! Compute: diff1 = dt*theta*(u_i-u_j) + dt*(1-theta)*(u0_i-u0_j)
        diff1 = u(i)-u(j); diff0 = u0(i)-u0(j)
        diff  = tstep*(theta*diff1+(1.0_DP-theta)*diff0)
        
        ! Compute antidiffusive flux f_ij=min(0,p_ij)*(u_j-u_i)
        if (abs(diff) < SYS_EPSREAL) then
          p_ij = 0
          f_ij = 0
        else
          p_ij = max(0.0_DP,m_ij*(diff1-diff0)/diff+d_ij)
          f_ij = p_ij*diff
        end if
        
        ! Prelimit the antidiffusive flux F'_IJ=MIN(-P_IJ,L_JI)(U_I-U_J)
        pf_ij = min(p_ij,l_ji)*diff; fluxExpl(iedge) = pf_ij
        
        ! Compute the remaining flux dF_IJ=F_IJ-F'_IJ
        df_ij = f_ij-pf_ij; fluxImpl(iedge) = df_ij
        
        ! Assemble P's accordingly
        pp(i) = pp(i)+max(0.0_DP,  f_ij); pm(i) = pm(i)+min(0.0_DP,  f_ij)
        pp(j) = pp(j)+max(0.0_DP,-df_ij); pm(j) = pm(j)+min(0.0_DP,-df_ij)
        
        q_ij = m_ij/tstep+l_ij; q_ji = m_ij/tstep+l_ji

        ! Assemble Q's
        qp(i) = qp(i)+q_ij*max(0.0_DP,-diff); qm(i) = qm(i)+q_ij*min(0.0_DP,-diff)
        qp(j) = qp(j)+q_ji*max(0.0_DP, diff); qm(j) = qm(j)+q_ji*min(0.0_DP, diff)
      end do

      ! Apply nodal limiter
      rp = afcstab_limit(pp, qp, 0.0_DP, 1.0_DP)
      rm = afcstab_limit(pm, qm, 0.0_DP, 1.0_DP)

      ! Apply limiter
      do iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed fluxes
        pf_ij = fluxExpl(iedge); df_ij = fluxImpl(iedge)
        
        ! Limit upwind contribution
        if (pf_ij > 0.0_DP) then
          pf_ij = rp(i)*pf_ij
        else
          pf_ij = rm(i)*pf_ij
        end if
        
        ! Limit symmetric contribution
        if (df_ij > 0.0_DP) then
          df_ij = min(rp(i), rm(j))*df_ij
        else
          df_ij = min(rm(i), rp(j))*df_ij
        end if
        
        f_ij = pf_ij+df_ij
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
    end subroutine doLimit_GP
    
  end subroutine gfsc_buildResScalarGP

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResBlockSymm(ru, dscale, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation by means of symmetric flux limiting for diffusion
    ! operators.  Note that this routine serves as a wrapper for block
    ! vectors. If there is only one block, then the corresponding
    ! scalar routine is called.  Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (ru%nblocks .ne. 1 .or. rres%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockSymm')
      call sys_halt()

    else

      call gfsc_buildResScalarSymm(ru%RvectorBlock(1), dscale,&
                                   rres%RvectorBlock(1), rafcstab)

    end if
  end subroutine gfsc_buildResBlockSymm

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResScalarSymm(ru, dscale, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies stabilisation
    ! by means of symmetric flux limiting for diffusion operators.
    !
    ! Yet, there is no publication available. This routine is based on
    ! private communication with D. Kuzmin.
    !
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! scaling parameter
    real(DP), intent(IN) :: dscale
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_u,p_res,p_flux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    
    ! Check if stabilisation is prepared
    if (rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_SYMMETRIC .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0    .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarSymm')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(rres, p_res)
    
    ! Perform symmetric flux limiting
    call doLimit_Symmetric(p_IverticesAtEdge, p_DcoefficientsAtEdge, p_u, dscale,&
                           rafcstab%NEDGE, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)
    
  contains
    
    ! Here, the working routine follows
    
    !**************************************************************
    ! Perform symmetric flux limiting
    
    subroutine doLimit_Symmetric(IverticesAtEdge, DcoefficientsAtEdge, u, dscale,&
                                 NEDGE, pp, pm, qp, qm, rp, rm, flux, res)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u
      real(DP), intent(IN) :: dscale
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT) :: pp,pm,qp,qm,rp,rm,flux,res

      ! local variables
      real(DP) :: d_ij,f_ij,s_ij,diff
      integer :: iedge,ij,i,j
      
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)
      
      ! Loop over edges
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        s_ij = DcoefficientsAtEdge(2,iedge)
        
        ! Determine fluxes
        diff = u(i)-u(j); f_ij = d_ij*diff
        flux(iedge) = f_ij
        
        ! Sums of raw positive/negative fluxes
        pp(i) = pp(i)+max(0.0_DP, f_ij); pp(j) = pp(j)+max(0.0_DP,-f_ij)
        pm(i) = pm(i)+min(0.0_DP, f_ij); pm(j) = pm(j)+min(0.0_DP,-f_ij)
        
        ! Upper/lower bounds
        f_ij = -s_ij*diff
        qp(i) = qp(i)+max(0.0_DP, f_ij); qp(j) = qp(j)+max(0.0_DP,-f_ij)
        qm(i) = qm(i)+min(0.0_DP, f_ij); qm(j) = qm(j)+min(0.0_DP,-f_ij)
      end do
      
      ! Apply the nodal limiter
      rp = afcstab_limit(pp, qp, 0.0_DP, 1.0_DP)
      rm = afcstab_limit(pm, qm, 0.0_DP, 1.0_DP)
      
      ! Apply limiter
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed raw antidiffusive flux
        f_ij = flux(iedge)
        
        if (f_ij > 0.0_DP) then
          f_ij = dscale*min(rp(i), rm(j))*f_ij
        else
          f_ij = dscale*min(rm(i), rp(j))*f_ij
        end if
        
        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
    end subroutine doLimit_Symmetric

  end subroutine gfsc_buildResScalarSymm

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianBlock(RcoeffMatrices, ru, fcb_calcConvection,&
                                         hstep, bStabilise, bclear, rjacobianMatrix)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective
    ! part of the discrete transport operator for a scalar convection
    ! equation.  Note that this routine serves as a wrapper for block
    ! vectors. If there is only one block, then the corresponding
    ! scalar routine is called.  Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru
    
    ! perturbation parameter
    real(DP), intent(IN) :: hstep

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN) :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (ru%nblocks .ne. 1) then
      
      call output_line('Solution vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvJacobianScalar(RcoeffMatrices, ru%RvectorBlock(1),&
                                        fcb_calcConvection, hstep, bStabilise,&
                                        bclear, rjacobianMatrix)
      
    end if

  end subroutine gfsc_buildConvJacobianBlock

   !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianScalar(RcoeffMatrices, ru, fcb_calcConvection,&
                                          hstep, bStabilise, bclear, rjacobianMatrix)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru
    
    ! perturbation parameter
    real(DP), intent(IN) :: hstep

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN) :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_Jac,p_u
    integer :: h_Ksep,ndim
    
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)
    
    ! Set pointers
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      
    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

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
      if (bStabilise) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep, rjacobianMatrix%NEQ,&
                               p_Cx, p_u, p_Jac)
        case (NDIM2D)
          call doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep, rjacobianMatrix%NEQ,&
                               p_Cx, p_Cy, p_u, p_Jac)
        case (NDIM3D)
          call doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep, rjacobianMatrix%NEQ,&
                               p_Cx, p_Cy, p_Cz, p_u, p_Jac)
        end select

      else   ! bStabilise

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep, rjacobianMatrix%NEQ,&
                                 p_Cx, p_u, p_Jac)
        case (NDIM2D)
          call doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep, rjacobianMatrix%NEQ,&
                                 p_Cx, p_Cy, p_u, p_Jac)
        case (NDIM3D)
          call doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep, rjacobianMatrix%NEQ,&
                                 p_Cx, p_Cy, p_Cz, p_u, p_Jac)
        end select

      end if   ! bStabilise

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
      if (bStabilise) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rjacobianMatrix%NEQ, p_Cx, p_u, p_Jac)
        case (NDIM2D)
          call doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rjacobianMatrix%NEQ, p_Cx, p_Cy, p_u, p_Jac)
        case (NDIM3D)
          call doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rjacobianMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_Jac)
        end select
      
      else   ! bStabilise

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rjacobianMatrix%NEQ, p_Cx, p_u, p_Jac)
        case (NDIM2D)
          call doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rjacobianMatrix%NEQ, p_Cx, p_Cy, p_u, p_Jac)
        case (NDIM3D)
          call doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rjacobianMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_Jac)
        end select

      end if   ! bStabilise

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

    subroutine doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
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
                    
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, k_ij ,k_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
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

    subroutine doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
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
                    
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, k_ij ,k_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
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

    subroutine doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables     
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
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
                    
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
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

    subroutine doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
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

    subroutine doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
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

    subroutine doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, k_ij, k_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, k_ij, k_ji)
          
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
    ! Assemble upwind Jacobian matrix for convective operator in 1D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

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
                    
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
    ! Assemble upwind Jacobian matrix for convective operator in 2D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

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
                    
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
    ! Assemble Jacobian matrix for convective operator in 3D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
    ! Assemble Jacobian matrix for convective operator in 1D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
    ! Assemble Jacobian matrix for convective operator in 2D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
    ! Assemble Jacobian matrix for convective operator in 3D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, Jac)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep

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
          
          ! Compute solution difference u_j-u_i
          diff = u(j)-u(i)

          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

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
          call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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
          call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0.0_DP, -l_ji)
          
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

  subroutine gfsc_buildJacLinearBlockFCT(ru, theta, tstep, hstep, bclear, rafcstab,&
                                         rjacobianMatrix,rconsistentMassMatrix)

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
    type(t_vectorBlock), intent(IN) :: ru

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(IN), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearBlockFCT')
      call sys_halt()

    else

      call gfsc_buildJacLinearScalarFCT(ru%RvectorBlock(1), theta, tstep, hstep, bclear,&
                                        rafcstab, rjacobianMatrix,rconsistentMassMatrix)

    end if
  end subroutine gfsc_buildJacLinearBlockFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearScalarFCT(ru, theta, tstep, hstep, bclear,&
                                          rafcstab, rjacobianMatrix, rconsistentMassMatrix)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  Note that the velocity is assumed
    ! to be linear.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(IN), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_fluxImpl,p_fluxExpl
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_u
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)        .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarFCT')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_fluxImpl)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_fluxExpl)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)
    

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT)
      
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
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_MC, p_u,&
              p_fluxImpl, p_fluxExpl, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_u,&
              p_fluxImpl, p_fluxExpl, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
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
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_MC, p_u,&
              p_fluxImpl, p_fluxExpl, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_u,&
              p_fluxImpl, p_fluxExpl, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        end if

        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarFCT')
        call sys_halt()
      end select
      
    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarFCT')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for semi-implicit FEM-FCT,
    ! whereby no mass antidiffusion is built into the matrix

    subroutine doJacobian_implFCTnoMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                        Kdiagonal, u, fluxImpl, fluxExpl,&
                                        theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,fluxImpl,fluxExpl
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac
      
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
        a_ij = theta*tstep*d_ij
        
        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Compute perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        ! Compute limited antidiffusive flux f(u_ij+h*e_i) 
        f_i = a_ij*diff_i+fluxExpl(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(fluxImpl(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(fluxImpl(iedge), 0.0_DP))
        end if
        
        ! Compute limited antidiffusive flux f(u_ij-h*e_j)
        f_j = a_ij*diff_j+fluxExpl(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(fluxImpl(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(fluxImpl(iedge), 0.0_DP))
        end if
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
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

    subroutine doJacobian_implFCTconsMass(IverticesAtEdge, DcoefficientsAtEdge,&
                                             Kdiagonal, MC, u, fluxImpl, fluxExpl,&
                                             theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u,fluxImpl,fluxExpl
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac
      
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
        a_ij = MC(ij)+theta*tstep*d_ij
        
        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Compute perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        ! Compute limited antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+fluxExpl(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(fluxImpl(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(fluxImpl(iedge), 0.0_DP))
        end if
        
        ! Compute limited antidiffusive flux f(u_ij-h*e_j)
        f_j = a_ij*diff_j+fluxExpl(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(fluxImpl(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(fluxImpl(iedge), 0.0_DP))
        end if
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)+f_ij
        Jac(jj) = Jac(jj)-f_ij
      end do

    end subroutine doJacobian_implFCTconsMass

  end subroutine gfsc_buildJacLinearScalarFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearBlockTVD(ru, tstep, hstep, bclear, rafcstab,&
                                         rjacobianMatrix, bextendedSparsity)

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
    type(t_vectorBlock), intent(IN) :: ru

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearBlockTVD')
      call sys_halt()

    else

      call gfsc_buildJacLinearScalarTVD(ru%RvectorBlock(1), tstep, hstep,&
                                        bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacLinearBlockTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearScalarTVD(ru, tstep, hstep, bclear, rafcstab,&
                                          rjacobianMatrix, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm
    real(DP), dimension(:), pointer :: p_Jac,p_u,p_flux
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagonalEdges
    integer, dimension(:), pointer :: p_IsubdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    
    integer :: h_Ksep
    logical :: bisExtended


    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_BOUNDS)          .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)          .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarTVD')
      call sys_halt()
    end if
      
    ! Check if subdiagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
        call afcstab_generateSubdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagonalEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)
    
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
      
      call doJacobianMat79_TVD(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                               p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                               p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                               p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
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
      
      call doJacobianMat79_TVD(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                               p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                               p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                               p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                               tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                               rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
        
      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarTVD')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    pure subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    pure subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k
      
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    pure subroutine doJacobianMat79_TVD(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                        IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                        u, flux, pp, pm, qp, qm, tstep, hstep,&
                                        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,flux,pp,pm,qp,qm
      real(DP), intent(IN) :: tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc
      integer, dimension(NNVEDGE) :: Kloc
      integer :: k,l,ild,iedge,iloc,nloc
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ

        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)

          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(IverticesAtEdge, DcoefficientsAtEdge, u,&
                                       pp, pm, qp, qm, tstep, hstep, iedge, iloc, k,&
                                       pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          call updateJacobianMat79_TVD(IverticesAtEdge, DcoefficientsAtEdge, u,&
                                       pp, pm, qp, qm, tstep, hstep, iedge, iloc, k,&
                                       pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(iloc)

          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
                                           hstep, iedge, iloc, k, l,&
                                           bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
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
    pure subroutine updateJacobianMat79_TVD(IverticesAtEdge, DcoefficientsAtEdge,&
                                            u, pp, pm, qp, qm, tstep, hstep,&
                                            iedge, iloc, k, pploc, pmloc,&
                                            qploc, qmloc, fluxloc, Kloc)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,pp,pm,qp,qm
      real(DP), intent(IN) :: tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: iedge,k,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT) :: pploc,pmloc,qploc,qmloc,fluxloc
      integer, dimension(:), intent(INOUT) :: Kloc

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
      diff = tstep*(u(i)-u(j))
      f_ij = min(d_ij, l_ji)*diff

      !-------------------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
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
        pploc(:,iloc) = pp(j)
        pmloc(:,iloc) = pm(j)
        qploc(:,iloc) = qp(j)-max(0.0_DP, f_ij)
        qmloc(:,iloc) = qm(j)-min(0.0_DP, f_ij)

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Compute perturbed antidiffusive flux
          f_ij = min(d_ij,l_ji)*(diff+tstep*dsign*hstep)
          fluxloc(iperturb,iloc) = f_ij

          ! For node k which is the upwind node
          pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP, f_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP, f_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,-f_ij)
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,-f_ij)
          
          ! For node l opposite to k which is the downwind node
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,f_ij)
        end do

      else

        ! Store global node number of the opposite node
        Kloc(iloc) = i
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0.0_DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0.0_DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0.0_DP,-f_ij)
        qmloc(:,iloc) = qm(i)-min(0.0_DP,-f_ij)

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Compute perturbed antidiffusive flux
          f_ij = min(d_ij,l_ji)*(diff-tstep*dsign*hstep)
          fluxloc(iperturb,iloc) = f_ij
          
          ! For node k which is the downwind node
          qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,f_ij)
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,f_ij)
          
          ! For node l opposite to k
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,-f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,-f_ij)
        end do
      end if
    end subroutine updateJacobianMat79_TVD


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    pure subroutine assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal,&
                                              flux, Kloc, rploc, rmloc, fluxloc,&
                                              hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:), intent(IN) :: flux
      real(DP), dimension(:,0:), intent(IN) :: rploc,rmloc,fluxloc
      real(DP), intent(IN) :: hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: Kdiagonal,Kloc
      integer, intent(IN) :: iedge,k,l,iloc
      logical, intent(IN) :: bisExtended

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      

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
            f_ij = fluxloc(iperturb,iloc)

            ! Limit flux 
            if (f_ij > 0.0_DP) then
              f_ij = rploc(iperturb,0)*f_ij
            else
              f_ij = rmloc(iperturb,0)*f_ij
            end if

            ! Adopt sign for perturbation direction
            df_ij = df_ij-(iperturb-1.5_DP)*f_ij/hstep
          end do
          
          ! Get corresponding matrix indices
          ik = Kdiagonal(i); jk = Ksep(j)
        else

          do iperturb = 1, 2
            
            ! Retrieve precomputed flux
            f_ij = fluxloc(iperturb,iloc)

            ! Limit flux
            if (f_ij > 0.0_DP) then
              f_ij = rploc(iperturb,iloc)*f_ij
            else
              f_ij = rmloc(iperturb,iloc)*f_ij
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

          if (flux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*(rploc(1,iloc)-rploc(2,iloc))*flux(iedge)/hstep
          else
            f_ij = 0.5_DP*(rmloc(1,iloc)-rmloc(2,iloc))*flux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end if
      end if
    end subroutine assembleJacobianMat79_TVD

  end subroutine gfsc_buildJacLinearScalarTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearBlockGP(rconsistentMassMatrix, ru, ru0, theta,&
                                        tstep, hstep, bclear, rafcstab,&
                                        rjacobianMatrix, bextendedSparsity)

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
    type(t_matrixScalar), intent(IN) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1 .or.&
        ru0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearBlockGP')
      call sys_halt()

    else

      call gfsc_buildJacLinearScalarGP(rconsistentMassMatrix, ru%RvectorBlock(1),&
                                       ru0%RvectorBlock(1), theta, tstep, hstep,&
                                       bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacLinearBlockGP

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearScalarGP(rconsistentMassMatrix, ru, ru0, theta,&
                                         tstep, hstep, bclear, rafcstab,&
                                         rjacobianMatrix, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorScalar), intent(IN) :: ru0

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_u,p_u0,p_flux,p_flux0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagonalEdges
    integer, dimension(:), pointer :: p_IsubdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep
    logical :: bisExtended

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_BOUNDS)          .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)          .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarGP')
      call sys_halt()
    end if
    
    ! Check if subdiagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
        call afcstab_generateSubdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagonalEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(ru0, p_u0)

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
      
      call doJacobianMat79_GP(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                              p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                              p_MC, p_u, p_u0, p_flux, p_flux0,&
                              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
                              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
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
      
      call doJacobianMat79_GP(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                              p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                              p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                              p_MC, p_u, p_u0, p_flux, p_flux0,&
                              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
                              theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                              rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      
      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarGP')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    pure subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    pure subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k
      
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    subroutine doJacobianMat79_GP(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                  IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                  DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                  MC, u, u0, flux, flux0, pp, pm, qp, qm, rp, rm,&
                                  theta, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
                                  bisExtended, bisMat7, Ksep, Jac)
    
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u,u0,flux,flux0,pp,pm,qp,qm,rp,rm
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc,fluxloc0
      integer, dimension(NNVEDGE) :: Kloc
      integer :: k,l,ild,iedge,iloc,nloc
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          call updateJacobianMat79_GP(IverticesAtEdge, DcoefficientsAtEdge, MC, u, u0,&
                                      flux, flux0, pp, pm, qp, qm, theta, tstep,&
                                      hstep, iedge, iloc, k, pploc, pmloc, qploc,&
                                      qmloc, fluxloc, fluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          call updateJacobianMat79_GP(IverticesAtEdge, DcoefficientsAtEdge, MC, u, u0,&
                                      flux, flux0, pp, pm, qp, qm, theta, tstep,&
                                      hstep, iedge, iloc, k, pploc, pmloc, qploc,&
                                      qmloc, fluxloc, fluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux, flux0,&
                                          rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                                          hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux, flux0,&
                                          rp, rm, Kloc, rploc, rmloc, fluxloc, fluxloc0,&
                                          hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
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
    subroutine updateJacobianMat79_GP(IverticesAtEdge, DcoefficientsAtEdge,&
                                      MC, u, u0, flux, flux0, pp, pm, qp, qm,&
                                      theta, tstep, hstep, iedge, iloc, k,&
                                      pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u,u0,flux,flux0,pp,pm,qp,qm
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, intent(IN) :: iedge,k,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT) :: pploc,pmloc,qploc,qmloc,fluxloc,fluxloc0
      integer, dimension(:), intent(INOUT)    :: Kloc

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
      diff1 = u(i)-u(j)
      diff0 = u0(i)-u0(j)

      ! Determine total solution difference
      diff = tstep*(theta*diff1+(1-theta)*diff0)

      ! Compute antidiffusive flux 
      if(abs(diff) < SYS_EPSREAL) then
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
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
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
        pploc(:,iloc) = pp(j)-max(0.0_DP,-df_ij)
        pmloc(:,iloc) = pm(j)-min(0.0_DP,-df_ij)
        qploc(:,iloc) = qp(j)-max(0.0_DP, diff)*q_ji
        qmloc(:,iloc) = qm(j)-min(0.0_DP, diff)*q_ji

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Update solution difference
          diff1 = u(i)-u(j)+dsign*hstep
          
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
          fluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij
          
          ! For node k which is the upwind node
          pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP, f_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP, f_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,-diff)*q_ij
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,-diff)*q_ij
          
          ! For node l opposite to k which is the downwind node
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP,-df_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP,-df_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP, diff)*q_ji
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP, diff)*q_ji
        end do

      else

        ! Store global node number of the opposite node
        Kloc(iloc) = i

        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0.0_DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0.0_DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0.0_DP,-diff)*q_ij
        qmloc(:,iloc) = qm(i)-min(0.0_DP,-diff)*q_ij

        do iperturb = 1, 2
        
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Update solution difference
          diff1 = u(i)-u(j)-dsign*hstep
          
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
          fluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij

          ! For node k which is the downwind node
          pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP,-df_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP,-df_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP, diff)*q_ji
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP, diff)*q_ji
          
          ! For node l opposite to k
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,-diff)*q_ij
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,-diff)*q_ij
        end do
      end if
      
    end subroutine updateJacobianMat79_GP
    

    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal,&
                                        flux, flux0, rp, rm, Kloc,&
                                        rploc, rmloc, fluxloc, fluxloc0,&
                                        hstep, iedge, iloc, k, l,&
                                        bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(IN) :: rploc,rmloc,fluxloc,fluxloc0
      real(DP), dimension(:), intent(IN) :: rp,rm,flux,flux0
      real(DP), intent(IN) :: hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: Kdiagonal,Kloc
      integer, intent(IN) :: iedge,iloc,k,l
      logical, intent(IN) :: bisExtended

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
          df_ij = fluxloc(iperturb,iloc)
          pf_ij = fluxloc0(iperturb,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = rploc(iperturb,0)*pf_ij
            else
              pf_ij = rmloc(iperturb,0)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
            else
              df_ij = min(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = rploc(iperturb,iloc)*pf_ij
            else
              pf_ij = rmloc(iperturb,iloc)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            else
              df_ij = min(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
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
          pf_ij = flux0(iedge)
          df_ij = flux(iedge)

          ! Limit upwind contribution
          if (pf_ij > 0.0_DP) then
            pf_ij = (rploc(1,iloc)-rploc(2,iloc))*pf_ij
          else
            pf_ij = (rmloc(1,iloc)-rmloc(2,iloc))*pf_ij
          end if

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(rploc(1,iloc), rm(j))-&
                     min(rploc(2,iloc), rm(j)))*df_ij
          else
            df_ij = (min(rmloc(1,iloc), rp(j))-&
                     min(rmloc(2,iloc), rp(j)))*df_ij
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
          df_ij = flux(iedge)

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(rp(i), rmloc(1,iloc))-&
                     min(rp(i), rmloc(2,iloc)))*df_ij
          else
            df_ij = (min(rm(i), rploc(1,iloc))-&
                     min(rm(i), rploc(2,iloc)))*df_ij
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
  end subroutine gfsc_buildJacLinearScalarGP

  !*****************************************************************************
  
!<subroutine>

  subroutine gfsc_buildJacobianBlockFCT(RcoeffMatrices, ru, fcb_calcConvection,&
                                        theta, tstep, hstep, bclear, rafcstab,&
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
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices
    
    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(IN), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockFCT')
      call sys_halt()

    else

      call gfsc_buildJacobianScalarFCT(RcoeffMatrices, ru%RvectorBlock(1),&
                                       fcb_calcConvection, theta, tstep, hstep,&
                                       bclear, rafcstab, rjacobianMatrix,&
                                       rconsistentMassMatrix)

    end if
  end subroutine gfsc_buildJacobianBlockFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarFCT(RcoeffMatrices, ru, fcb_calcConvection,&
                                         theta, tstep, hstep, bclear, rafcstab,&
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
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear
    
    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
    
    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(IN), optional :: rconsistentMassMatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_flux,p_flux0,p_u,p_Cx,p_Cy,p_Cz,p_MC,p_Jac
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer :: ndim
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)        .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarFCT')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarFCT')
      call sys_halt()
    end select
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT)
      
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
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Cx,&
                p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Cx,&
                p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM2D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Cx, p_Cy,&
                p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE,  p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Cx, p_Cy,&
                p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM3D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Cx, p_Cy, p_Cz,&
                p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kld, p_Cx, p_Cy, p_Cz,&
                p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
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
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Cx,&
                p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Cx,&
                p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
            
        case (NDIM2D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Cx, p_Cy,&
                p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Cx, p_Cy,&
                p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM3D)
          if (present(rconsistentMassMatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Cx, p_Cy, p_Cz,&
                p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IverticesAtEdge, p_DcoefficientsAtEdge, p_Kdiagonal, p_Cx, p_Cy, p_Cz,&
                p_u, p_flux, p_flux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
        end select
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarFCT')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarFCT')
      call sys_halt()
    end select

  contains
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 1D,
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_1D(IverticesAtEdge, DcoefficientsAtEdge,&
                                           Kdiagonal, Cx, u, flux, flux0,&
                                           theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,u,flux,flux0
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac

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
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        
        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_j+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
        f_i = a_ij*diff_j+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_i+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTnoMass_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 1D,
    ! whereby consistent mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTconsMass_1D(IverticesAtEdge, DcoefficientsAtEdge,&
                                                Kdiagonal, Cx, MC, u, flux, flux0,&
                                                theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,MC,u,flux,flux0
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac

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
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        
        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_j+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
        f_i = a_ij*diff_j+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_i+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTconsMass_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 2D,
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_2D(IverticesAtEdge, DcoefficientsAtEdge,&
                                           Kdiagonal, Cx, Cy, u, flux, flux0,&
                                           theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,u,flux,flux0
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN)   :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac

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
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
        
        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_j+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
        f_i = a_ij*diff_j+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_i+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTnoMass_2D

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 2D,
    ! whereby consistent mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTconsMass_2D(IverticesAtEdge, DcoefficientsAtEdge,&
                                             Kdiagonal, Cx, Cy, MC, u, flux, flux0,&
                                             theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,MC,u,flux,flux0
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN)   :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac

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
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_j+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
        f_i = a_ij*diff_j+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_i+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do

    end subroutine doJacobian_implFCTconsMass_2D
    
    
    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 3D
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_3D(IverticesAtEdge, DcoefficientsAtEdge,&
                                           Kdiagonal, Cx, Cy, Cz, u, flux, flux0,&
                                           theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u,flux,flux0
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: NEDGE

      real(DP), dimension(:), intent(INOUT) :: Jac

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
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
        C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)

        ! Compute solution difference
        diff = u(i)-u(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep


        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if


        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)

        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_j+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if


        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep

        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij

        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------

        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)

        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
        f_i = a_ij*diff_j+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_i+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTnoMass_3D
    

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 3D
    ! whereby consistent mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTconsMass_3D(IverticesAtEdge, DcoefficientsAtEdge,&
                                                Kdiagonal, Cx, Cy, Cz, MC, u, flux, flux0,&
                                                theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,MC,u,flux,flux0
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: Jac

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
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
        C_ij(3) = Cz(ij); C_ji(3) = Cz(ji)
        
        ! Compute solution difference
        diff = u(i)-u(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)+hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
        f_i = a_ij*diff_i+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_j+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)+hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
        f_i = a_ij*diff_j+flux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(flux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute perturbed velocity
        call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = MC(ij)+theta*tstep*max(-l_ij, 0.0_DP, -l_ji)
        
        ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
        f_j = b_ij*diff_i+flux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(flux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(flux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
    
    end subroutine doJacobian_implFCTconsMass_3D

  end subroutine gfsc_buildJacobianScalarFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianBlockTVD(RcoeffMatrices, ru, fcb_calcConvection,&
                                        tstep, hstep, bclear, rafcstab,&
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
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
    
    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then
      
      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockTVD')
      call sys_halt()

    else

      call gfsc_buildJacobianScalarTVD(RcoeffMatrices, ru%RvectorBlock(1),&
                                       fcb_calcConvection, tstep, hstep,&
                                       bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianBlockTVD
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarTVD(RcoeffMatrices, ru, fcb_calcConvection,&
                                         tstep, hstep, bclear, rafcstab,&
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
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear
    
    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
    
    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_flux
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_Jac,p_u
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagonalEdges
    integer, dimension(:), pointer :: p_IsubdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisExtended
    

    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_BOUNDS)          .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)          .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarTVD')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagonalEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarTVD')
      call sys_halt()
    end select

    ! Check if subdiagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
        call afcstab_generateSubdiagEdges(rafcstab)
    
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
        call doJacobianMat79_TVD_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                    p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                    p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                    p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                    tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                    rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                    p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                    p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                    p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                    tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                    rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                    p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                    p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                    p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
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
      call lsyssc_getbase_Kcol(rjacobianMatrix,   p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobianMatrix, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobianMatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobianMatrix%NEQ+1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_TVD_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                    p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                    p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                    p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                    tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                    rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                    p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                    p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                    p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                    tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                    rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                    p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                    p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                    p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                    tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                    rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      end select
        
      ! Free storage
      call storage_free(h_Ksep)

    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarTVD')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    subroutine doJacobianMat79_TVD_1D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                      IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                      DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                      Cx, u, flux, pp, pm, qp, qm, tstep, hstep,&
                                      NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,u,flux,pp,pm,qp,qm
      real(DP), intent(IN) :: tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc
      real(DP), dimension(NDIM1D) :: c_ij, c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc

      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                       c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
                                       iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)

          ! Update local coefficients
          call updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                       c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
                                       iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
                                           hstep, iedge, iloc, k, l,&
                                           bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
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
    subroutine doJacobianMat79_TVD_2D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                      IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                      DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                      Cx, Cy, u, flux, pp, pm, qp, qm, tstep, hstep,&
                                      NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,u,flux,pp,pm,qp,qm
      real(DP), intent(IN) :: tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc
      real(DP), dimension(NDIM2D) :: c_ij, c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc      
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                       c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
                                       iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)

          ! Update local coefficients
          call updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                       c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
                                       iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
                                           hstep, iedge, iloc, k, l,&
                                           bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
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
    end subroutine doJacobianMat79_TVD_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 3D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD_3D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                      IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                      DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                      Cx, Cy, Cz, u, flux, pp, pm, qp, qm, tstep, hstep,&
                                      NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u,flux,pp,pm,qp,qm
      real(DP), intent(IN) :: tstep,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc
      real(DP), dimension(NDIM3D) :: c_ij, c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc


      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                       c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
                                       iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)

          ! Update local coefficients
          call updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                       c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
                                       iloc, k, pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
                                           hstep, iedge, iloc, k, l,&
                                           bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                           Kloc, rploc, rmloc, fluxloc,&
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
    end subroutine doJacobianMat79_TVD_3D


    !**************************************************************
    ! Update the local coefficients for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.    
    subroutine updateJacobianMat79_TVD(DcoefficientsAtEdge, u, pp, pm, qp, qm, &
                                       c_ij, c_ji, tstep, hstep, iedge,&
                                       i, j, ij, ji, iloc, k,&
                                       pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,pp,pm,qp,qm,C_ij,C_ji
      real(DP), intent(IN) :: tstep,hstep
      integer, intent(IN) :: iedge,i,j,k,ij,ji,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT) :: pploc,pmloc,qploc,qmloc,fluxloc
      integer, dimension(:,:), intent(INOUT) :: Kloc

      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff,hstep_ik,hstep_jk,dsign
      integer  :: iperturb
      

      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Determine prelimited antidiffusive flux
      diff = tstep*(u(i)-u(j))
      f_ij = min(d_ij,l_ji)*diff
      
      if (i .eq. k) then
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0.0_DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)
        pmloc(:,iloc) = pm(j)
        qploc(:,iloc) = qp(j)-max(0.0_DP, f_ij)
        qmloc(:,iloc) = qm(j)-min(0.0_DP, f_ij)

      else

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0.0_DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0.0_DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0.0_DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0.0_DP,-f_ij)
        qmloc(:,iloc) = qm(i)-min(0.0_DP,-f_ij)
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
        call fcb_calcConvection(u(i)+dsign*hstep_ik, u(j)+dsign*hstep_jk,&
                                C_ij, C_ji, i, j, l_ij, l_ji)
        
        ! Compute diffusion coefficient
        d_ij = max(-l_ij, 0.0_DP, -l_ji)

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
          fluxloc(iperturb,iloc) = f_ij
        
          if (i .eq. k) then

            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,-f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,-f_ij)
            
            ! For node l opposite to k which is the downwind node
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,f_ij)

          else

            ! For node k which is the downwind node
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,f_ij)
            
            ! For node l opposite to k
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,-f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,-f_ij)

          end if
          
        else
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! In this case the orientation of edge ij needs to be
          ! reverted so as to let i denote the 'upwind' node
          f_ij = -min(d_ij,l_ij)*(diff+tstep*dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij
          
          if (j .eq. k) then

            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,-f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,-f_ij)
            
            ! For node "l" opposite to k which is the downwind node
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,f_ij)

          else

            ! For node k which is the downwind node
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,f_ij)
            
            ! For node "l" opposite to k which is the upwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,-f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,-f_ij)
          end if
        end if
      end do
    end subroutine updateJacobianMat79_TVD


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal, flux,&
                                         Kloc, rploc, rmloc, fluxloc, hstep,&
                                         iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(IN) :: rploc,rmloc,fluxloc
      real(DP), dimension(:), intent(IN) :: flux
      real(DP), intent(IN) :: hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge,Kloc
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: iedge,iloc,k,l
      logical, intent(IN) :: bisExtended

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
          f_ij = fluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit flux 
            if (f_ij > 0.0_DP) then
              f_ij = rploc(iperturb,0)*f_ij
            else
              f_ij = rmloc(iperturb,0)*f_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            if (f_ij > 0.0_DP) then
              f_ij = rploc(iperturb,iloc)*f_ij
            else
              f_ij = rmloc(iperturb,iloc)*f_ij
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

          if (flux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*(rploc(1,iloc)-rploc(2,iloc))*flux(iedge)/hstep
          else
            f_ij = 0.5_DP*(rmloc(1,iloc)-rmloc(2,iloc))*flux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end if
      end if
    end subroutine assembleJacobianMat79_TVD
    
  end subroutine gfsc_buildJacobianScalarTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianBlockGP(RcoeffMatrices, rconsistentMassMatrix,&
                                       ru, ru0, fcb_calcConvection, theta, tstep,&
                                       hstep, bclear, rafcstab, rjacobianMatrix,&
                                       bextendedSparsity)

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
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! consistent mass matrix
    type(t_matrixScalar), intent(IN) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorBlock), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1 .or.&
        ru0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockGP')
      call sys_halt()
      
    else
      
      call gfsc_buildJacobianScalarGP(RcoeffMatrices, rconsistentMassMatrix,&
                                      ru%RvectorBlock(1), ru0%RvectorBlock(1),&
                                      fcb_calcConvection, theta, tstep, hstep,&
                                      bclear, rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianBlockGP

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarGP(RcoeffMatrices, rconsistentMassMatrix,&
                                        ru, ru0, fcb_calcConvection, theta, tstep,&
                                        hstep, bclear, rafcstab, rjacobianMatrix,&
                                        bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary. 
    ! This routine will also work for linear velocities but then it is inefficient
    ! since the solution perturbation does not affect the velocity.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! consistent mass matrix
    type(t_matrixScalar), intent(IN) :: rconsistentMassMatrix

    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! initial solution vector
    type(t_vectorScalar), intent(IN) :: ru0

    ! implicitness parameter
    real(DP), intent(IN) :: theta

    ! time step size
    real(DP), intent(IN) :: tstep

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm,p_flux,p_flux0
    real(DP), dimension(:), pointer :: p_Cx,p_Cy,p_Cz,p_MC,p_Jac,p_u,p_u0
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagonalEdges
    integer, dimension(:), pointer :: p_IsubdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisExtended
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_BOUNDS)          .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)          .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                        OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarGP')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagonalEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    call lsyssc_getbase_double(rconsistentMassMatrix, p_MC)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)
    call lsyssc_getbase_double(ru0, p_u0)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarGP')
      call sys_halt()
    end select

    ! Check if subdiagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
        call afcstab_generateSubdiagEdges(rafcstab)
    
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
        call doJacobianMat79_GP_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                   p_Cx, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                   p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta,&
                                   tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_GP_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                   p_Cx, p_Cy, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                   p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta,&
                                   tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_GP_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                   p_Cx, p_Cy, p_Cz, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                   p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta,&
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
        call doJacobianMat79_GP_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                   p_Cx, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                   p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta,&
                                   tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_GP_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                   p_Cx, p_Cy, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                   p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta,&
                                   tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_GP_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                   p_Cx, p_Cy, p_Cz, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                   p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta,&
                                   tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      end select

      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarGP')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    subroutine doJacobianMat79_GP_1D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                     IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                     DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                     Cx, MC, u, u0, flux, flux0, pp, pm,&
                                     qp, qm, rp, rm, theta, tstep, hstep,&
                                     NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,MC,u,u0,flux,flux0,pp,pm,qp,qm,rp,rm
      real(DP), intent(IN) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc,fluxloc0
      real(DP), dimension(NDIM1D) :: c_ij,c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)
          
          ! Update local coefficients
          call updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                      flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                      theta, tstep, hstep, iedge, i, j, ij, ji,&
                                      iloc, k,  pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = Cx(ij)
          c_ji = Cx(ji)
   
          ! Update local coefficients
          call updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                      flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                      theta, tstep, hstep, iedge, i, j, ij, ji,&
                                      iloc, k, pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                          flux0, rp, rm, Kloc, rploc, rmloc,&
                                          fluxloc, fluxloc0, hstep, iedge,&
                                          iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                          flux0, rp, rm, Kloc, rploc, rmloc,&
                                          fluxloc, fluxloc0, hstep, iedge,&
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
    subroutine doJacobianMat79_GP_2D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                     IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                     DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                     Cx, Cy, MC, u, u0, flux, flux0, pp, pm,&
                                     qp, qm, rp, rm, theta, tstep, hstep,&
                                     NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,MC,u,u0,flux,flux0,pp,pm,qp,qm,rp,rm
      real(DP), intent(IN) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc,fluxloc0
      real(DP), dimension(NDIM2D) :: c_ij,c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc

      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                      flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                      theta, tstep, hstep, iedge, i, j, ij, ji,&
                                      iloc, k, pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij)/)
          c_ji = (/Cx(ji),Cy(ji)/)
   
          ! Update local coefficients
          call updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                      flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                      theta, tstep, hstep, iedge, i, j, ij, ji,&
                                      iloc, k, pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                          flux0, rp, rm, Kloc, rploc, rmloc,&
                                          fluxloc, fluxloc0, hstep, iedge,&
                                          iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                          flux0, rp, rm, Kloc, rploc, rmloc,&
                                          fluxloc, fluxloc0, hstep, iedge,&
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
    subroutine doJacobianMat79_GP_3D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                     IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                     DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                     Cx, Cy, Cz, MC, u, u0, flux, flux0, pp, pm,&
                                     qp, qm, rp, rm, theta, tstep, hstep,&
                                     NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,MC,u,u0,flux,flux0,pp,pm,qp,qm,rp,rm
      real(DP), intent(IN) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc,fluxloc0
      real(DP), dimension(NDIM3D) :: c_ij,c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc


      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                      flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                      theta, tstep, hstep, iedge, i, j, ij, ji,&
                                      iloc, k, pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Determine matrix indices
          ij = IverticesAtEdge(3,iedge)
          ji = IverticesAtEdge(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/Cx(ij),Cy(ij),Cz(ij)/)
          c_ji = (/Cx(ji),Cy(ji),Cz(ji)/)
   
          ! Update local coefficients
          call updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                      flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                      theta, tstep, hstep, iedge, i, j, ij, ji,&
                                      iloc, k, pploc, pmloc, qploc, qmloc,&
                                      fluxloc, fluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                          flux0, rp, rm, Kloc, rploc, rmloc,&
                                          fluxloc, fluxloc0, hstep, iedge,&
                                          iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                          flux0, rp, rm, Kloc, rploc, rmloc,&
                                          fluxloc, fluxloc0, hstep, iedge,&
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
    pure subroutine updateJacobianMat79_GP(DcoefficientsAtEdge, MC, u, u0, flux,&
                                           flux0, pp, pm, qp, qm, c_ij, c_ji,&
                                           theta, tstep, hstep, iedge, i, j, ij, ji,&
                                           iloc, k, pploc, pmloc, qploc, qmloc,&
                                           fluxloc, fluxloc0, Kloc)
      
      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: MC,u,u0,flux,flux0,pp,pm,qp,qm,C_ij,C_ji
      real(DP), intent(IN) :: theta,tstep,hstep
      integer, intent(IN) :: iedge,i,j,k,ij,ji,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT) :: pploc,pmloc,qploc,qmloc,fluxloc,fluxloc0
      integer, dimension(:,:), intent(INOUT)  :: Kloc

      ! local variables
      real(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji
      real(DP) :: diff,diff1,diff0,hstep_ik,hstep_jk,dsign
      integer :: iperturb

      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
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
      diff1 = u(i)-u(j)
      diff0 = u0(i)-u0(j)

      ! Determine total solution difference
      diff = tstep*(theta*diff1+(1-theta)*diff0)

      ! Compute antidiffusive flux 
      if(abs(diff) < SYS_EPSREAL) then
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
        pploc(:,iloc) = pp(j)-max(0.0_DP,-df_ij)
        pmloc(:,iloc) = pm(j)-min(0.0_DP,-df_ij)
        qploc(:,iloc) = qp(j)-max(0.0_DP, diff)*q_ji
        qmloc(:,iloc) = qm(j)-min(0.0_DP, diff)*q_ji

      else

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0.0_DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0.0_DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0.0_DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0.0_DP,-diff)*q_ij
        qmloc(:,iloc) = qm(i)-min(0.0_DP,-diff)*q_ij
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
        call fcb_calcConvection(u(i)+dsign*hstep_ik,&
            u(j)+dsign*hstep_jk, C_ij, C_ji, i, j, l_ij, l_ji)

        ! Compute diffusion coefficient
        d_ij = max(-l_ij, 0.0_DP, -l_ji)

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
          diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)

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
          fluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij
        
          if (i .eq. k) then

            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,-diff)*q_ij
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,-diff)*q_ij
            
            ! For node l opposite to k which is the downwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP,-df_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP,-df_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP, diff)*q_ji
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP, diff)*q_ji

          else

            ! For node k which is the downwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP,-df_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP,-df_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP, diff)*q_ji
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP, diff)*q_ji
            
            ! For node l opposite to k
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,-diff)*q_ij
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,-diff)*q_ij

          end if
          
        else
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! Update solution difference
          diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)
          
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
          fluxloc0(iperturb,iloc) = pf_ij

          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij
          
          if (j .eq. k) then
            
            ! For node k which is the upwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP, diff)*q_ij
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP, diff)*q_ij
                       
            ! For node "l" opposite to k which is the downwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP,-df_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP,-df_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP,-diff)*q_ji
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP,-diff)*q_ji

          else

            ! For node k which is the downwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0.0_DP,-df_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0.0_DP,-df_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0.0_DP,-diff)*q_ji
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0.0_DP,-diff)*q_ji
            
            ! For node "l" opposite to k which is the upwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP, diff)*q_ij
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP, diff)*q_ij

          end if
        end if
      end do
    end subroutine updateJacobianMat79_GP


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    pure subroutine assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal, flux,&
                                             flux0, rp, rm, Kloc, rploc, rmloc,&
                                             fluxloc, fluxloc0, hstep, iedge,&
                                             iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(IN) :: rploc,rmloc,fluxloc,fluxloc0
      real(DP), dimension(:), intent(IN) :: flux,flux0,rp,rm
      real(DP), intent(IN) :: hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge,Kloc
      integer, dimension(:), intent(IN) :: Kdiagonal
      integer, intent(IN) :: iedge,iloc,k,l
      logical, intent(IN) :: bisExtended

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
          df_ij = fluxloc(iperturb,iloc)
          pf_ij = fluxloc0(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = rploc(iperturb,0)*pf_ij
            else
              pf_ij = rmloc(iperturb,0)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
            else
              df_ij = min(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = rploc(iperturb,iloc)*pf_ij
            else
              pf_ij = rmloc(iperturb,iloc)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(rmloc(iperturb,0), rploc(iperturb,iloc))*df_ij
            else
              df_ij = min(rploc(iperturb,0), rmloc(iperturb,iloc))*df_ij
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
          pf_ij = flux0(iedge)
          df_ij = flux(iedge)

          ! Limit upwind contribution
          if (pf_ij > 0.0_DP) then
            pf_ij = (rploc(1,iloc)-rploc(2,iloc))*pf_ij
          else
            pf_ij = (rmloc(1,iloc)-rmloc(2,iloc))*pf_ij
          end if

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(rploc(1,iloc), rm(j))-&
                     min(rploc(2,iloc), rm(j)))*df_ij
          else
            df_ij = (min(rmloc(1,iloc), rp(j))-&
                     min(rmloc(2,iloc), rp(j)))*df_ij
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
          df_ij = flux(iedge)

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(rp(i), rmloc(1,iloc))-&
                     min(rp(i), rmloc(2,iloc)))*df_ij
          else
            df_ij = (min(rm(i), rploc(1,iloc))-&
                     min(rm(i), rploc(2,iloc)))*df_ij
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
    
  end subroutine gfsc_buildJacobianScalarGP

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianBlockSymm(ru, dscale, hstep, bclear, rafcstab,&
                                         rjacobianMatrix, bextendedSparsity)

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
    type(t_vectorBlock), intent(IN) :: ru

    ! scaling parameter
    real(DP), intent(IN) :: dscale

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    if (ru%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockSymm')
      call sys_halt()

    else

      call gfsc_buildJacobianScalarSymm(ru%RvectorBlock(1), dscale, hstep, bclear, &
                                        rafcstab, rjacobianMatrix, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianBlockSymm

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarSymm(ru, dscale, hstep, bclear, rafcstab,&
                                          rjacobianMatrix, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete diffusion operator for a scalar convection equation.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(IN) :: ru

    ! scaling parameter
    real(DP), intent(IN) :: dscale

    ! perturbation parameter
    real(DP), intent(IN) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(IN), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rjacobianMatrix   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_flux,p_u,p_Jac
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IsuperdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagonalEdges
    integer, dimension(:), pointer :: p_IsubdiagonalEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep
    logical :: bisExtended

    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)      .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)   .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_BOUNDS)          .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)          .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarSymm')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobianMatrix)
    
    ! Check if subdiagonal edges need to be generated
    if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
        call afcstab_generateSubdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagonalEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rjacobianMatrix, p_Jac)
    call lsyssc_getbase_double(ru, p_u)
    
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
      
      call doJacobianMat79_Symm(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                p_u, p_flux, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
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
      
      call doJacobianMat79_Symm(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                p_u, p_flux, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
                                dscale, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)

      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarSymm')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    ! is "moved" to the given column "k". For efficiency reasons, only
    ! those entries are considered which are present in column "k".
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: k

      integer, dimension(:), intent(INOUT) :: Ksep
      
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
    subroutine doJacobianMat79_Symm(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                    IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                    DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                    u, flux, pp, pm, qp, qm, rp, rm, dscale, hstep,&
                                    NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,flux,pp,pm,qp,qm,rp,rm
      real(DP), intent(IN) :: dscale,hstep
      integer, dimension(:,:), intent(IN) :: IverticesAtEdge
      integer, dimension(:), intent(IN) :: IsuperdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdgesIdx
      integer, dimension(:), intent(IN) :: IsubdiagonalEdges
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE,NNVEDGE
      logical, intent(IN) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc,qploc,qmloc,rploc,rmloc,fluxloc
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ild,iedge,i,j,k,l,iloc,nloc
      

      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        pploc(:,0) = 0; pmloc(:,0) = 0
        qploc(:,0) = 0; qmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagonalEdgesIdx(k), IsubdiagonalEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagonalEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1

          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)
          
          ! Update local coefficients
          call updateJacobianMat79_Symm(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                        hstep, iedge, i, j, iloc, k,&
                                        pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
             
          ! Determine indices
          i = IverticesAtEdge(1,iedge)
          j = IverticesAtEdge(2,iedge)

          ! Update local coefficients
          call updateJacobianMat79_Symm(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                        hstep, iedge, i, j, iloc, k,&
                                        pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Save total number of local neighbors
        nloc = iloc
        

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit(pploc(:,0:nloc), qploc(:,0:nloc), 0.0_DP, 1.0_DP)
        rmloc(:,0:nloc) = afcstab_limit(pmloc(:,0:nloc), qmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagonalEdgesIdx(l), IsubdiagonalEdgesIdx(l+1)-1
            
            ! Get edge number
            iedge = IsubdiagonalEdges(ild)
            
            call assembleJacobianMat79_Symm(IverticesAtEdge, Kld, Kcol, flux, rp, rm,&
                                            Kloc, rploc, rmloc, fluxloc, dscale, hstep,&
                                            iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
          
          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_Symm(IverticesAtEdge, Kld, Kcol, flux, rp, rm,&
                                            Kloc, rploc, rmloc, fluxloc, dscale, hstep,&
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
    subroutine updateJacobianMat79_Symm(DcoefficientsAtEdge, u, pp, pm, qp, qm,&
                                        hstep, iedge, i, j, iloc, k,&
                                        pploc, pmloc, qploc, qmloc, fluxloc, Kloc)

      real(DP), dimension(:,:), intent(IN) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN) :: u,pp,pm,qp,qm
      real(DP), intent(IN) :: hstep
      integer, intent(IN) :: iedge,i,j,k,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT) :: pploc,pmloc,qploc,qmloc,fluxloc
      integer, dimension(:,:), intent(INOUT) :: Kloc

      ! local variables
      real(DP) :: d_ij,f_ij,s_ij,diff,hstep_ik,hstep_jk,dsign
      integer :: iperturb

      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     "eliminate" the contribution of the edge IJ for the
      !     unperturbed solution values u_i and u_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      s_ij = DcoefficientsAtEdge(2,iedge)

      ! Determine solution difference
      diff = u(i)-u(j)

      if (i .eq. k) then
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0.0_DP
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        pploc(:,iloc) = pp(j)-max(0.0_DP, -f_ij)
        pmloc(:,iloc) = pm(j)-max(0.0_DP, -f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        qploc(:,iloc) = qp(j)-max(0.0_DP, -f_ij)
        qmloc(:,iloc) = qm(j)-min(0.0_DP, -f_ij)
        
      else
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = i
        
        ! Compute signed perturbation parameters
        hstep_ik = 0.0_DP; hstep_jk = hstep
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        pploc(:,iloc) = pp(i)-max(0.0_DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0.0_DP, f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        qploc(:,iloc) = qp(i)-max(0.0_DP, f_ij)
        qmloc(:,iloc) = qm(i)-min(0.0_DP, f_ij)
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
          fluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          pploc(iperturb,0)    = pploc(iperturb,0)+max(0.0_DP, f_ij)
          pmloc(iperturb,0)    = pmloc(iperturb,0)+min(0.0_DP, f_ij)
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, -f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))

          ! Update upper/lower bounds
          qploc(iperturb,0)    = qploc(iperturb,0)+max(0.0_DP, f_ij)
          qmloc(iperturb,0)    = qmloc(iperturb,0)+min(0.0_DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP, -f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP, -f_ij)
          
        else
          
          ! Compute raw antidiffusve flux
          f_ij = d_ij*(diff+dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0.0_DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          pploc(iperturb,0)    = pploc(iperturb,0)+max(0.0_DP, -f_ij)
          pmloc(iperturb,0)    = pmloc(iperturb,0)+min(0.0_DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))
          
          ! Update upper/lower bounds
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0.0_DP, f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          qploc(iperturb,0)    = qploc(iperturb,0)+max(0.0_DP, -f_ij)
          qmloc(iperturb,0)    = qmloc(iperturb,0)+min(0.0_DP, -f_ij)
        end if
      end do
    end subroutine updateJacobianMat79_Symm

    
    !**************************************************************
    ! Assemble the given column of the Jacobian for symmetric flux limiting
    subroutine assembleJacobianMat79_Symm(IverticesAtEdge, Kdiagonal, Kcol, flux,&
                                          rp, rm, Kloc, rploc, rmloc, fluxloc,&
                                          dscale, hstep, iedge, iloc, k, l,&
                                          bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(IN) :: rploc,rmloc,fluxloc
      real(DP), dimension(:), intent(IN) :: rp,rm,flux
      real(DP), intent(IN) :: dscale,hstep
      integer, dimension(:,:), intent(IN)  :: IverticesAtEdge,Kloc
      integer, dimension(:), intent(IN) :: Kdiagonal,Kcol
      integer, intent(IN) :: iedge,iloc,k,l
      logical, intent(IN) :: bisExtended

      real(DP), dimension(:), intent(INOUT) :: Jac
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
          f_ij = fluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit flux 
            if (f_ij > 0.0_DP) then
              f_ij = dscale*min(rploc(iperturb,0), rmloc(iperturb,iloc))*f_ij
            else
              f_ij = dscale*min(rmloc(iperturb,0), rploc(iperturb,iloc))*f_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            if (f_ij > 0.0_DP) then
              f_ij = dscale*min(rploc(iperturb,iloc), rmloc(iperturb,0))*f_ij
            else
              f_ij = dscale*min(rmloc(iperturb,iloc), rploc(iperturb,0))*f_ij
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

          if (flux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*dscale*(min(rploc(1,iloc), rm(j))-&
                                  min(rploc(2,iloc), rm(j)))*flux(iedge)/hstep
          else
            f_ij = 0.5_DP*dscale*(min(rmloc(1,iloc), rp(j))-&
                                  min(rmloc(2,iloc), rp(j)))*flux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        else

          if (flux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*dscale*(min(rp(i), rmloc(1,iloc))-&
                                  min(rp(i), rmloc(2,iloc)))*flux(iedge)/hstep
          else
            f_ij = 0.5_DP*dscale*(min(rm(i), rploc(1,iloc))-&
                                  min(rm(i), rploc(2,iloc)))*flux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        end if

      end if
    end subroutine assembleJacobianMat79_Symm
  end subroutine gfsc_buildJacobianScalarSymm

  !*****************************************************************************

!<function>

  pure function gfsc_hasOrientation(rafcstab) result(bhasOrientation)

!<description>
    ! This function returns .TRUE. if the given stabilisation structure
    ! requires an oriented edge data structure. Otherwise it returns .FALSE.
!</description>

!<input>
    ! stabilisation structure
    type(t_afcstab), intent(IN) :: rafcstab
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

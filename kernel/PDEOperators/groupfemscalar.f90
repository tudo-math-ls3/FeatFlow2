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
!#     -> assemble the residual vector for AFC stabilisation of FCT type
!#
!# 7.) gfsc_buildResidualTVD = gfsc_buildResScalarTVD /
!#                             gfsc_buildResBlockTVD
!#     -> assemble the residual vector for AFC stabilisation of TVD type
!#
!# 8.) gfsc_buildResidualSymm = gfsc_buildResScalarSymm /
!#                                    gfsc_buildResBlockSymm
!#     -> assemble the residual vector for stabilisation by means of
!#        symmetric flux limiting for diffusion operators
!#
!# 9.) gfsc_buildConvectionJacobian = gfsc_buildConvJacobianScalar /
!#                                    gfsc_buildConvJacobianBlock
!#     -> assemble the Jacobian matrix for the convective part of
!#        the transport operator for a scalar convection equation
!#
!# 10.) gfsc_buildJacobianFCT = gfsc_buildJacLinearScalarFCT /
!#                              gfsc_buildJacLinearBlockFCT /
!#                              gfsc_buildJacobianScalarFCT /
!#                              gfsc_buildJacobianBlockFCT
!#      -> assemble the Jacobian matrix for the stabilisation part of FCT type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 11.) gfsc_buildJacobianTVD = gfsc_buildJacLinearScalarTVD /
!#                              gfsc_buildJacLinearBlockTVD /
!#                              gfsc_buildJacobianScalarTVD /
!#                              gfsc_buildJacobianBlockTVD
!#      -> assemble the Jacobian matrix for the stabilisation part of TVD type
!#         and/or for the general purpose limiter; 
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 12.) gfsc_buildJacobianSymm = gfsc_buildJacobianScalarSymm /
!#                               gfsc_buildJacobianBlockSymm
!#      -> assemble the Jacobian matrix for the stabilisation part of symmetric
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
  public :: gfsc_buildResidualSymm
  public :: gfsc_buildConvectionJacobian
  public :: gfsc_buildJacobianFCT
  public :: gfsc_buildJacobianTVD
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
          AFCSTAB_FEMFCT_EXP,&
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

  subroutine gfsc_buildConvOperatorBlock(RcoeffMatrices, rsolution, fcb_calcConvection,&
                                         bZeroRowsum, bStabilise, bclear, rconvMatrix, rafcstab)
    
!<description>
    ! This subroutine assembles the discrete transport operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !     $$ \nabla\cdot({\bf v}u)$$
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
    type(t_vectorBlock), intent(IN) :: rsolution
    
    ! Switch for row-sum property
    ! TRUE  : assume zero row-sum, that is, the diegonal coefficient
    !         is computed as the sum of negative off-diagonal coefficients
    ! FALSE : compute the diagonal coefficient by standard finite elements
    logical, intent(IN) :: bZeroRowsum

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
    ! The transport operator
    type(t_matrixScalar), intent(INOUT) :: rconvMatrix
    
    ! OPTIONAL: the stabilisation structure
    type(t_afcstab), intent(INOUT), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (rsolution%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvOperatorBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvOperatorScalar(RcoeffMatrices, rsolution%RvectorBlock(1),&
                                        fcb_calcConvection, bZeroRowsum, bStabilise,&
                                        bclear, rconvMatrix, rafcstab)

    end if
  end subroutine gfsc_buildConvOperatorBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvOperatorScalar(RcoeffMatrices, rsolution, fcb_calcConvection,&
                                          bZeroRowsum, bStabilise, bclear, rconvMatrix, rafcstab)

!<description>
    ! This subroutine assembles the discrete transport operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !     $$ \nabla\cdot({\bf v}u)$$
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
    type(t_vectorScalar), intent(IN) :: rsolution

    ! Switch for row-sum property
    ! TRUE  : assume zero row-sum, that is, the diegonal coefficient
    !         is computed as the sum of negative off-diagonal coefficients
    ! FALSE : compute the diagonal coefficient by standard finite elements
    logical, intent(IN) :: bZeroRowsum

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
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal,p_IsuperdiagonalEdgesIdx
    integer :: h_Ksep
    integer :: idim,ndim
    
    ! Check if all matrices/vectors are compatible
    call lsyssc_isMatrixCompatible(rsolution, rconvMatrix, .false.)

    ndim = size(RcoeffMatrices,1)
    do idim = 1, ndim
      call lsyssc_isMatrixCompatible(RcoeffMatrices(idim), rconvMatrix)
    end do

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rconvMatrix)

    ! Set pointer
    call lsyssc_getbase_Kld(rconvMatrix, p_Kld)
    call lsyssc_getbase_Kcol(rconvMatrix, p_Kcol)
    call lsyssc_getbase_double(rconvMatrix, p_ConvOp)
    call lsyssc_getbase_double(rsolution, p_u)

    ! How many dimensions do we have?
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
        
        ! Check if matrix/vector and stabilisation
        ! structure are compatible to each other
        call gfsc_isMatrixCompatible(rafcstab, rconvMatrix)
        call gfsc_isVectorCompatible(rafcstab, rsolution)

        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        ! Do we need edge orientation?
        if (gfsc_hasOrientation(rafcstab)) then
          
          ! Adopt orientation convention IJ, such that L_ij < L_ji
          ! and generate edge structure for the flux limiter
          
          if (bZeroRowsum) then
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindZRSMat7_ordAFC_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                             p_Cx, p_u, p_ConvOp,&
                                             p_IsuperdiagonalEdgesIdx,&
                                             p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindZRSMat7_ordAFC_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                             p_Cx, p_Cy, p_u, p_ConvOp,&
                                             p_IsuperdiagonalEdgesIdx,&
                                             p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindZRSMat7_ordAFC_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                             p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                             p_IsuperdiagonalEdgesIdx,&
                                             p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else   ! bZeroRowsum
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindMat7_ordAFC_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindMat7_ordAFC_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindMat7_ordAFC_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if   ! bZeroRowsum

          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION)

        else   ! bhasOrientation

          ! Adopt no orientation convention and generate edge structure

          if (bZeroRowsum) then
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindZRSMat7_AFC_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindZRSMat7_AFC_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindZRSMat7_AFC_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                          p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else   ! bZeroRowsum

            select case(ndim)
            case (NDIM1D)
              call doUpwindMat7_AFC_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                       p_Cx, p_u, p_ConvOp,&
                                       p_IsuperdiagonalEdgesIdx,&
                                       p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindMat7_AFC_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                       p_Cx, p_Cy, p_u, p_ConvOp,&
                                       p_IsuperdiagonalEdgesIdx,&
                                       p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindMat7_AFC_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                       p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                       p_IsuperdiagonalEdgesIdx,&
                                       p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if   ! bZeroRowsum
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

        end if   ! bhasOrientation
        
      elseif (bStabilise) then   ! present(rafcstab)
        
        ! Perform discrete upwinding but do not generate the edge data structure

        if (bZeroRowsum) then
          
          select case(ndim)
          case (NDIM1D)
            call doUpwindZRSMat7_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindZRSMat7_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindZRSMat7_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                    p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select
          
        else   ! bZeroRowsum

          select case(ndim)
          case (NDIM1D)
            call doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                 p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                 p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                 p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        end if   !bZeroRowsum

      else   ! present(rafcstab) and bStabilise
        
        ! Apply standard Galerkin discretisation

        if (bZeroRowsum) then
          
          select case(ndim)
          case (NDIM1D)
            call doGalerkinZRSMat7_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinZRSMat7_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinZRSMat7_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                      p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        else   ! bZeroRowsum

          select case(ndim)
          case (NDIM1D)
            call doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                   p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                   p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep, rconvMatrix%NEQ,&
                                   p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select

        end if   ! bZeroRowsum

      end if   ! present(rafcstab)


      ! Release diagonal separator
      call storage_free(h_Ksep)


    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointer to diagonal
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
        
        ! Check if matrix/vector and stabilisation
        ! structure are compatible to each other
        call gfsc_isMatrixCompatible(rafcstab, rconvMatrix)
        call gfsc_isVectorCompatible(rafcstab, rsolution)

        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)

        ! Do we need edge orientation?
        if (gfsc_hasOrientation(rafcstab)) then

          ! Adopt orientation convention IJ, such that L_ij < L_ji 
          ! and generate edge structure for the flux limiter
          if (bZeroRowsum) then
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindZRSMat9_ordAFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                             rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                             p_IsuperdiagonalEdgesIdx,&
                                             p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindZRSMat9_ordAFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                             rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                             p_IsuperdiagonalEdgesIdx,&
                                             p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindZRSMat9_ordAFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                             rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                             p_IsuperdiagonalEdgesIdx,&
                                             p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select
            
          else   ! bZeroRowsum

            select case(ndim)
            case (NDIM1D)
              call doUpwindMat9_ordAFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindMat9_ordAFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindMat9_ordAFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if   ! bZeroRowsum
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEORIENTATION)

        else   ! bhasOrientation

          ! Adopt no orientation convention and generate edge structure
          
          if (bZeroRowsum) then
            
            select case(ndim)
            case (NDIM1D)
              call doUpwindZRSMat9_AFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindZRSMat9_AFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindZRSMat9_AFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                          rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                          p_IsuperdiagonalEdgesIdx,&
                                          p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          else   ! bZeroRowsum

            select case(ndim)
            case (NDIM1D)
              call doUpwindMat9_AFC_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                       rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp,&
                                       p_IsuperdiagonalEdgesIdx,&
                                       p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM2D)
              call doUpwindMat9_AFC_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                       rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp,&
                                       p_IsuperdiagonalEdgesIdx,&
                                       p_IverticesAtEdge, p_DcoefficientsAtEdge)
            case (NDIM3D)
              call doUpwindMat9_AFC_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                       rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp,&
                                       p_IsuperdiagonalEdgesIdx,&
                                       p_IverticesAtEdge, p_DcoefficientsAtEdge)
            end select

          end if   ! bZeroRowsum
          
          ! Set state of stabilisation
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE)
          rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_EDGEVALUES)

        end if   ! bhasOrientation
          
      elseif (bStabilise) then   ! present(rafcstab)
        
        ! Perform discrete upwinding but do not generate the edge data structure
        
        if (bZeroRowsum) then
          
          select case(ndim)
          case (NDIM1D)
            call doUpwindZRSMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindZRSMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindZRSMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select
          
        else   ! bZeroRowsum
          
          select case(ndim)
          case (NDIM1D)
            call doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select
          
        end if   ! bZeroRowsum
        
      else   ! present(rafcstab) and bStabilise
        
        ! Apply standard Galerkin discretisation
        
        if (bZeroRowsum) then
          
          select case(ndim)
          case (NDIM1D)
            call doGalerkinZRSMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinZRSMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinZRSMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                      rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select
          
        else   ! bZeroRowsum
          
          select case(ndim)
          case (NDIM1D)
            call doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                   rconvMatrix%NEQ, p_Cx, p_u, p_ConvOp)
          case (NDIM2D)
            call doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                   rconvMatrix%NEQ, p_Cx, p_Cy, p_u, p_ConvOp)
          case (NDIM3D)
            call doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                   rconvMatrix%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_ConvOp)
          end select
          
        end if   ! bZeroRowsum
        
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
    ! Assemble high-order Galerkin operator K in 1D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7
    
    subroutine doGalerkinZRSMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
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
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        end do
      end do
    end subroutine doGalerkinZRSMat7_1D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 7
    
    subroutine doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, K)

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
        K(ii) = K(ii)+k_ii

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
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        end do
      end do
    end subroutine doGalerkinMat7_1D

    
    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7
    
    subroutine doGalerkinZRSMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
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
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        end do
      end do
    end subroutine doGalerkinZRSMat7_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 7
    
    subroutine doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, K)

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
        K(ii) = K(ii)+k_ii
        
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
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        end do
      end do
    end subroutine doGalerkinMat7_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! and assume zero row-sums
    ! All matrices are stored in matrix format 7
    
    subroutine doGalerkinZRSMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
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
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        end do
      end do
    end subroutine doGalerkinZRSMat7_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 7
    
    subroutine doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, K)

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
        K(ii) = K(ii)+k_ii

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
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        end do
      end do
    end subroutine doGalerkinMat7_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D
    ! and assume zero row-sums
    ! All matrices are stored in matrix format 9
    
    subroutine doGalerkinZRSMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
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
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        end do
      end do
    end subroutine doGalerkinZRSMat9_1D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 1D.
    ! All matrices are stored in matrix format 9
    
    subroutine doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, K)

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
        K(ii) = K(ii)+k_ii
        
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
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        end do
      end do
    end subroutine doGalerkinMat9_1D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D
    ! and assume zero row-sums
    ! All matrices are stored in matrix format 9
    
    subroutine doGalerkinZRSMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
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
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        end do
      end do
    end subroutine doGalerkinZRSMat9_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 2D.
    ! All matrices are stored in matrix format 9
    
    subroutine doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, K)

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
        K(ii) = K(ii)+k_ii

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
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        end do
      end do
    end subroutine doGalerkinMat9_2D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9
    
    subroutine doGalerkinZRSMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, K)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,Cz,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ

      real(DP), dimension(:), intent(INOUT) :: K
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji
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
          
          ! Assemble the global operator
          K(ij) = K(ij)+k_ij; K(ii) = K(ii)-k_ij
          K(ji) = K(ji)+k_ji; K(jj) = K(jj)-k_ji
        end do
      end do
    end subroutine doGalerkinZRSMat9_3D


    !**************************************************************
    ! Assemble high-order Galerkin operator K in 3D.
    ! All matrices are stored in matrix format 9
    
    subroutine doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, K)

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
        K(ii) = K(ii)+k_ii

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
          K(ij) = K(ij)+k_ij; K(ji) = K(ji)+k_ji
        end do
      end do
    end subroutine doGalerkinMat9_3D


    !**************************************************************
    ! Assemble low-order operator L in 1D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindZRSMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L)

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        end do
      end do
    end subroutine doUpwindZRSMat7_1D


    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindMat7_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L)

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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doUpwindMat7_1D


    !**************************************************************
    ! Assemble low-order operator L in 2D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindZRSMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L)

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        end do
      end do
    end subroutine doUpwindZRSMat7_2D


    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindMat7_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L)

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
        L(ii) = L(ii)+k_ii

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doUpwindMat7_2D


    !**************************************************************
    ! Assemble low-order operator L in 3D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindZRSMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L)

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        end do
      end do
    end subroutine doUpwindZRSMat7_3D


    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindMat7_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L)

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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doUpwindMat7_3D

    
    !**************************************************************
    ! Assemble low-order operator L in 1D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 9
    
    subroutine doUpwindZRSMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L)

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        end do
      end do
    end subroutine doUpwindZRSMat9_1D


    !**************************************************************
    ! Assemble low-order operator L in 1D.
    ! All matrices are stored in matrix format 9
    
    subroutine doUpwindMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L)

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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doUpwindMat9_1D


    !**************************************************************
    ! Assemble low-order operator L in 2D
    ! and assume zero row-sum.
    ! All matrices are stored in matrix format 9
    
    subroutine doUpwindZRSMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L)

      real(DP), dimension(:), intent(IN) :: Cx,Cy,u
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        end do
      end do
    end subroutine doUpwindZRSMat9_2D


    !**************************************************************
    ! Assemble low-order operator L in 2D.
    ! All matrices are stored in matrix format 9
    
    subroutine doUpwindMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L)

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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doUpwindMat9_2D


    !**************************************************************
    ! Assemble low-order operator L in 3D
    ! and assume zero row-sums.
    ! All matrices are stored in matrix format 9
    
    subroutine doUpwindZRSMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L)

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
        end do
      end do
    end subroutine doUpwindZRSMat9_3D


    !**************************************************************
    ! Assemble low-order operator L in 3D.
    ! All matrices are stored in matrix format 9
    
    subroutine doUpwindMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L)

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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
        end do
      end do
    end subroutine doUpwindMat9_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindZRSMat7_AFC_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L, &
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindZRSMat7_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7
    
    subroutine doUpwindMat7_AFC_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,& 
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindMat7_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindZRSMat7_AFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindZRSMat7_AFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_AFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindMat7_AFC_2D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindZRSMat7_AFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindZRSMat7_AFC_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_AFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindMat7_AFC_3D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindZRSMat9_AFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
                                      Cx, u, L, IsuperdiagonalEdgesIdx,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindZRSMat9_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_AFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindMat9_AFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindZRSMat9_AFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindZRSMat9_AFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_AFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindMat9_AFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge 
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindZRSMat9_AFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindZRSMat9_AFC_3D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data w/o edge 
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_AFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
          ! Increase edge counter
          iedge = iedge+1
          
          ! AFC w/o edge orientation
          IverticesAtEdge(:,iedge)     = (/i, j, ij, ji/)
          DcoefficientsAtEdge(:,iedge) = (/d_ij, k_ij, k_ji/)          
        end do
      end do
      
      ! Set index for last entry
      IsuperdiagonalEdgesIdx(NEQ+1) = iedge+1
    end subroutine doUpwindMat9_AFC_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindZRSMat7_ordAFC_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
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
    end subroutine doUpwindZRSMat7_ordAFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_ordAFC_1D(Kld, Kcol, Ksep, NEQ, Cx, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
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
    end subroutine doUpwindMat7_ordAFC_1D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindZRSMat7_ordAFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
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
    end subroutine doUpwindZRSMat7_ordAFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_ordAFC_2D(Kld, Kcol, Ksep, NEQ, Cx, Cy, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
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
    end subroutine doUpwindMat7_ordAFC_2D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindZRSMat7_ordAFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
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
    end subroutine doUpwindZRSMat7_ordAFC_3D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_ordAFC_3D(Kld, Kcol, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
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
    end subroutine doUpwindMat7_ordAFC_3D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindZRSMat9_ordAFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L, &
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
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
    end subroutine doUpwindZRSMat9_ordAFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 1D.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_ordAFC_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
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
    end subroutine doUpwindMat9_ordAFC_1D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindZRSMat9_ordAFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
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
    end subroutine doUpwindZRSMat9_ordAFC_2D


    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 2D.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_ordAFC_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, u, L,&
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
        L(ii) = L(ii)+k_ii

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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
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
    end subroutine doUpwindMat9_ordAFC_2D

    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindZRSMat9_ordAFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-k_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-k_ji
          
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
    end subroutine doUpwindZRSMat9_ordAFC_3D
    
    
    !**************************************************************
    ! Assemble low-order operator L and AFC data with edge
    ! orientation in 3D.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_ordAFC_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ, Cx, Cy, Cz, u, L,&
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
        L(ii) = L(ii)+k_ii
        
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
          d_ij = max(-k_ij, 0._DP, -k_ji)
          
          ! Apply artificial diffusion
          k_ij = k_ij+d_ij
          k_ji = k_ji+d_ij
          
          ! Assemble the global operator
          L(ij) = L(ij)+k_ij; L(ii) = L(ii)-d_ij
          L(ji) = L(ji)+k_ji; L(jj) = L(jj)-d_ij
          
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
    end subroutine doUpwindMat9_ordAFC_3D

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
    real(DP), dimension(:), pointer :: p_S,p_L
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal,p_IsuperdiagonalEdgesIdx
    integer :: h_Ksep

    
    ! Check if matrices are compatible
    call lsyssc_isMatrixCompatible(rcoeffMatrix, rdiffMatrix)
    
    ! Should matrix be cleared?
    if (bclear) call lsyssc_clearMatrix(rdiffMatrix)

    ! Set pointers
    call lsyssc_getbase_double(rcoeffMatrix, p_S)
    call lsyssc_getbase_double(rdiffMatrix, p_L)      

    
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

        ! Check if matrix and stabilisation
        ! structure are compatible to each other
        call gfsc_isMatrixCompatible(rafcstab, rdiffMatrix)

        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    
        call doLoworderMat7_AFC(p_Kld, p_Kcol, p_Ksep, rdiffMatrix%NEQ,&
                                p_S, p_L, p_IsuperdiagonalEdgesIdx,&
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

        call doLoworderMat7(p_Kld, p_Kcol, p_Ksep, rdiffMatrix%NEQ, p_S, p_L)

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

        ! Check if matrix and stabilisation
        ! structure are compatible to each other
        call gfsc_isMatrixCompatible(rafcstab, rdiffMatrix)

        ! Set additional pointers
        call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
        call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        call doLoworderMat9_AFC(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rdiffMatrix%NEQ,&
                                p_S, p_L, p_IsuperdiagonalEdgesIdx,&
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

        call doLoworderMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, rdiffMatrix%NEQ, p_S, p_L)

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
          d_ij = max(0._DP, -S(ij)) 
          
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
          d_ij = max(0._DP, -S(ij)) 
          
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
          d_ij = max(0._DP, -S(ij))
          s_ij = max(0._DP,  S(ij))
          
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
          d_ij = max(0._DP, -S(ij))
          s_ij = max(0._DP,  S(ij))
          
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

  subroutine gfsc_buildResBlockFCT(rmatrixMC, rmatrixML, ru,&
                                   theta, tstep, binit, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-FCT type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)   :: rmatrixMC

    ! lumped mass matrix
    type(t_matrixScalar), intent(IN)   :: rmatrixML

    ! solution vector
    type(t_vectorBlock), intent(IN)    :: ru

    ! implicitness parameter
    real(DP), intent(IN)               :: theta

    ! time step size
    real(DP), intent(IN)               :: tstep
    
    ! Switch for residual
    ! TRUE  : build the initial residual
    ! FALSE : build an intermediate residual
    logical, intent(IN)                :: binit
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT)     :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (ru%nblocks .ne. 1 .or. rres%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockFCT')
      call sys_halt()

    else

      call gfsc_buildResScalarFCT(rmatrixMC, rmatrixML, ru%RvectorBlock(1),&
          theta, tstep, binit, rres%RvectorBlock(1), rafcstab)

    end if
  end subroutine gfsc_buildResBlockFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResScalarFCT(rmatrixMC, rmatrixML, ru,&
                                    theta, tstep, binit, rres, rafcstab)

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
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)    :: rmatrixMC

    ! lumped mass matrix
    type(t_matrixScalar), intent(IN)    :: rmatrixML

    ! solution vector
    type(t_vectorScalar), intent(IN)    :: ru

    ! implicitness parameter
    real(DP), intent(IN)                :: theta

    ! time step size
    real(DP), intent(IN)                :: tstep
    
    ! Switch for residual
    ! TRUE  : build the initial residual
    ! FALSE : build an intermediate residual
    logical, intent(IN)                 :: binit
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorScalar), pointer                 :: ruLow
    integer(PREC_VECIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_flux,p_flux0
    real(DP), dimension(:), pointer :: p_u,p_ulow,p_res
    real(DP), dimension(:), pointer :: p_MC,p_ML

    ! Check if vectors are compatible
    call lsyssc_isVectorCompatible(ru, rres)
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
      call sys_halt()
    end if
 
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
    call lsyssc_getbase_double(rmatrixMC, p_MC)
    call lsyssc_getbase_double(rmatrixML, p_ML)
    call lsyssc_getbase_double(ru,   p_u)
    call lsyssc_getbase_double(rres, p_res)

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT)

      ! Should we build up the initial residual?
      if (binit) then
        
        ! Do we have a fully implicit time discretisation?
        if (theta < 1.0_DP) then

          ! Compute the low-order predictor
          !
          ! $$\tilde u=u^n+(1-\theta)\Delta tM_L^{-1}Lu^n$$
          ! 
          ! whereby the residual of the 0-th iteration is assumed to
          ! be $r^{(0)}=\Delta tLu^n$.
          ruLow => rafcstab%RnodalVectors(5)
          call lsyssc_invertedDiagMatVec(rmatrixML, rres, 1._DP, ruLow)
          call lsyssc_vectorLinearComb(ru, ruLow, 1._DP, 1._DP-theta)
          call lsyssc_getbase_double(ruLow, p_ulow)
          
          ! Initialise the flux limiter
          call doInit_implFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_MC, p_ML, p_u, p_ulow, theta, tstep, rafcstab%NEDGE,&
              (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS),&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0)
          
        else
          
          ! The low-order predictor is simply given by
          !
          ! $$\tilde u=u^n$$
          !
          ! Initialise the flux limiter with u in leu of ulow
          call doInit_implFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_MC, p_ML, p_u, p_u, theta, tstep, rafcstab%NEDGE,&
              (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS),&
              p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0)
          
        end if
        
        ! Set specifier
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)        
      end if
      
      ! Check if correction factors and fluxes are available
      if (iand(rafcstab%iSpec, AFCSTAB_LIMITER) .eq. 0 .or.&
          iand(rafcstab%iSpec, AFCSTAB_FLUXES)  .eq. 0) then
        call output_line('Stabilisation does not provide precomputed fluxes &
            &and/or nodal correction factors',OU_CLASS_ERROR,OU_MODE_STD,&
            'gfsc_buildResScalarFCT')
        call sys_halt()
      end if
      
      ! Apply the limited antidiffusion
      call doLimit_implFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
          p_MC, p_u, p_flux, p_flux0, theta, tstep, rafcstab%NEDGE,&
          (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_res)


    case (AFCSTAB_FEMFCT_EXP)

      ! Should we build up the initial residual?
      if (binit) then
        ! Initialise the flux limiter
        call doInit_explFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_MC, p_u, theta, tstep, rafcstab%NEDGE,&
            (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_flux0)

        ! Set specifier
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
        rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES) 
      end if

      ! Check if correction factors and fluxes are available
      if (iand(rafcstab%iSpec, AFCSTAB_LIMITER) .eq. 0 .or.&
          iand(rafcstab%iSpec, AFCSTAB_FLUXES)  .eq. 0) then
        call output_line('Stabilisation does not provide precomputed fluxes &
            &and/or nodal correction factors',OU_CLASS_ERROR,OU_MODE_STD,&
            'gfsc_buildResScalarFCT')
        call sys_halt()
      end if

      ! Do we have a fully implicit time discretisation?
      if (theta < 1.0_DP) then

        ! Compute the low-order predictor
        !
        ! $$\tilde u=u^n+(1-\theta)\Delta tM_L^{-1}Lu^n$$
        ! 
        ! whereby the residual of the 0-th iteration is assumed to
        ! be $r^{(0)}=\Delta tLu^n$.
        ruLow => rafcstab%RnodalVectors(5)
        call lsyssc_invertedDiagMatVec(rmatrixML, rres, 1._DP, ruLow)
        call lsyssc_vectorLinearComb(ru, ruLow, 1._DP, 1._DP-theta)
        call lsyssc_getbase_double(ruLow, p_ulow)
        
        ! Apply the flux limiter
        call doLimit_explFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_MC, p_ML, p_u, p_ulow, p_flux0, theta, tstep, rafcstab%NEDGE,&
            (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS),&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)
        
      else
        
        ! The low-order predictor is simply given by
        !
        ! $$\tilde u=u^n$$
        !
        ! Apply the flux limiter with u in leu of ulow
        call doLimit_explFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_MC, p_ML, p_u, p_u, p_flux0, theta, tstep, rafcstab%NEDGE,&
            (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS),&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

      end if


    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarFCT')
      call sys_halt()
    end select

  contains
      
    ! Here, the working routine follow

    !**************************************************************
    ! Initialisation of the semi-implicit FEM-FCT procedure
    
    subroutine doInit_implFCT(IverticesAtEdge, DcoefficientsAtEdge,&
                              MC, ML, u, ulow, theta, tstep, NEDGE, bmass,&
                              pp, pm, qp, qm, rp, rm, flux, flux0)
      
      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: MC,ML
      real(DP), dimension(:), intent(IN)               :: u,ulow
      real(DP), intent(IN)                             :: theta,tstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: pp,pm
      real(DP), dimension(:), intent(INOUT)            :: qp,qm
      real(DP), dimension(:), intent(INOUT)            :: rp,rm
      real(DP), dimension(:), intent(INOUT)            :: flux,flux0
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,m_ij,diff

      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)

      ! Should we apply the consistent mass matrix?
      if (bmass) then

        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
          
          ! Determine fluxes
          diff = u(i)-u(j);   f_ij = tstep*d_ij*diff
          flux(iedge) = f_ij; flux0(iedge) = -m_ij*diff
          
          ! Sum of positive/negative fluxes
          pp(i) = pp(i)+max(0._DP,f_ij); pp(j) = pp(j)+max(0._DP,-f_ij)
          pm(i) = pm(i)+min(0._DP,f_ij); pm(j) = pm(j)+min(0._DP,-f_ij)
          
          ! Upper/lower bounds
          diff = ulow(j)-ulow(i)
          qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
          qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
        end do

        ! Adopt the explicit part (if required)
        if (theta < 1.0_DP) then
          call lalg_vectorLinearCombDble(flux, flux0, 1.0_DP-theta, 1.0_DP)
        end if
        
      else

        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          
          ! Determine fluxes
          diff = u(i)-u(j);   f_ij = tstep*d_ij*diff
          flux(iedge) = f_ij
          
          ! Sum of positive/negative fluxes
          pp(i) = pp(i)+max(0._DP,f_ij); pp(j) = pp(j)+max(0._DP,-f_ij)
          pm(i) = pm(i)+min(0._DP,f_ij); pm(j) = pm(j)+min(0._DP,-f_ij)
          
          ! Upper/lower bounds
          diff = ulow(j)-ulow(i)
          qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
          qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
        end do

        ! Adopt the explicit part (if required)
        if (theta < 1.0_DP) then
          call lalg_copyVectorDble(flux, flux0)
          call lalg_scaleVectorDble(flux0, 1.0_DP-theta)
        else
          call lalg_clearVectorDble(flux0)
        end if

      end if

      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit( pp, rp, 0._DP)
      rm =-ML*qm; rm = afcstab_limit(-pm, rm, 0._DP)
      
      ! Limiting procedure
      do iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        if (flux(iedge) > 0.0_DP) then
          flux(iedge) = min(rp(i), rm(j))*flux(iedge)
        else
          flux(iedge) = min(rm(i), rp(j))*flux(iedge)
        end if
      end do
    end subroutine doInit_implFCT
    

    !**************************************************************
    ! The semi-implicit FEM-FCT limiting procedure
    
    subroutine doLimit_implFCT(IverticesAtEdge, DcoefficientsAtEdge,&
                               MC, u, flux, flux0, theta, tstep, NEDGE, bmass, res)

      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: MC
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), dimension(:), intent(IN)               :: flux,flux0
      real(DP), intent(IN)                             :: theta,tstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: res
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,m_ij

      ! Should we apply the consistent mass matrix
      if (bmass) then
        
        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij)
          
          ! Determine fluxes
          f_ij = (m_ij+theta*tstep*d_ij)*(u(i)-u(j))+flux0(iedge)
          
          if (f_ij > 0.0_DP) then
            f_ij = min(f_ij,max(flux(iedge),0._DP))
          else
            f_ij = max(f_ij,min(flux(iedge),0._DP))
          end if
          
          ! Update the defect vector
          res(i) = res(i)+f_ij
          res(j) = res(j)-f_ij
        end do
        
      else
        
        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          
          ! Determine fluxes
          f_ij = (theta*tstep*d_ij)*(u(i)-u(j))+flux0(iedge)
          
          if (f_ij > 0.0_DP) then
            f_ij = min(f_ij,max(flux(iedge),0._DP))
          else
            f_ij = max(f_ij,min(flux(iedge),0._DP))
          end if
          
          ! Update the defect vector
          res(i) = res(i)+f_ij
          res(j) = res(j)-f_ij
        end do
        
      end if
    end subroutine doLimit_implFCT


    !**************************************************************
    ! Initialisation of the semi-explicit FEM-FCT procedure
    
    subroutine doInit_explFCT(IverticesAtEdge, DcoefficientsAtEdge,&
                              MC, u, theta, tstep, NEDGE, bmass, flux0)
      
      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: MC
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: theta,tstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: flux0
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,m_ij,diff

      ! Should we apply the consistent mass matrix?
      if (bmass) then

        ! Should we use semi-implicit scheme?
        if (theta < 1._DP) then
          
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
            flux0(iedge) = -m_ij*diff+(1-theta)*tstep*d_ij*diff
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
            flux0(iedge) = -m_ij*diff
          end do
          
        end if

      else

        ! Should we use semi-implicit scheme?
        if (theta < 1._DP) then
          
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
            flux0(iedge) = (1-theta)*tstep*d_ij*diff
          end do

        else
          
          ! Initialise explicit fluxes by zero
          call lalg_clearVectorDble(flux0)
          
        end if
        
      end if
    end subroutine doInit_explFCT


    !**************************************************************
    ! The semi-explicit FEM-FCT limiting procedure
    
    subroutine doLimit_explFCT(IverticesAtEdge, DcoefficientsAtEdge,&
                              MC, ML, u, ulow, flux0, theta, tstep, NEDGE, bmass,&
                              pp, pm, qp, qm, rp, rm, flux, res)

      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: MC,ML
      real(DP), dimension(:), intent(IN)               :: u,ulow
      real(DP), dimension(:), intent(IN)               :: flux0
      real(DP), intent(IN)                             :: theta,tstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: pp,pm
      real(DP), dimension(:), intent(INOUT)            :: qp,qm
      real(DP), dimension(:), intent(INOUT)            :: rp,rm
      real(DP), dimension(:), intent(INOUT)            :: flux
      real(DP), dimension(:), intent(INOUT)            :: res
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: diff,d_ij,f_ij,m_ij
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(pp)
      call lalg_clearVectorDble(pm)
      call lalg_clearVectorDble(qp)
      call lalg_clearVectorDble(qm)


      ! Should we apply the consistent mass matrix
      if (bmass) then

        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          ij = IverticesAtEdge(3,iedge)
          
          ! Determine coefficients and solution difference
          d_ij = DcoefficientsAtEdge(1,iedge); m_ij = MC(ij); diff=u(i)-u(j)
          
          ! Determine antidiffusive flux
          f_ij = flux0(iedge)+m_ij*diff+theta*tstep*d_ij*diff

          ! Determine low-order solution difference
          diff = ulow(j)-ulow(i)

          ! Perform prelimiting
          if (f_ij*diff .ge. 0) f_ij = 0._DP         
          flux(iedge) = f_ij
          
          ! Sum of positive/negative fluxes
          pp(i) = pp(i)+max(0._DP,f_ij); pp(j) = pp(j)+max(0._DP,-f_ij)
          pm(i) = pm(i)+min(0._DP,f_ij); pm(j) = pm(j)+min(0._DP,-f_ij)
          
          ! Upper/lower bounds
          qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
          qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
        end do
         
      else

        ! Loop over edges
        do iedge = 1, NEDGE
          
          ! Determine indices
          i  = IverticesAtEdge(1,iedge)
          j  = IverticesAtEdge(2,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge); diff=u(i)-u(j)
          
          ! Determine antidiffusive flux
          f_ij = flux0(iedge)+theta*tstep*d_ij*diff

          ! Determine low-order solution difference
          diff = ulow(j)-ulow(i)

          ! Perform prelimiting
          if (f_ij*diff .ge. 0) f_ij = 0._DP         
          flux(iedge) = f_ij
          
          ! Sum of positive/negative fluxes
          pp(i) = pp(i)+max(0._DP,f_ij); pp(j) = pp(j)+max(0._DP,-f_ij)
          pm(i) = pm(i)+min(0._DP,f_ij); pm(j) = pm(j)+min(0._DP,-f_ij)
          
          ! Upper/lower bounds
          qp(i) = max(qp(i),diff); qp(j) = max(qp(j),-diff)
          qm(i) = min(qm(i),diff); qm(j) = min(qm(j),-diff)
        end do
        
      end if

      ! Apply the nodal limiter
      rp = ML*qp; rp = afcstab_limit( pp, rp, 0._DP, 1._DP)
      rm =-ML*qm; rm = afcstab_limit(-pm, rm, 0._DP, 1._DP)

      ! Limiting procedure
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        if (flux(iedge) > 0.0_DP) then
          f_ij = min(rp(i), rm(j))*flux(iedge)
        else
          f_ij = min(rm(i), rp(j))*flux(iedge)
        end if

        ! Update the defect vector
        res(i) = res(i)+f_ij
        res(j) = res(j)-f_ij
      end do
    end subroutine doLimit_explFCT
  end subroutine gfsc_buildResScalarFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResBlockTVD(rmatrixMC, ru, ru0,&
                                   theta, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector
    ! and applies stabilisation of FEM-TVD type.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)   :: rmatrixMC

    ! solution vector
    type(t_vectorBlock), intent(IN)    :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN)    :: ru0

    ! implicitness parameter
    real(DP), intent(IN)               :: theta

    ! time step size
    real(DP), intent(IN)               :: tstep
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT)     :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (ru%nblocks   .ne. 1 .or.&
        ru0%nblocks  .ne. 1 .or.&
        rres%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResBlockTVD')
      call sys_halt()

    else

      call gfsc_buildResScalarTVD(rmatrixMC, ru%RvectorBlock(1),&
          ru0%RvectorBlock(1), theta, tstep, rres%RvectorBlock(1), rafcstab)
      
    end if
  end subroutine gfsc_buildResBlockTVD
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResScalarTVD(rmatrixMC, ru, ru0,&
                                    theta, tstep, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies
    ! stabilisation of FEM-TVD type and/or employes the general
    ! purpose limiter.
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
    !
    ! For time-dependent problems, the lack of a consistent mass
    ! matrix deteriorates the phase accuracy significnatly. It is thus
    ! that symmetric limiting strategies taylored to mass
    ! antidiffusion have been considered. If both convective and mass
    ! antidiffusion has to be limited then the so-called
    ! general-purpose (GP) flux limiter comes into play:
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
    type(t_matrixScalar), intent(IN)    :: rmatrixMC

    ! solution vector
    type(t_vectorScalar), intent(IN)    :: ru

    ! initial solution vector
    type(t_vectorScalar), intent(IN)    :: ru0

    ! implicitness parameter
    real(DP), intent(IN)                :: theta

    ! time step size
    real(DP), intent(IN)                :: tstep
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_VECIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_flux,p_flux0
    real(DP), dimension(:), pointer :: p_u,p_u0,p_res
    real(DP), dimension(:), pointer :: p_MC

    ! Check if vectors are compatible
    call lsyssc_isVectorCompatible(ru, rres)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMTVD)
      
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
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_double(ru,   p_u)
      call lsyssc_getbase_double(rres, p_res)

      call doLimit_TVD(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
          p_u, tstep, rafcstab%NEDGE, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

      ! Set specifier
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_BOUNDS)
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)
      

    case (AFCSTAB_FEMGP)
      
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
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
      call lsyssc_getbase_double(rmatrixMC, p_MC)
      call lsyssc_getbase_double(ru,   p_u)
      call lsyssc_getbase_double(ru0,  p_u0)
      call lsyssc_getbase_double(rres, p_res)

      if (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS) then
        call doLimit_GP(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_MC, p_u, p_u0, theta, tstep, rafcstab%NEDGE,&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_flux0, p_res)
      else
        call doLimit_TVD(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_u, tstep, rafcstab%NEDGE,&
            p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)
      end if
      
      ! Set specifier
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_BOUNDS)
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
      rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)

    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarTVD')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! The FEM-TVD limiting procedure
    
    subroutine doLimit_TVD(IverticesAtEdge, DcoefficientsAtEdge,&
                           u, tstep, NEDGE, pp, pm, qp, qm, rp, rm, flux, res)

      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: tstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      real(DP), dimension(:), intent(INOUT)            :: pp,pm
      real(DP), dimension(:), intent(INOUT)            :: qp,qm
      real(DP), dimension(:), intent(INOUT)            :: rp,rm
      real(DP), dimension(:), intent(INOUT)            :: flux
      real(DP), dimension(:), intent(INOUT)            :: res
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff

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
        pp(i) = pp(i)+max(0._DP,f_ij); pm(i) = pm(i)+min(0._DP,f_ij)
        
        ! Assemble Q's
        qp(i) = qp(i)+max(0._DP,-f_ij); qp(j) = qp(j)+max(0._DP, f_ij)
        qm(i) = qm(i)+min(0._DP,-f_ij); qm(j) = qm(j)+min(0._DP, f_ij)
      end do
      
      ! Apply the nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP, 1._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP, 1._DP)

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


    !**************************************************************
    ! The FEM-GP limiting procedure

    subroutine doLimit_GP(IverticesAtEdge, DcoefficientsAtEdge,&
                          MC, u, u0, theta, tstep, NEDGE,&
                          pp, pm, qp, qm, rp, rm, flux, flux0, res)

      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: MC
      real(DP), dimension(:), intent(IN)               :: u,u0
      real(DP), intent(IN)                             :: theta,tstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      real(DP), dimension(:), intent(INOUT)            :: pp,pm
      real(DP), dimension(:), intent(INOUT)            :: qp,qm
      real(DP), dimension(:), intent(INOUT)            :: rp,rm
      real(DP), dimension(:), intent(INOUT)            :: flux,flux0
      real(DP), dimension(:), intent(INOUT)            :: res
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,l_ij,l_ji,m_ij,p_ij,pf_ij,df_ij,q_ij,q_ji
      real(DP) :: diff,diff0,diff1

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
        diff  = tstep*(theta*diff1+(1._DP-theta)*diff0)
        
        ! Compute antidiffusive flux f_ij=min(0,p_ij)*(u_j-u_i)
        if (abs(diff) < SYS_EPSREAL) then
          p_ij = 0
          f_ij = 0
        else
          p_ij = max(0._DP,m_ij*(diff1-diff0)/diff+d_ij)
          f_ij = p_ij*diff
        end if
        
        ! Prelimit the antidiffusive flux F'_IJ=MIN(-P_IJ,L_JI)(U_I-U_J)
        pf_ij = min(p_ij,l_ji)*diff; flux0(iedge) = pf_ij
        
        ! Compute the remaining flux dF_IJ=F_IJ-F'_IJ
        df_ij = f_ij-pf_ij; flux(iedge) = df_ij
        
        ! Assemble P's accordingly
        pp(i) = pp(i)+max(0._DP,  f_ij); pm(i) = pm(i)+min(0._DP,  f_ij)
        pp(j) = pp(j)+max(0._DP,-df_ij); pm(j) = pm(j)+min(0._DP,-df_ij)
        
        q_ij = m_ij/tstep+l_ij; q_ji = m_ij/tstep+l_ji

        ! Assemble Q's
        qp(i) = qp(i)+q_ij*max(0._DP,-diff); qm(i) = qm(i)+q_ij*min(0._DP,-diff)
        qp(j) = qp(j)+q_ji*max(0._DP, diff); qm(j) = qm(j)+q_ji*min(0._DP, diff)
      end do

      ! Apply nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP, 1._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP, 1._DP)

      ! Apply limiter
      do iedge = 1, NEDGE

        ! Determine indices
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        
        ! Get precomputed fluxes
        pf_ij = flux0(iedge); df_ij = flux(iedge)
        
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
  end subroutine gfsc_buildResScalarTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildResBlockSymm(ru, dscale, rres, rafcstab)

!<description>
    ! This subroutine assembles the residual vector and applies stabilisation
    ! by means of symmetric flux limiting for diffusion operators.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN)    :: ru

    ! scaling parameter
    real(DP), intent(IN)               :: dscale
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT)     :: rafcstab
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
    type(t_vectorScalar), intent(IN)    :: ru

    ! scaling parameter
    real(DP), intent(IN)                :: dscale
!</input>

!<inputoutput>
    ! residual vector
    type(t_vectorScalar), intent(INOUT) :: rres

    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_VECIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer :: p_flux
    real(DP), dimension(:), pointer :: p_u,p_res
    

    ! Check if stabilisation is prepared
    if (rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_SYMMETRIC .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0    .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildResScalarSymm')
      call sys_halt()
    end if

    ! Check if vectors are compatible
    call lsyssc_isVectorCompatible(ru, rres)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
    call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
    call lsyssc_getbase_double(ru,   p_u)
    call lsyssc_getbase_double(rres, p_res)

    call doLimit_Symmetric(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
        p_u, dscale, rafcstab%NEDGE,&
        p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, p_flux, p_res)

    ! Set specifier
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_BOUNDS)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_ANTIDIFFUSION)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_LIMITER)
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_FLUXES)

  contains
    
    ! Here, the working routine follow
    
    !**************************************************************
    ! Perform symmetric flux limiting
    
    subroutine doLimit_Symmetric(IverticesAtEdge, DcoefficientsAtEdge, u, dscale,&
                                 NEDGE, pp, pm, qp, qm, rp, rm, flux, res)
      
      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: dscale
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      real(DP), dimension(:), intent(INOUT)            :: pp,pm
      real(DP), dimension(:), intent(INOUT)            :: qp,qm
      real(DP), dimension(:), intent(INOUT)            :: rp,rm
      real(DP), dimension(:), intent(INOUT)            :: flux
      real(DP), dimension(:), intent(INOUT)            :: res
      
      integer(PREC_MATIDX) :: iedge,ij
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,s_ij,diff
      
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
        pp(i) = pp(i)+max(0._DP, f_ij); pp(j) = pp(j)+max(0._DP,-f_ij)
        pm(i) = pm(i)+min(0._DP, f_ij); pm(j) = pm(j)+min(0._DP,-f_ij)
        
        ! Upper/lower bounds
        f_ij = -s_ij*diff
        qp(i) = qp(i)+max(0._DP, f_ij); qp(j) = qp(j)+max(0._DP,-f_ij)
        qm(i) = qm(i)+min(0._DP, f_ij); qm(j) = qm(j)+min(0._DP,-f_ij)
      end do
      
      ! Apply the nodal limiter
      rp = afcstab_limit( pp, qp, 0._DP, 1._DP)
      rm = afcstab_limit(-pm,-qm, 0._DP, 1._DP)
      
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
                                         hstep, bStabilise, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(IN)                :: ru
    
    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN)                            :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (ru%nblocks .ne. 1) then
      
      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvJacobianScalar(RcoeffMatrices, ru%RvectorBlock(1),&
          fcb_calcConvection, hstep, bStabilise, bclear, rmatrixJ)
      
    end if
  end subroutine gfsc_buildConvJacobianBlock

   !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianScalar(RcoeffMatrices, ru, fcb_calcConvection,&
                                          hstep, bStabilise, bclear, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(IN)               :: ru
    
    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(IN)                            :: bStabilise

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

    ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kld,p_Ksep
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kdiagonal
    integer(PREC_VECIDX), dimension(:), pointer   :: p_Kcol
    real(DP), dimension(:), pointer               :: p_Cx,p_Cy,p_Cz,p_J,p_u
    integer :: h_Ksep
    integer :: idim,ndim
    
    ! Check if all matrices are compatible
    call lsyssc_isMatrixCompatible(ru, rmatrixJ, .false.)
    ndim = size(RcoeffMatrices,1)
    do idim = 1, ndim
      call lsyssc_isMatrixCompatible(RcoeffMatrices(idim), rmatrixJ)
    end do

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrixJ)

    ! Set pointers
    call lsyssc_getbase_Kld(rmatrixJ,    p_Kld)
    call lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
    call lsyssc_getbase_double(rmatrixJ, p_J)
    call lsyssc_getbase_double(ru,       p_u)

    ! How many dimensions do we have?
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
    select case(rmatrixJ%cmatrixFormat)
    case(LSYSSC_MATRIX7)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrixJ%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

      ! Do we have to build the upwind Jacobian?
      if (bStabilise) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ,&
                               p_Cx, p_u, p_J)
        case (NDIM2D)
          call doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ,&
                               p_Cx, p_Cy, p_u, p_J)
        case (NDIM3D)
          call doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ,&
                               p_Cx, p_Cy, p_Cz, p_u, p_J)
        end select

      else

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ,&
                                 p_Cx, p_u, p_J)
        case (NDIM2D)
          call doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ,&
                                 p_Cx, p_Cy, p_u, p_J)
        case (NDIM3D)
          call doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep, rmatrixJ%NEQ,&
                                 p_Cx, p_Cy, p_Cz, p_u, p_J)
        end select

      end if

      ! Release diagonal separator
      call storage_free(h_Ksep)

      
    case(LSYSSC_MATRIX9)
      
      ! Set pointers
      call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
    
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrixJ%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bStabilise) then
        
        ! Build Jacobian
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rmatrixJ%NEQ, p_Cx, p_u, p_J)
        case (NDIM2D)
          call doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rmatrixJ%NEQ, p_Cx, p_Cy, p_u, p_J)
        case (NDIM3D)
          call doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                               rmatrixJ%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_J)
        end select
      
      else

        ! Build Jacobian
        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rmatrixJ%NEQ, p_Cx, p_u, p_J)
        case (NDIM2D)
          call doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rmatrixJ%NEQ, p_Cx, p_Cy, p_u, p_J)
        case (NDIM3D)
          call doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                 rmatrixJ%NEQ, p_Cx, p_Cy, p_Cz, p_u, p_J)
        end select

      end if

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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP)                    :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP)                    :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,Cz,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP)                    :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP)                    :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP)                    :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,Cz,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP)                    :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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
      
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP)                    :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j
      
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP)                    :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,Cz,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP)                    :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP)                    :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP)                    :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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

      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(INOUT) :: Ksep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,Cz,u
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP)                    :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer(PREC_MATIDX)        :: ii,ij,ji,jj
      integer(PREC_VECIDX)        :: i,j

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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)
          d_ij = max(-l_ij, 0._DP, -l_ji)
          
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

  subroutine gfsc_buildJacLinearBlockFCT(rmatrixMC, ru, theta, tstep, hstep,&
                                         bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)    :: rmatrixMC

    ! solution vector
    type(t_vectorBlock), intent(IN)     :: ru

    ! implicitness parameter
    real(DP), intent(IN)                :: theta

    ! time step size
    real(DP), intent(IN)                :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearBlockFCT')
      call sys_halt()

    else

      call gfsc_buildJacLinearScalarFCT(rmatrixMC, ru%RvectorBlock(1),&
          theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    end if
  end subroutine gfsc_buildJacLinearBlockFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearScalarFCT(rmatrixMC, ru, theta, tstep, hstep,&
                                          bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)    :: rmatrixMC

    ! solution vector
    type(t_vectorScalar), intent(IN)    :: ru

    ! implicitness parameter
    real(DP), intent(IN)                :: theta

    ! time step size
    real(DP), intent(IN)                :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_VECIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kld
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kdiagonal
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer               :: p_flux,p_flux0
    real(DP), dimension(:), pointer               :: p_u
    real(DP), dimension(:), pointer               :: p_MC,p_Jac
    
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_EDGEVALUES)    .eq. 0 .or.&
        iand(rafcstab%iSpec, AFCSTAB_FLUXES)        .eq. 0) then
      call output_line('Stabilisation does not provide required structures',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarFCT')
      call sys_halt()
    end if
    
    ! Check if matrices/vectors are compatible
    call lsyssc_isMatrixCompatible(ru, rmatrixMC, .false.)
    call lsyssc_isMatrixCompatible(ru, rmatrixJ,  .false.)
    call gfsc_isMatrixCompatible(rafcstab, rmatrixMC)
    call gfsc_isVectorCompatible(rafcstab, ru)
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrixJ)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    call lsyssc_getbase_double(rmatrixMC, p_MC)
    call lsyssc_getbase_double(rmatrixJ,  p_Jac)
    call lsyssc_getbase_double(ru,        p_u)
    

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT)
      
      ! What kind of matrix are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call doJacobian_implFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kld, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
        
      case(LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        call doJacobian_implFCT(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
            p_Kdiagonal, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
            rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
        
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
    ! Assemble the Jacobian matrix for semi-implicit FEM-FCT

    subroutine doJacobian_implFCT(IverticesAtEdge, DcoefficientsAtEdge,&
                                  Kdiagonal, MC, u, flux, flux0,&
                                  theta, tstep, hstep, NEDGE, bmass, Jac)

      integer(PREC_VECIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)   :: Kdiagonal
      real(DP), dimension(:), intent(IN)               :: MC
      real(DP), dimension(:), intent(IN)               :: flux,flux0
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: theta,tstep,hstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: Jac

      ! local variables
      integer(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      integer(PREC_VECIDX) :: i,j
      real(DP)             :: f_i,f_j,f_ij,d_ij,a_ij
      real(DP)             :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      if (bmass) then

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
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if
          
          ! Compute limited antidiffusive flux f(u_ij-h*e_j)
          f_j = a_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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

      else

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
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if
          
          ! Compute limited antidiffusive flux f(u_ij-h*e_j)
          f_j = a_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
      end if     
    end subroutine doJacobian_implFCT
  end subroutine gfsc_buildJacLinearScalarFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearBlockTVD(rmatrixMC, ru, ru0, theta, tstep, hstep,&
                                         bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)    :: rmatrixMC

    ! solution vector
    type(t_vectorBlock), intent(IN)     :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN)     :: ru0

    ! implicitness parameter
    real(DP), intent(IN)                :: theta

    ! time step size
    real(DP), intent(IN)                :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1 .or.&
        ru0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearBlockTVD')
      call sys_halt()

    else

      call gfsc_buildJacLinearScalarTVD(rmatrixMC,&
          ru%RvectorBlock(1), ru0%RvectorBlock(1),&
          theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    end if
  end subroutine gfsc_buildJacLinearBlockTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearScalarTVD(rmatrixMC, ru, ru0, theta, tstep, hstep,&
                                          bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)    :: rmatrixMC

    ! solution vector
    type(t_vectorScalar), intent(IN)    :: ru

    ! initial solution vector
    type(t_vectorScalar), intent(IN)    :: ru0

    ! implicitness parameter
    real(DP), intent(IN)                :: theta

    ! time step size
    real(DP), intent(IN)                :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                 :: bclear
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)      :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT) :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_VECIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    integer(PREC_VECIDX), dimension(:), pointer   :: p_IsuperdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdges
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kld,p_Ksep
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kdiagonal
    integer(PREC_VECIDX), dimension(:), pointer   :: p_Kcol
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer               :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer               :: p_flux,p_flux0
    real(DP), dimension(:), pointer               :: p_u,p_u0
    real(DP), dimension(:), pointer               :: p_MC,p_Jac
    integer :: h_Ksep
    logical :: bisExtended

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrixJ)

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case(AFCSTAB_FEMTVD)
      
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
      
      ! Check if matrices/vectors are compatible
      call lsyssc_isMatrixCompatible(ru, rmatrixJ, .false.)
      call gfsc_isVectorCompatible(rafcstab, ru)

      ! Check if subdiagonal edges need to be generated
      if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
          call afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      call afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_double(rmatrixJ, p_Jac)
      call lsyssc_getbase_double(ru,       p_u)

      ! What kind of matrix format are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)

        ! Create diagonal separator and increase it by one
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        call lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)
        call doJacobianMat79_TVD(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                 p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                 p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                 p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                 theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                 rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
        
        ! Free storage
        call storage_free(h_Ksep)

      case(LSYSSC_MATRIX9)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)
        call doJacobianMat79_TVD(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                 p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                 p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                 p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                 theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                 rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
        
        ! Free storage
        call storage_free(h_Ksep)
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarTVD')
        call sys_halt()
      end select
      

    case(AFCSTAB_FEMGP)

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

      ! Check if matrices/vectors are compatible
      call lsyssc_isMatrixCompatible(ru, rmatrixMC, .false.)
      call lsyssc_isMatrixCompatible(ru, rmatrixJ,  .false.)
      call gfsc_isMatrixCompatible(rafcstab, rmatrixMC)
      call gfsc_isVectorCompatible(rafcstab, ru)
      call lsyssc_isVectorCompatible(ru, ru0)
      
      ! Check if subdiagonal edges need to be generated
      if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
          call afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      call afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
      call lsyssc_getbase_double(rmatrixJ, p_Jac)
      call lsyssc_getbase_double(ru,       p_u)
      call lsyssc_getbase_double(ru0,      p_u0)

      ! What kind of matrix format are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        call lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)
        if (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS) then
          ! Set pointers
          call lsyssc_getbase_double(rmatrixMC, p_MC)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
          
          call doJacobianMat79_GP(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                  p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                  p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                  p_MC, p_u, p_u0, p_flux, p_flux0,&
                                  p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
                                  theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                  rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
          
        else
          call doJacobianMat79_TVD(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                   p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                   theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
        end if

        ! Free storage
        call storage_free(h_Ksep)

      case(LSYSSC_MATRIX9)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)
        if (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS) then
          ! Set pointers
          call lsyssc_getbase_double(rmatrixMC, p_MC)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

          call doJacobianMat79_GP(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                  p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                  p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                  p_MC, p_u, p_u0, p_flux, p_flux0,&
                                  p_pp, p_pm, p_qp, p_qm, p_rp, p_rm,&
                                  theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                  rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
        else
          call doJacobianMat79_TVD(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                   p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                   p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                   p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                   theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                   rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
        end if

        ! Free storage
        call storage_free(h_Ksep)
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearScalarTVD')
        call sys_halt()
      end select
      
    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
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
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_VECIDX), intent(IN)                  :: k
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      integer(PREC_MATIDX) :: ild,isep
      integer(PREC_VECIDX) :: l
      integer :: iloc
      
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
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_VECIDX), intent(IN)                  :: k
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      integer(PREC_MATIDX) :: ild,isep
      integer(PREC_VECIDX) :: l
      integer :: iloc
      
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
                                        u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
                                        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                :: u
      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7

      real(DP), dimension(:), intent(INOUT)             :: Jac
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      integer(PREC_VECIDX), dimension(NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc
           
      integer(PREC_MATIDX) :: ild,iedge
      integer(PREC_VECIDX) :: k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

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
      
      real(DP), dimension(:,:), intent(IN)                 :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                   :: u
      real(DP), dimension(:), intent(IN)                   :: pp,pm,qp,qm
      real(DP), intent(IN)                                 :: tstep,hstep
      integer(PREC_MATIDX), dimension(:,:), intent(IN)     :: IverticesAtEdge
      integer(PREC_MATIDX), intent(IN)                     :: iedge
      integer(PREC_VECIDX), intent(IN)                     :: k
      integer, intent(IN)                                  :: iloc

      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT)             :: pploc,pmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: qploc,qmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: fluxloc
      integer(PREC_VECIDX), dimension(:), intent(INOUT)    :: Kloc

      ! local variables
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff,dsign
      integer  :: iperturb

      
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
        qploc(:,iloc) = qp(j)-max(0._DP, f_ij)
        qmloc(:,iloc) = qm(j)-min(0._DP, f_ij)

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Compute perturbed antidiffusive flux
          f_ij = min(d_ij,l_ji)*(diff+tstep*dsign*hstep)
          fluxloc(iperturb,iloc) = f_ij

          ! For node k which is the upwind node
          pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP, f_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP, f_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,-f_ij)
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,-f_ij)
          
          ! For node l opposite to k which is the downwind node
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,f_ij)
        end do

      else

        ! Store global node number of the opposite node
        Kloc(iloc) = i
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0._DP,-f_ij)
        qmloc(:,iloc) = qm(i)-min(0._DP,-f_ij)

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Compute perturbed antidiffusive flux
          f_ij = min(d_ij,l_ji)*(diff-tstep*dsign*hstep)
          fluxloc(iperturb,iloc) = f_ij
          
          ! For node k which is the downwind node
          qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,f_ij)
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,f_ij)
          
          ! For node l opposite to k
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,-f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,-f_ij)
        end do
      end if
    end subroutine updateJacobianMat79_TVD


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    pure subroutine assembleJacobianMat79_TVD(IverticesAtEdge, Kdiagonal,&
                                              flux, Kloc, rploc, rmloc, fluxloc,&
                                              hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:,0:), intent(IN)             :: rploc,rmloc
      real(DP), dimension(:,0:), intent(IN)             :: fluxloc
      real(DP), intent(IN)                              :: hstep
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kloc
      integer(PREC_MATIDX), intent(IN)                  :: iedge
      integer(PREC_VECIDX), intent(IN)                  :: k,l
      integer, intent(IN)                               :: iloc
      logical, intent(IN)                               :: bisExtended

      real(DP), dimension(:), intent(INOUT)             :: Jac
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      

      ! local variables
      integer(PREC_MATIDX) :: ik,jk
      integer(PREC_VECIDX) :: i,j,m
      real(DP)             :: f_ij,df_ij
      integer              :: iperturb
      
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

    
    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                  IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                  DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                  MC, u, u0, flux, flux0, pp, pm, qp, qm, rp, rm,&
                                  theta, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
                                  bisExtended, bisMat7, Ksep, Jac)
    
      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                :: MC
      real(DP), dimension(:), intent(IN)                :: u,u0
      real(DP), dimension(:), intent(IN)                :: flux,flux0
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm,rp,rm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7

      real(DP), dimension(:), intent(INOUT)             :: Jac
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      integer(PREC_VECIDX), dimension(NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc,fluxloc0

      integer(PREC_MATIDX) :: ild,iedge
      integer(PREC_VECIDX) :: k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


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
      
      real(DP), dimension(:,:), intent(IN)                 :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                   :: MC
      real(DP), dimension(:), intent(IN)                   :: u,u0
      real(DP), dimension(:), intent(IN)                   :: flux,flux0
      real(DP), dimension(:), intent(IN)                   :: pp,pm,qp,qm
      real(DP), intent(IN)                                 :: theta,tstep,hstep
      integer(PREC_MATIDX), dimension(:,:), intent(IN)     :: IverticesAtEdge
      integer(PREC_MATIDX), intent(IN)                     :: iedge
      integer(PREC_VECIDX), intent(IN)                     :: k
      integer, intent(IN)                                  :: iloc

      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT)             :: pploc,pmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: qploc,qmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: fluxloc,fluxloc0
      integer(PREC_VECIDX), dimension(:), intent(INOUT)    :: Kloc

      ! local variables
      integer(PREC_VECIDX) :: i,j,ij
      real(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji
      real(DP) :: diff,diff1,diff0,dsign
      integer  :: iperturb

      
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
        p_ij = max(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
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
        pploc(:,iloc) = pp(j)-max(0._DP,-df_ij)
        pmloc(:,iloc) = pm(j)-min(0._DP,-df_ij)
        qploc(:,iloc) = qp(j)-max(0._DP, diff)*q_ji
        qmloc(:,iloc) = qm(j)-min(0._DP, diff)*q_ji

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
            p_ij = max(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          end if
          
          ! Prelimit the antidiffusive flux
          pf_ij = min(p_ij, l_ji)*diff
          fluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij
          
          ! For node k which is the upwind node
          pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP, f_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP, f_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,-diff)*q_ij
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,-diff)*q_ij
          
          ! For node l opposite to k which is the downwind node
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP,-df_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP,-df_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP, diff)*q_ji
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP, diff)*q_ji
        end do

      else

        ! Store global node number of the opposite node
        Kloc(iloc) = i

        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0._DP,-diff)*q_ij
        qmloc(:,iloc) = qm(i)-min(0._DP,-diff)*q_ij

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
            p_ij = max(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          end if
          
          ! Prelimit the antidiffusive flux
          pf_ij = min(p_ij, l_ji)*diff
          fluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          fluxloc(iperturb,iloc) = df_ij

          ! For node k which is the downwind node
          pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP,-df_ij)
          pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP,-df_ij)
          qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP, diff)*q_ji
          qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP, diff)*q_ji
          
          ! For node l opposite to k
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,-diff)*q_ij
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,-diff)*q_ij
        end do
      end if
     
        
!!$      DO iperturb = 1, 2
!!$        
!!$        ! Compute correct sign of perturbation
!!$        dsign = -2*iperturb+3
!!$        
!!$        ! The velocity is assumed to be linear.
!!$        ! Hence the orientation convention for the edge ij is preserved

!!$        ! Save oriented node numbers
!!$        Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
        
!!$        ! Update solution difference
!!$        diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)
!!$        
!!$        ! Update total solution difference
!!$        diff = tstep*(theta*diff1+(1-theta)*diff0)
!!$        
!!$        ! Compute antidiffusive flux
!!$        IF (ABS(diff) < SYS_EPSREAL) THEN
!!$          p_ij = 0
!!$          f_ij = 0
!!$        ELSE
!!$          p_ij = MAX(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
!!$          f_ij = p_ij*diff
!!$        END IF
!!$        
!!$        ! Prelimit the antidiffusive flux
!!$        pf_ij = MIN(p_ij, l_ji)*diff
!!$        fluxloc0(iperturb,iloc) = pf_ij
!!$        
!!$        ! Compute the remaining flux
!!$        df_ij = f_ij-pf_ij
!!$        fluxloc(iperturb,iloc) = df_ij
!!$        
!!$        IF (i .EQ. k) THEN
!!$
!!$          ! For node k which is the upwind node
!!$          pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP, f_ij)
!!$          pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP, f_ij)
!!$          qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP,-diff)*q_ij
!!$          qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP,-diff)*q_ij
!!$          
!!$          ! For node l opposite to k which is the downwind node
!!$          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP,-df_ij)
!!$          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP,-df_ij)
!!$          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP, diff)*q_ji
!!$          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP, diff)*q_ji
!!$
!!$        ELSE
!!$
!!$          ! For node k which is the downwind node
!!$          pploc(iperturb,0) = pploc(iperturb,0)+MAX(0._DP,-df_ij)
!!$          pmloc(iperturb,0) = pmloc(iperturb,0)+MIN(0._DP,-df_ij)
!!$          qploc(iperturb,0) = qploc(iperturb,0)+MAX(0._DP, diff)*q_ji
!!$          qmloc(iperturb,0) = qmloc(iperturb,0)+MIN(0._DP, diff)*q_ji
!!$          
!!$          ! For node l opposite to k
!!$          pploc(iperturb,iloc) = pploc(iperturb,iloc)+MAX(0._DP, f_ij)
!!$          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+MIN(0._DP, f_ij)
!!$          qploc(iperturb,iloc) = qploc(iperturb,iloc)+MAX(0._DP,-diff)*q_ij
!!$          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+MIN(0._DP,-diff)*q_ij
!!$
!!$        END IF        
!!$      END DO
    end subroutine updateJacobianMat79_GP


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_GP(IverticesAtEdge, Kdiagonal,&
                                        flux, flux0, rp, rm, Kloc,&
                                        rploc, rmloc, fluxloc, fluxloc0,&
                                        hstep, iedge, iloc, k, l,&
                                        bisExtended, Ksep, Jac)

      real(DP), dimension(:), intent(IN)                :: flux,flux0
      real(DP), dimension(:), intent(IN)                :: rp,rm
      real(DP), dimension(:,0:), intent(IN)             :: rploc,rmloc
      real(DP), dimension(:,0:), intent(IN)             :: fluxloc,fluxloc0
      real(DP), intent(IN)                              :: hstep
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kloc
      integer(PREC_MATIDX), intent(IN)                  :: iedge
      integer, intent(IN)                               :: iloc
      integer(PREC_VECIDX), intent(IN)                  :: k,l
      logical, intent(IN)                               :: bisExtended

      real(DP), dimension(:), intent(INOUT)             :: Jac
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep


      ! local variables
      integer(PREC_MATIDX) :: ik,jk
      integer(PREC_VECIDX) :: i,j,m
      real(DP)             :: f_ij,pf_ij,df_ij
      integer              :: iperturb
      
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

!!$          ! Adjust edge orientation
!!$          i = Kloc(2*iperturb,iloc)
!!$          j = Kloc(2*iperturb+1,iloc)
          
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
  end subroutine gfsc_buildJacLinearScalarTVD

  !*****************************************************************************
  
!<subroutine>

  subroutine gfsc_buildJacobianBlockFCT(RcoeffMatrices, rmatrixMC, ru, &
                                        fcb_calcConvection, theta, tstep, hstep,&
                                        bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)               :: rmatrixMC

    ! solution vector
    type(t_vectorBlock), intent(IN)                :: ru

    ! implicitness parameter
    real(DP), intent(IN)                           :: theta

    ! time step size
    real(DP), intent(IN)                           :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)                 :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockFCT')
      call sys_halt()

    else

      call gfsc_buildJacobianScalarFCT(RcoeffMatrices, rmatrixMC, ru%RvectorBlock(1),&
          fcb_calcConvection, theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    end if
  end subroutine gfsc_buildJacobianBlockFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarFCT(RcoeffMatrices, rmatrixMC, ru,&
                                         fcb_calcConvection, theta, tstep, hstep,&
                                         bclear, rafcstab, rmatrixJ)

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
    type(t_matrixScalar), intent(IN)               :: rmatrixMC

    ! solution vector
    type(t_vectorScalar), intent(IN)               :: ru

    ! implicitness parameter
    real(DP), intent(IN)                           :: theta

    ! time step size
    real(DP), intent(IN)                           :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)                 :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_MATIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kld
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kdiagonal
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer               :: p_flux,p_flux0
    real(DP), dimension(:), pointer               :: p_u
    real(DP), dimension(:), pointer               :: p_Cx,p_Cy,p_Cz,p_MC,p_Jac
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
    if (bclear) call lsyssc_clearMatrix(rmatrixJ)

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(1), p_flux)
    call lsyssc_getbase_double(rafcstab%RedgeVectors(2), p_flux0)
    call lsyssc_getbase_double(rmatrixMC, p_MC)
    call lsyssc_getbase_double(rmatrixJ,  p_Jac)
    call lsyssc_getbase_double(ru,        p_u)
    
    ! Set spatial dimensions
    ndim = size(RcoeffMatrices,1)

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case (AFCSTAB_FEMFCT)
      
      ! What kind of matrix format are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        
        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          
          call doJacobian_implFCT_1D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_Kld, p_Cx, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
              rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
          
        case (NDIM2D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
          
          call doJacobian_implFCT_2D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_Kld, p_Cx, p_Cy, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
              rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
          
        case (NDIM3D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)        
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
          call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)
          
          call doJacobian_implFCT_3D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_Kld, p_Cx, p_Cy, p_Cz, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
              rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
        end select
        
        
      case(LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          
          call doJacobian_implFCT_1D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Cx, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
              rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
          
        case (NDIM2D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
          
          call doJacobian_implFCT_2D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Cx, p_Cy, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
              rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
          
        case (NDIM3D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)        
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
          call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)
          
          call doJacobian_implFCT_3D(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
              p_Kdiagonal, p_Cx, p_Cy, p_Cz, p_MC, p_u, p_flux, p_flux0, theta, tstep, hstep,&
              rafcstab%NEDGE, (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS), p_Jac)
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
    ! Assemble the Jacobian matrix for FEM-FCT in 1D
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCT_1D(IverticesAtEdge, DcoefficientsAtEdge,&
                                     Kdiagonal, Cx, MC, u, flux, flux0,&
                                     theta, tstep, hstep, NEDGE, bmass, Jac)

      integer(PREC_MATIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)   :: Kdiagonal
      real(DP), dimension(:), intent(IN)               :: Cx,MC
      real(DP), dimension(:), intent(IN)               :: flux,flux0
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: theta,tstep,hstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: Jac

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      integer(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      integer(PREC_VECIDX) :: i,j
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji
      real(DP) :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      if (bmass) then

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
          a_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if

          
          ! Compute perturbed velocity
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          b_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
          a_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
          f_i = a_ij*diff_j+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          b_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
          end if
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        end do

      else

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
          a_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          b_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
        
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
          a_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
          f_i = a_ij*diff_j+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          b_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
          end if
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        end do
        
      end if
    end subroutine doJacobian_implFCT_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 2D
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCT_2D(IverticesAtEdge, DcoefficientsAtEdge,&
                                     Kdiagonal, Cx, Cy, MC, u, flux, flux0,&
                                     theta, tstep, hstep, NEDGE, bmass, Jac)

      integer(PREC_MATIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)   :: Kdiagonal
      real(DP), dimension(:), intent(IN)               :: Cx,Cy,MC
      real(DP), dimension(:), intent(IN)               :: flux,flux0
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: theta,tstep,hstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: Jac

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      integer(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      integer(PREC_VECIDX) :: i,j
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji
      real(DP) :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      if (bmass) then

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
          a_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if

          
          ! Compute perturbed velocity
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          b_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
          a_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
          f_i = a_ij*diff_j+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          b_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
          end if
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        end do

      else

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
          a_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          b_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
        
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
          a_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
          f_i = a_ij*diff_j+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          b_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
          end if
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        end do
        
      end if
    end subroutine doJacobian_implFCT_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 3D
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCT_3D(IverticesAtEdge, DcoefficientsAtEdge,&
                                     Kdiagonal, Cx, Cy, Cz, MC, u, flux, flux0,&
                                     theta, tstep, hstep, NEDGE, bmass, Jac)

      integer(PREC_MATIDX), dimension(:,:), intent(IN) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)             :: DcoefficientsAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)   :: Kdiagonal
      real(DP), dimension(:), intent(IN)               :: Cx,Cy,Cz,MC
      real(DP), dimension(:), intent(IN)               :: flux,flux0
      real(DP), dimension(:), intent(IN)               :: u
      real(DP), intent(IN)                             :: theta,tstep,hstep
      integer(PREC_MATIDX), intent(IN)                 :: NEDGE
      logical, intent(IN)                              :: bmass
      real(DP), dimension(:), intent(INOUT)            :: Jac

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      integer(PREC_MATIDX) :: iedge,ij,ji,ii,jj
      integer(PREC_VECIDX) :: i,j
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji
      real(DP) :: diff,diff_i,diff_j

      ! Should we apply the consistent mass matrix?
      if (bmass) then

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
          a_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if

          
          ! Compute perturbed velocity
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          b_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
          a_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j) 
          f_i = a_ij*diff_j+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          b_ij = MC(ij)+theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
          end if
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        end do

      else

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
          a_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_i)
          f_i = a_ij*diff_i+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i)-hstep, u(j), C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_i)
          b_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
        
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_j+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
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
          a_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij+h*e_j)
          f_i = a_ij*diff_j+flux0(iedge)
          if (f_i > 0.0_DP) then
            f_i = min(f_i, max(flux(iedge), 0._DP))
          else
            f_i = max(f_i, min(flux(iedge), 0._DP))
          end if


          ! Compute perturbed velocity
          call fcb_calcConvection(u(i), u(j)-hstep, C_ij, C_ji, i, j, l_ij, l_ji)

          ! Compute perturbed coefficient b_ij(u-hstep*e_j)
          b_ij = theta*tstep*max(-l_ij, 0._DP, -l_ji)
          
          ! Compute and limit raw antidiffusive flux f(u_ij-h*e_j)
          f_j = b_ij*diff_i+flux0(iedge)
          if (f_j > 0.0_DP) then
            f_j = min(f_j, max(flux(iedge), 0._DP))
          else
            f_j = max(f_j, min(flux(iedge), 0._DP))
          end if
          

          ! Compute divided differences of fluxes
          f_ij = 0.5_DP*(f_i-f_j)/hstep
          
          ! Apply j-th column
          Jac(ij) = Jac(ij)-f_ij
          Jac(jj) = Jac(jj)+f_ij
        end do
        
      end if
    end subroutine doJacobian_implFCT_3D
  end subroutine gfsc_buildJacobianScalarFCT

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianBlockTVD(RcoeffMatrices, rmatrixMC, ru, ru0,&
                                        fcb_calcConvection, theta, tstep, hstep,&
                                        bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(IN) :: RcoeffMatrices

    ! consistent mass matrix
    type(t_matrixScalar), intent(IN)               :: rmatrixMC

    ! solution vector
    type(t_vectorBlock), intent(IN)                :: ru

    ! initial solution vector
    type(t_vectorBlock), intent(IN)                :: ru0

    ! implicitness parameter
    real(DP), intent(IN)                           :: theta

    ! time step size
    real(DP), intent(IN)                           :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)                 :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1 .or.&
        ru0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockTVD')
      call sys_halt()

    else

      call gfsc_buildJacobianScalarTVD(RcoeffMatrices, rmatrixMC,&
          ru%RvectorBlock(1), ru0%RvectorBlock(1), fcb_calcConvection,&
          theta, tstep, hstep, bclear, rafcstab, rmatrixJ)

    end if
  end subroutine gfsc_buildJacobianBlockTVD
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarTVD(RcoeffMatrices, rmatrixMC, ru, ru0,&
                                         fcb_calcConvection, theta, tstep, hstep,&
                                         bclear, rafcstab, rmatrixJ)

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
    type(t_matrixScalar), intent(IN)               :: rmatrixMC

    ! solution vector
    type(t_vectorScalar), intent(IN)               :: ru

    ! initial solution vector
    type(t_vectorScalar), intent(IN)               :: ru0

    ! implicitness parameter
    real(DP), intent(IN)                           :: theta

    ! time step size
    real(DP), intent(IN)                           :: tstep

    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)                 :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_MATIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    integer(PREC_VECIDX), dimension(:), pointer   :: p_IsuperdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdges
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kld,p_Ksep
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kdiagonal
    integer(PREC_VECIDX), dimension(:), pointer   :: p_Kcol
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer               :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer               :: p_flux,p_flux0
    real(DP), dimension(:), pointer               :: p_u,p_u0
    real(DP), dimension(:), pointer               :: p_Cx,p_Cy,p_Cz,p_MC,p_Jac
    integer :: h_Ksep
    integer :: ndim
    logical :: bisExtended
    
    ! Set spatial dimensions
    ndim = size(RcoeffMatrices,1)
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrixJ)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)
      
    case(AFCSTAB_FEMTVD)
      
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

      ! Check if subdiagonal edges need to be generated
      if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
          call afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      call afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_double(rmatrixJ, p_Jac)
      call lsyssc_getbase_double(ru,       p_u)

      ! What kind of matrix format are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        call lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)

        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

          call doJacobianMat79_TVD_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                      p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                      p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                      p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                      theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                      rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)

        case (NDIM2D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

          call doJacobianMat79_TVD_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                      p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                      p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                      p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                      theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                      rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
          
        case (NDIM3D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
          call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

          call doJacobianMat79_TVD_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                      p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                      p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                      p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                      theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                      rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
        end select
               
        ! Free storage
        call storage_free(h_Ksep)

        
      case(LSYSSC_MATRIX9)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)

        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          
          call doJacobianMat79_TVD_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                      p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                      p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                      p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                      theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                      rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)

        case (NDIM2D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

          call doJacobianMat79_TVD_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                      p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                      p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                      p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                      theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                      rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)

        case (NDIM3D)
          call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
          call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
          call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

          call doJacobianMat79_TVD_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                      p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                      p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                      p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                      theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                      rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
        end select
        
        ! Free storage
        call storage_free(h_Ksep)

        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarTVD')
        call sys_halt()
      end select

      
    case(AFCSTAB_FEMGP)
      
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

      ! Check if subdiagonal edges need to be generated
      if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
          call afcstab_generateSubdiagEdges(rafcstab)

      ! Set pointers
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      call afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(2),  p_flux0)
      call lsyssc_getbase_double(rmatrixJ, p_Jac)
      call lsyssc_getbase_double(ru,       p_u)
      call lsyssc_getbase_double(ru0,      p_u0)

      ! What kind of matrix format are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        call lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)
        if (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrixMC, p_MC)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)

            call doJacobianMat79_GP_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                       p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                       p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                       p_Cx, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                       p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                                       rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
                                       bisExtended, .true., p_Ksep, p_Jac)

          case (NDIM2D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

            call doJacobianMat79_GP_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                       p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                       p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                       p_Cx, p_Cy, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                       p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                                       rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
                                       bisExtended, .true., p_Ksep, p_Jac)
            
          case (NDIM3D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
            call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

            call doJacobianMat79_GP_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                       p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                       p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                       p_Cx, p_Cy, p_Cz, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                       p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                                       rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
                                       bisExtended, .true., p_Ksep, p_Jac)
          end select

        else

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            
            call doJacobianMat79_TVD_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                        p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                        p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                        p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                        theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                        rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
            
          case (NDIM2D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

            call doJacobianMat79_TVD_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                        p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                        p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                        p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                        theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                        rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
            
          case (NDIM3D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
            call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

            call doJacobianMat79_TVD_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                        p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                        p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
                                        p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                        theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                        rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
          end select

        end if

        ! Free storage
        call storage_free(h_Ksep)
        
      case(LSYSSC_MATRIX9)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        call lsyssc_getbase_Kcol(rmatrixJ, p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrixJ, p_Kdiagonal)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)
        if (rafcstab%imass .eq. AFCSTAB_CONSISTENTMASS) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrixMC, p_MC)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
          call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            
            call doJacobianMat79_GP_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                       p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                       p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                       p_Cx, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                       p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                                       rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
                                       bisExtended, .false., p_Ksep, p_Jac)

          case (NDIM2D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)

            call doJacobianMat79_GP_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                       p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                       p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                       p_Cx, p_Cy, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                       p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                                       rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
                                       bisExtended, .false., p_Ksep, p_Jac)
            
          case (NDIM3D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
            call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

            call doJacobianMat79_GP_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                       p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                       p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                       p_Cx, p_Cy, p_Cz, p_MC, p_u, p_u0, p_flux, p_flux0,&
                                       p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, theta, tstep, hstep,&
                                       rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
                                       bisExtended, .false., p_Ksep, p_Jac)
          end select

        else

          ! How many dimensions do we have?
          select case(ndim)
          case (NDIM1D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            
            call doJacobianMat79_TVD_1D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                        p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                        p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                        p_Cx, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                        theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                        rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
            
          case (NDIM2D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
            
            call doJacobianMat79_TVD_2D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                        p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                        p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                        p_Cx, p_Cy, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                        theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                        rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
            
          case (NDIM3D)
            call lsyssc_getbase_double(RcoeffMatrices(1), p_Cx)
            call lsyssc_getbase_double(RcoeffMatrices(2), p_Cy)
            call lsyssc_getbase_double(RcoeffMatrices(3), p_Cz)

            call doJacobianMat79_TVD_3D(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
                                        p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges,&
                                        p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
                                        p_Cx, p_Cy, p_Cz, p_u, p_flux, p_pp, p_pm, p_qp, p_qm,&
                                        theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
                                        rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
          end select

        end if

        ! Free storage
        call storage_free(h_Ksep)
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarTVD')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarTVD')
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
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_VECIDX), intent(IN)                  :: k
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      integer(PREC_MATIDX) :: ild,isep
      integer(PREC_VECIDX) :: l
      integer :: iloc
      
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
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_VECIDX), intent(IN)                  :: k
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      integer(PREC_MATIDX) :: ild,isep
      integer(PREC_VECIDX) :: l
      integer :: iloc
      
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
                                      Cx, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
                                      NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                :: Cx
      real(DP), dimension(:), intent(IN)                :: u
      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7

      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc
      real(DP), dimension(NDIM1D)      :: c_ij, c_ji

      integer(PREC_MATIDX) :: ij,ji,ild,iedge
      integer(PREC_VECIDX) :: i,j,k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

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
                                      Cx, Cy, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
                                      NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                :: Cx,Cy
      real(DP), dimension(:), intent(IN)                :: u
      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc
      real(DP), dimension(NDIM2D)      :: c_ij, c_ji

      integer(PREC_MATIDX) :: ij,ji,ild,iedge
      integer(PREC_VECIDX) :: i,j,k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

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
                                      Cx, Cy, Cz, u, flux, pp, pm, qp, qm, theta, tstep, hstep,&
                                      NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,Cz
      real(DP), dimension(:), intent(IN)                :: u
      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc
      real(DP), dimension(NDIM3D)      :: c_ij, c_ji

      integer(PREC_MATIDX) :: ij,ji,ild,iedge
      integer(PREC_VECIDX) :: i,j,k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

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
      
      real(DP), dimension(:,:), intent(IN)                 :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                   :: u
      real(DP), dimension(:), intent(IN)                   :: pp,pm,qp,qm
      real(DP), dimension(:), intent(IN)                   :: C_ij,C_ji      
      real(DP), intent(IN)                                 :: tstep,hstep
      integer(PREC_MATIDX), intent(IN)                     :: iedge
      integer(PREC_VECIDX), intent(IN)                     :: i,j
      integer(PREC_MATIDX), intent(IN)                     :: ij,ji
      integer, intent(IN)                                  :: iloc
      integer(PREC_VECIDX), intent(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT)             :: pploc,pmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: qploc,qmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: fluxloc
      integer(PREC_VECIDX), dimension(:,:), intent(INOUT)  :: Kloc

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
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)
        pmloc(:,iloc) = pm(j)
        qploc(:,iloc) = qp(j)-max(0._DP, f_ij)
        qmloc(:,iloc) = qm(j)-min(0._DP, f_ij)

      else

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0._DP,-f_ij)
        qmloc(:,iloc) = qm(i)-min(0._DP,-f_ij)
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
        d_ij = max(-l_ij, 0._DP, -l_ji)

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
            pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,-f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,-f_ij)
            
            ! For node l opposite to k which is the downwind node
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,f_ij)

          else

            ! For node k which is the downwind node
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,f_ij)
            
            ! For node l opposite to k
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,-f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,-f_ij)

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
            pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,-f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,-f_ij)
            
            ! For node "l" opposite to k which is the downwind node
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,f_ij)

          else

            ! For node k which is the downwind node
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,f_ij)
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,f_ij)
            
            ! For node "l" opposite to k which is the upwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,-f_ij)
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,-f_ij)
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

      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:,0:), intent(IN)             :: rploc,rmloc
      real(DP), dimension(:,0:), intent(IN)             :: fluxloc
      real(DP), intent(IN)                              :: hstep
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:,:), intent(IN)  :: Kloc
      integer(PREC_MATIDX), intent(IN)                  :: iedge
      integer, intent(IN)                               :: iloc
      integer(PREC_VECIDX), intent(IN)                  :: k,l
      logical, intent(IN)                               :: bisExtended

      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac

      ! local variables
      integer(PREC_MATIDX) :: ik,jk
      integer(PREC_VECIDX) :: i,j,m
      real(DP)             :: f_ij
      integer              :: iperturb
      
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


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 1D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP_1D(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                                     IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                                     DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
                                     Cx, MC, u, u0, flux, flux0,&
                                     pp, pm, qp, qm, rp, rm,&
                                     theta, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
                                     bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                :: Cx,MC
      real(DP), dimension(:), intent(IN)                :: u,u0
      real(DP), dimension(:), intent(IN)                :: flux,flux0
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm,rp,rm
      real(DP), intent(IN)                              :: theta,tstep,hstep  
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc,fluxloc0
      real(DP), dimension(NDIM1D)      :: c_ij,c_ji
      
      integer(PREC_MATIDX) :: ij,ji,ild,iedge
      integer(PREC_VECIDX) :: i,j,k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


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
                                     Cx, Cy, MC, u, u0, flux, flux0,&
                                     pp, pm, qp, qm, rp, rm,&
                                     theta, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
                                     bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge    
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,MC
      real(DP), dimension(:), intent(IN)                :: u,u0
      real(DP), dimension(:), intent(IN)                :: flux,flux0
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm,rp,rm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc,fluxloc0
      real(DP), dimension(NDIM2D)      :: c_ij,c_ji
      
      integer(PREC_MATIDX) :: ij,ji,ild,iedge
      integer(PREC_VECIDX) :: i,j,k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


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
                                     Cx, Cy, Cz, MC, u, u0, flux, flux0,&
                                     pp, pm, qp, qm, rp, rm,&
                                     theta, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
                                     bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge    
      real(DP), dimension(:), intent(IN)                :: Cx,Cy,Cz,MC
      real(DP), dimension(:), intent(IN)                :: u,u0
      real(DP), dimension(:), intent(IN)                :: flux,flux0
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm,rp,rm
      real(DP), intent(IN)                              :: theta,tstep,hstep
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      logical, intent(IN)                               :: bisMat7
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc,fluxloc0
      real(DP), dimension(NDIM3D)      :: c_ij,c_ji
      
      integer(PREC_MATIDX) :: ij,ji,ild,iedge
      integer(PREC_VECIDX) :: i,j,k,l
      integer :: iloc,nloc

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
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)


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
      
      real(DP), dimension(:,:), intent(IN)                 :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                   :: MC
      real(DP), dimension(:), intent(IN)                   :: u,u0
      real(DP), dimension(:), intent(IN)                   :: flux,flux0
      real(DP), dimension(:), intent(IN)                   :: pp,pm,qp,qm
      real(DP), dimension(:), intent(IN)                   :: C_ij,C_ji
      real(DP), intent(IN)                                 :: theta,tstep,hstep
      integer(PREC_MATIDX), intent(IN)                     :: iedge
      integer(PREC_VECIDX), intent(IN)                     :: i,j
      integer(PREC_MATIDX), intent(IN)                     :: ij,ji
      integer, intent(IN)                                  :: iloc
      integer(PREC_VECIDX), intent(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT)             :: pploc,pmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: qploc,qmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: fluxloc,fluxloc0
      integer(PREC_VECIDX), dimension(:,:), intent(INOUT)  :: Kloc

      ! local variables
      real(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji
      real(DP) :: diff,diff1,diff0,hstep_ik,hstep_jk,dsign
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
        p_ij = max(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
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
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        pploc(:,iloc) = pp(j)-max(0._DP,-df_ij)
        pmloc(:,iloc) = pm(j)-min(0._DP,-df_ij)
        qploc(:,iloc) = qp(j)-max(0._DP, diff)*q_ji
        qmloc(:,iloc) = qm(j)-min(0._DP, diff)*q_ji

      else

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        pploc(:,iloc) = pp(i)-max(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0._DP, f_ij)
        qploc(:,iloc) = qp(i)-max(0._DP,-diff)*q_ij
        qmloc(:,iloc) = qm(i)-min(0._DP,-diff)*q_ij
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
        d_ij = max(-l_ij, 0._DP, -l_ji)

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
            p_ij = max(0._DP,m_ij*(diff1-diff0)/diff+d_ij)
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
            pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,-diff)*q_ij
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,-diff)*q_ij
            
            ! For node l opposite to k which is the downwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP,-df_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP,-df_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP, diff)*q_ji
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP, diff)*q_ji

          else

            ! For node k which is the downwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP,-df_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP,-df_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP, diff)*q_ji
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP, diff)*q_ji
            
            ! For node l opposite to k
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,-diff)*q_ij
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,-diff)*q_ij

          end if
          
        else
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! Update solution difference
          diff1 = u(i)-u(j)+dsign*(hstep_ik-hstep_jk)
          
          ! Update total solution difference
          diff = tstep*(theta*diff1+(1._DP-theta)*diff0)

          ! Compute antidiffusive flux
          if (abs(diff) < SYS_EPSREAL) then
            p_ij = 0
            f_ij = 0
          else
            p_ij = max(0._DP, m_ij*(diff1-diff0)/diff+d_ij)
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
            pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP, f_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP, f_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP, diff)*q_ij
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP, diff)*q_ij
                       
            ! For node "l" opposite to k which is the downwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP,-df_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP,-df_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP,-diff)*q_ji
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP,-diff)*q_ji

          else

            ! For node k which is the downwind node
            pploc(iperturb,0) = pploc(iperturb,0)+max(0._DP,-df_ij)
            pmloc(iperturb,0) = pmloc(iperturb,0)+min(0._DP,-df_ij)
            qploc(iperturb,0) = qploc(iperturb,0)+max(0._DP,-diff)*q_ji
            qmloc(iperturb,0) = qmloc(iperturb,0)+min(0._DP,-diff)*q_ji
            
            ! For node "l" opposite to k which is the upwind node
            pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
            pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
            qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP, diff)*q_ij
            qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP, diff)*q_ij

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

      real(DP), dimension(:), intent(IN)                :: flux,flux0
      real(DP), dimension(:), intent(IN)                :: rp,rm
      real(DP), dimension(:,0:), intent(IN)             :: rploc,rmloc
      real(DP), dimension(:,0:), intent(IN)             :: fluxloc,fluxloc0
      real(DP), intent(IN)                              :: hstep
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:,:), intent(IN)  :: Kloc
      integer(PREC_MATIDX), intent(IN)                  :: iedge
      integer, intent(IN)                               :: iloc
      integer(PREC_VECIDX), intent(IN)                  :: k,l
      logical, intent(IN)                               :: bisExtended

      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac

      ! local variables
      integer(PREC_MATIDX) :: ik,jk
      integer(PREC_VECIDX) :: i,j,m
      real(DP)             :: f_ij,pf_ij,df_ij
      integer              :: iperturb
      
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
  end subroutine gfsc_buildJacobianScalarTVD

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianBlockSymm(ru, dscale, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete diffusion operator for a scalar convection equation.
    ! Note that this routine serves as a wrapper for block vectors. If there
    ! is only one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN)                :: ru

    ! scaling parameter
    real(DP), intent(IN)                           :: dscale

    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)                 :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    if (ru%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlockSymm')
      call sys_halt()

    else

      call gfsc_buildJacobianScalarSymm(ru%RvectorBlock(1), dscale,&
          hstep, bclear, rafcstab, rmatrixJ)

    end if
  end subroutine gfsc_buildJacobianBlockSymm

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalarSymm(ru, dscale, hstep, bclear, rafcstab, rmatrixJ)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete diffusion operator for a scalar convection equation.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(IN)               :: ru

    ! scaling parameter
    real(DP), intent(IN)                           :: dscale

    ! perturbation parameter
    real(DP), intent(IN)                           :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(IN)                            :: bclear

     ! callback functions to compute velocity
    include 'intf_gfsccallback.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)                 :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(INOUT)            :: rmatrixJ   
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_MATIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    integer(PREC_VECIDX), dimension(:), pointer   :: p_IsuperdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdges
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_Kld,p_Ksep
    integer(PREC_VECIDX), dimension(:), pointer   :: p_Kcol
    real(DP), dimension(:,:), pointer             :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer               :: p_pp,p_pm,p_qp,p_qm,p_rp,p_rm
    real(DP), dimension(:), pointer               :: p_flux
    real(DP), dimension(:), pointer               :: p_u
    real(DP), dimension(:), pointer               :: p_Jac
    integer :: h_Ksep
    logical :: bisExtended

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrixJ)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case(AFCSTAB_SYMMETRIC)

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
      
      ! Check if subdiagonal edges need to be generated
      if (iand(rafcstab%iSpec, AFCSTAB_SUBDIAGONALEDGES) .eq. 0)&
          call afcstab_generateSubdiagEdges(rafcstab)
      
      ! Set pointers
      call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
      call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)
      call afcstab_getbase_IsubdiagEdge(rafcstab,    p_IsubdiagonalEdges)
      call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)
      call afcstab_getbase_DcoeffsAtEdge(rafcstab,   p_DcoefficientsAtEdge)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(1), p_pp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(2), p_pm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(3), p_qp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(4), p_qm)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(5), p_rp)
      call lsyssc_getbase_double(rafcstab%RnodalVectors(6), p_rm)
      call lsyssc_getbase_double(rafcstab%RedgeVectors(1),  p_flux)
      call lsyssc_getbase_Kcol(rmatrixJ,   p_Kcol)
      call lsyssc_getbase_double(rmatrixJ, p_Jac)
      call lsyssc_getbase_double(ru,       p_u)
      
      ! What kind of matrix format are we?
      select case(rmatrixJ%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixJ, p_Kld)
        
        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrixJ%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrixJ%NEQ+1)
        call lalg_vectorAddScalarInt(p_Ksep, 1)

        ! Assembled extended Jacobian matrix
        bisExtended = (rafcstab%iextendedJacobian .ne. 0)

        call doJacobian_Symm(p_IsuperdiagonalEdgesIdx, p_IverticesAtEdge,&
            p_IsubdiagonalEdgesIdx, p_IsubdiagonalEdges, p_DcoefficientsAtEdge,&
            p_Kld, p_Kcol, p_u, p_flux, p_pp, p_pm, p_qp, p_qm, p_rp, p_rm, dscale,&
            hstep, rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE, bisExtended, p_Ksep, p_Jac)

        ! Free storage
        call storage_free(h_Ksep)
        
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarSymm')
        call sys_halt()
      end select
      

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalarSymm')
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
    subroutine adjustKsep(Kld, Kcol, k, Ksep)
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      integer(PREC_VECIDX), intent(IN)                  :: k
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      
      integer(PREC_MATIDX) :: ild,isep
      integer(PREC_VECIDX) :: l
      integer :: iloc

      ! Loop over all entries of the k-th row
      do ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number and the position of the separator
        l = Kcol(ild); isep = Ksep(l)

        ! If the separator does not point to the k-th column
        ! it must be adjusted accordingly
        if (Kcol(isep) < k) Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsep


    !**************************************************************
    ! Assemble the Jacobian matrix for symmetric flux limiting
    subroutine doJacobian_Symm(IsuperdiagonalEdgesIdx, IverticesAtEdge,&
                               IsubdiagonalEdgesIdx, IsubdiagonalEdges,&
                               DcoefficientsAtEdge, Kld, Kcol, u, flux,&
                               pp, pm, qp, qm, rp, rm, dscale, hstep,&
                               NEQ, NEDGE, NNVEDGE, bisExtended, Ksep, Jac)
      
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: IsuperdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdgesIdx
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: IsubdiagonalEdges
      real(DP), dimension(:,:), intent(IN)              :: DcoefficientsAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kld
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      real(DP), dimension(:), intent(IN)                :: u
      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:), intent(IN)                :: pp,pm,qp,qm,rp,rm
      real(DP), intent(IN)                              :: dscale,hstep
      integer(PREC_VECIDX), intent(IN)                  :: NEQ
      integer(PREC_MATIDX), intent(IN)                  :: NEDGE
      integer, intent(IN)                               :: NNVEDGE
      logical, intent(IN)                               :: bisExtended
      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac

      ! local variables
      integer(PREC_VECIDX), dimension(5,NNVEDGE) :: Kloc
      real(DP), dimension(2,0:NNVEDGE) :: pploc,pmloc
      real(DP), dimension(2,0:NNVEDGE) :: qploc,qmloc
      real(DP), dimension(2,0:NNVEDGE) :: rploc,rmloc
      real(DP), dimension(2,0:NNVEDGE) :: fluxloc

      integer(PREC_MATIDX) :: ild,iedge
      integer(PREC_VECIDX) :: k,l
      integer :: iloc,nloc

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
          call updateJacobian_Symm(IverticesAtEdge, DcoefficientsAtEdge,&
              u, pp, pm, qp, qm, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagonalEdgesIdx(k), IsuperdiagonalEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          call updateJacobian_Symm(IverticesAtEdge, DcoefficientsAtEdge,&
              u, pp, pm, qp, qm, hstep, iedge, iloc, k,&
              pploc, pmloc, qploc, qmloc, fluxloc, Kloc)
        end do

        ! Save total number of local neighbors
        nloc = iloc
        
        ! Adjust the diagonal separator
        call adjustKsep(Kld, Kcol, k, Ksep)

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        rploc(:,0:nloc) = afcstab_limit( pploc(:,0:nloc), qploc(:,0:nloc), 0._DP, 1._DP)
        rmloc(:,0:nloc) = afcstab_limit(-pmloc(:,0:nloc),-qmloc(:,0:nloc), 0._DP, 1._DP)

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
            
            call assembleJacobian_Symm(IverticesAtEdge, Kld, Kcol, flux, rp, rm,&
                Kloc, rploc, rmloc, fluxloc, dscale, hstep,&
                iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
          
          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagonalEdgesIdx(l), IsuperdiagonalEdgesIdx(l+1)-1
            
            call assembleJacobian_Symm(IverticesAtEdge, Kld, Kcol, flux, rp, rm,&
                Kloc, rploc, rmloc, fluxloc, dscale, hstep,&
                iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do
      end do   ! end-of k-loop
    end subroutine doJacobian_Symm

    
    !**************************************************************
    ! Update the local coefficients for symmetric flux limiting
    subroutine updateJacobian_Symm(IverticesAtEdge, DcoefficientsAtEdge,&
                                   u, pp, pm, qp, qm, hstep, iedge, iloc, k,&
                                   pploc, pmloc, qploc, qmloc, fluxloc, Kloc)

      integer(PREC_MATIDX), dimension(:,:), intent(IN)     :: IverticesAtEdge
      real(DP), dimension(:,:), intent(IN)                 :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(IN)                   :: u
      real(DP), dimension(:), intent(IN)                   :: pp,pm,qp,qm
      real(DP), intent(IN)                                 :: hstep
      integer(PREC_MATIDX), intent(IN)                     :: iedge
      integer, intent(IN)                                  :: iloc
      integer(PREC_VECIDX), intent(IN)                     :: k

      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(INOUT)             :: pploc,pmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: qploc,qmloc
      real(DP), dimension(:,0:), intent(INOUT)             :: fluxloc
      integer(PREC_VECIDX), dimension(:,:), intent(INOUT)  :: Kloc

      ! local variables
      integer(PREC_VECIDX) :: i,j
      real(DP) :: d_ij,f_ij,s_ij,diff,hstep_ik,hstep_jk,dsign
      integer  :: iperturb


      ! Determine indices. Obviously, either i or j must be equal
      ! to k. Otherwise, the edge ij would not be present in the
      ! list of incident edges for node k.
      i = IverticesAtEdge(1,iedge)
      j = IverticesAtEdge(2,iedge)

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
        hstep_ik = hstep; hstep_jk = 0._DP
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        pploc(:,iloc) = pp(j)-max(0._DP, -f_ij)
        pmloc(:,iloc) = pm(j)-max(0._DP, -f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        qploc(:,iloc) = qp(j)-max(0._DP, -f_ij)
        qmloc(:,iloc) = qm(j)-min(0._DP, -f_ij)
        
      else
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = i
        
        ! Compute signed perturbation parameters
        hstep_ik = 0._DP; hstep_jk = hstep
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        pploc(:,iloc) = pp(i)-max(0._DP, f_ij)
        pmloc(:,iloc) = pm(i)-min(0._DP, f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        qploc(:,iloc) = qp(i)-max(0._DP, f_ij)
        qmloc(:,iloc) = qm(i)-min(0._DP, f_ij)
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
          pploc(iperturb,0)    = pploc(iperturb,0)+max(0._DP, f_ij)
          pmloc(iperturb,0)    = pmloc(iperturb,0)+min(0._DP, f_ij)
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, -f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))

          ! Update upper/lower bounds
          qploc(iperturb,0)    = qploc(iperturb,0)+max(0._DP, f_ij)
          qmloc(iperturb,0)    = qmloc(iperturb,0)+min(0._DP, f_ij)
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP, -f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP, -f_ij)
          
        else
          
          ! Compute raw antidiffusve flux
          f_ij = d_ij*(diff+dsign*(hstep_ik-hstep_jk))
          fluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          pploc(iperturb,iloc) = pploc(iperturb,iloc)+max(0._DP, f_ij)
          pmloc(iperturb,iloc) = pmloc(iperturb,iloc)+min(0._DP, f_ij)
          pploc(iperturb,0)    = pploc(iperturb,0)+max(0._DP, -f_ij)
          pmloc(iperturb,0)    = pmloc(iperturb,0)+min(0._DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))
          
          ! Update upper/lower bounds
          qploc(iperturb,iloc) = qploc(iperturb,iloc)+max(0._DP, f_ij)
          qmloc(iperturb,iloc) = qmloc(iperturb,iloc)+min(0._DP, f_ij)
          qploc(iperturb,0)    = qploc(iperturb,0)+max(0._DP, -f_ij)
          qmloc(iperturb,0)    = qmloc(iperturb,0)+min(0._DP, -f_ij)
        end if
      end do
    end subroutine updateJacobian_Symm

    
    !**************************************************************
    ! Assemble the given column of the Jacobian for symmetric flux limiting
    subroutine assembleJacobian_Symm(IverticesAtEdge, Kdiagonal, Kcol,&
                                     flux, rp, rm, Kloc, rploc, rmloc, fluxloc,&
                                     dscale, hstep, iedge, iloc, k, l,&
                                     bisExtended, Ksep, Jac)
      
      integer(PREC_MATIDX), dimension(:,:), intent(IN)  :: IverticesAtEdge
      integer(PREC_MATIDX), dimension(:), intent(IN)    :: Kdiagonal
      integer(PREC_VECIDX), dimension(:), intent(IN)    :: Kcol
      real(DP), dimension(:), intent(IN)                :: flux
      real(DP), dimension(:), intent(IN)                :: rp,rm
      integer(PREC_VECIDX), dimension(:,:), intent(IN)  :: Kloc
      real(DP), dimension(:,0:), intent(IN)             :: rploc,rmloc
      real(DP), dimension(:,0:), intent(IN)             :: fluxloc
      real(DP), intent(IN)                              :: dscale,hstep
      integer(PREC_MATIDX), intent(IN)                  :: iedge
      integer, intent(IN)                               :: iloc
      integer(PREC_VECIDX), intent(IN)                  :: k,l
      logical, intent(IN)                               :: bisExtended

      integer(PREC_MATIDX), dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT)             :: Jac
      
      ! local variables
      integer(PREC_MATIDX) :: ik,jk
      integer(PREC_VECIDX) :: i,j,m
      real(DP)             :: f_ij
      integer              :: iperturb

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
    end subroutine assembleJacobian_Symm
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
      bhasOrientation = .TRUE.
      
    case DEFAULT
      bhasOrientation = .FALSE.
    end select
  end function gfsc_hasOrientation

end module groupfemscalar

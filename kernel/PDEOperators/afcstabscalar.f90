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
!# 1.) afcsc_initStabilisation = afcsc_initStabByMatrix /
!#                               afcsc_initStabByGroupFEMSet
!#     -> Initialises the stabilisation structure
!#
!# 2.) afcsc_initPerfConfig
!#     -> Initialises the global performance configuration
!#
!# 3.) afcsc_renderOperatorLED
!#     -> Renders the given operator local extremum diminishing LED
!#
!# </purpose>
!##############################################################################

module afcstabscalar

!$use omp_lib
  use afcstabbase
  use fsystem
  use genoutput
  use groupfembase
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use spatialdiscretisation
  use storage

  implicit none

  private
  public :: afcsc_perfconfig
  public :: afcsc_initStabilisation
  public :: afcsc_initPerfConfig
  public :: afcsc_renderOperatorLED

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

  subroutine afcsc_initStabByMatrix(rmatrix, rafcstab,&
      rblockDiscretisation, rdiscretisation)

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

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary 1-block vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation
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

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    call afcstab_allocEdgeStructure(rafcstab,6)
    call afcstab_genEdgeList(rmatrix, rafcstab)

    ! Allocate internal data structures and auxiliary coefficients
    call afcstab_allocInternalData(rafcstab, .true.,&
        rblockDiscretisation, rdiscretisation)

  end subroutine afcsc_initStabByMatrix

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_initStabByGroupFEMSet(rgroupFEMSet, rafcstab,&
      rblockDiscretisation, rdiscretisation)

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

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary 1-block vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation
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

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
      rafcstab%h_IedgeListIdx     = rgroupFEMSet%h_IedgeListIdx
      rafcstab%h_IedgeList        = rgroupFEMSet%h_IedgeList
      rafcstab%iduplicationFlag   = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
      rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
                                        AFCSTAB_HAS_EDGELIST)
    else
      call output_line('Group finite element set does not provide edge structure',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
      call sys_halt()
    end if

    ! Allocate internal data structures and auxiliary coefficients
    call afcstab_allocInternalData(rafcstab, .true.,&
        rblockDiscretisation, rdiscretisation)

  end subroutine afcsc_initStabByGroupFEMSet

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_renderOperatorLED(rafcstab, rmatrix, rperfconfig)

!<description>
    ! This subroutine renders the given operator rmatrix local
    ! extremum diminishing (LED). This is achieved by adding a
    ! symmetric operator which has zero row and column sums and
    ! eliminates all negative off-diagonal entries.
!</description>

!<input>
    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    real(SP), dimension(:,:), pointer :: p_Fcoefficients

    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if

    ! Check if stabilisation structure is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST) .eq. 0) then
      call output_line('Stabilisation structure does not provide required data!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_renderOperatorLED')
      call sys_halt()
    end if

    ! Set pointers
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)

    ! What data types are we?
    select case(rmatrix%cdataType)
    case (ST_DOUBLE)
      ! Set pointers
      call lsyssc_getbase_double(rmatrix, p_Ddata)

      ! Check if coefficients should be stored in stabilisation
      if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then

        ! Check if stabilisation has the same data type
        if (rafcstab%cdataType .ne. ST_DOUBLE) then
          call output_line('Stabilisation must have double precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_renderOperatorLED')
          call sys_halt()
        end if

        ! Set additional pointers
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)

        !-----------------------------------------------------------------------
        ! Assemble operator with stabilisation and generate coefficients
        !-----------------------------------------------------------------------
        call doOperatorEdgeAFCDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%bisSymmetricOperator, p_Ddata, p_Dcoefficients)

        ! Set state of stabilisation
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)

        ! Do we need edge orientation?
        if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
          call afcstab_upwindOrientation(p_Dcoefficients, p_IedgeList, 2, 3)
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
        else
          rafcstab%istabilisationSpec =&
              iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
        end if
      else

        !-----------------------------------------------------------------------
        ! Assemble operator with stabilisation but do not generate coeffs
        !-----------------------------------------------------------------------
        call doOperatorEdgeDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%bisSymmetricOperator, p_Ddata)
      end if

    case (ST_SINGLE)
      ! Set pointers
      call lsyssc_getbase_single(rmatrix, p_Fdata)

      ! Check if coefficients should be stored in stabilisation
      if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then

        ! Check if stabilisation has the same data type
        if (rafcstab%cdataType .ne. ST_SINGLE) then
          call output_line('Stabilisation must have single precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_renderOperatorLED')
          call sys_halt()
        end if

        ! Set additional pointers
        call afcstab_getbase_FcoeffsAtEdge(rafcstab, p_Fcoefficients)

        !-----------------------------------------------------------------------
        ! Assemble operator with stabilisation and generate coefficients
        !-----------------------------------------------------------------------
        call doOperatorEdgeAFCSP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%bisSymmetricOperator, p_Fdata, p_Fcoefficients)

        ! Set state of stabilisation
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)       

        ! Do we need edge orientation?
        if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
          call afcstab_upwindOrientation(p_Fcoefficients, p_IedgeList, 2, 3)
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
        else
          rafcstab%istabilisationSpec =&
              iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
        end if
      else

        !-----------------------------------------------------------------------
        ! Assemble operator with stabilisation but do not generate coeffs
        !-----------------------------------------------------------------------
        call doOperatorEdgeSP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%bisSymmetricOperator, p_Fdata)
      end if

    end select

  contains

    ! Here, the real working routines start

    !**************************************************************
    ! Render the matrix local extremum without storing the artificial
    ! diffusion coefficient for later use

    subroutine doOperatorEdgeDP(IedgeListIdx, IedgeList, bsymm, Ddata)

      ! input parameters
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      logical, intent(in) :: bsymm

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      real(DP) :: d_ij
      integer :: iedge,igroup,ii,jj,ij,ji

      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      if (bsymm) then

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Symmetric artificial diffusion coefficient
            d_ij = max(0.0_DP, -Ddata(ij))

            ! Update the global operator
            Ddata(ii) = Ddata(ii) - d_ij
            Ddata(jj) = Ddata(jj) - d_ij
            Ddata(ij) = Ddata(ij) + d_ij
            Ddata(ji) = Ddata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      else

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Non-symmetric artificial diffusion coefficient
            d_ij = max(0.0_DP, -Ddata(ij), -Ddata(ji))

            ! Update the global operator
            Ddata(ii) = Ddata(ii) - d_ij
            Ddata(jj) = Ddata(jj) - d_ij
            Ddata(ij) = Ddata(ij) + d_ij
            Ddata(ji) = Ddata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      end if

      !$omp end parallel

    end subroutine doOperatorEdgeDP

    !**************************************************************
    ! Render the matrix local extremum without storing the artificial
    ! diffusion coefficient for later use

    subroutine doOperatorEdgeAFCDP(IedgeListIdx, IedgeList, bsymm,&
        Ddata, Dcoefficients)

      ! input parameters
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      logical, intent(in) :: bsymm

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! output parameters
      real(DP), dimension(:,:), intent(out) :: Dcoefficients

      ! local variables
      real(DP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,jj,ij,ji

      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      if (bsymm) then

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Symmetric artificial diffusion coefficient
            d_ij = max(0.0_DP, -Ddata(ij))
            k_ij = max(0.0_DP,  Ddata(ij))

            ! Symmetric AFC w/o edge orientation
            Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)

            ! Update the global operator
            Ddata(ii) = Ddata(ii) - d_ij
            Ddata(jj) = Ddata(jj) - d_ij
            Ddata(ij) = Ddata(ij) + d_ij
            Ddata(ji) = Ddata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      else

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Non-symmetric artificial diffusion coefficient
            d_ij = max(0.0_DP, -Ddata(ij), -Ddata(ji))
            k_ij = Ddata(ij) + d_ij
            k_ji = Ddata(ji) + d_ij

            ! Non-symmetric AFC w/o edge orientation
            Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)

            ! Update the global operator
            Ddata(ii) = Ddata(ii) - d_ij
            Ddata(jj) = Ddata(jj) - d_ij
            Ddata(ij) = Ddata(ij) + d_ij
            Ddata(ji) = Ddata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      end if

      !$omp end parallel

    end subroutine doOperatorEdgeAFCDP

    !**************************************************************
    ! Render the matrix local extremum without storing the artificial
    ! diffusion coefficient for later use

    subroutine doOperatorEdgeSP(IedgeListIdx, IedgeList, bsymm, Fdata)

      ! input parameters
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      logical, intent(in) :: bsymm

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      real(SP) :: d_ij
      integer :: iedge,igroup,ii,jj,ij,ji

      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      if (bsymm) then

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Symmetric artificial diffusion coefficient
            d_ij = max(0.0_SP, -Fdata(ij))

            ! Update the global operator
            Fdata(ii) = Fdata(ii) - d_ij
            Fdata(jj) = Fdata(jj) - d_ij
            Fdata(ij) = Fdata(ij) + d_ij
            Fdata(ji) = Fdata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      else

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Non-symmetric artificial diffusion coefficient
            d_ij = max(0.0_SP, -Fdata(ij), -Fdata(ji))

            ! Update the global operator
            Fdata(ii) = Fdata(ii) - d_ij
            Fdata(jj) = Fdata(jj) - d_ij
            Fdata(ij) = Fdata(ij) + d_ij
            Fdata(ji) = Fdata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      end if

      !$omp end parallel

    end subroutine doOperatorEdgeSP

    !**************************************************************
    ! Render the matrix local extremum without storing the artificial
    ! diffusion coefficient for later use

    subroutine doOperatorEdgeAFCSP(IedgeListIdx, IedgeList, bsymm,&
        Fdata, Fcoefficients)

      ! input parameters
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      logical, intent(in) :: bsymm

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! output parameters
      real(SP), dimension(:,:), intent(out) :: Fcoefficients

      ! local variables
      real(SP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,jj,ij,ji

      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      if (bsymm) then

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Symmetric artificial diffusion coefficient
            d_ij = max(0.0_SP, -Fdata(ij))
            k_ij = max(0.0_SP,  Fdata(ij))

            ! Symmetric AFC w/o edge orientation
            Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)

            ! Update the global operator
            Fdata(ii) = Fdata(ii) - d_ij
            Fdata(jj) = Fdata(jj) - d_ij
            Fdata(ij) = Fdata(ij) + d_ij
            Fdata(ji) = Fdata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      else

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            ii = IedgeList(5,iedge)
            jj = IedgeList(6,iedge)

            ! Non-symmetric artificial diffusion coefficient
            d_ij = max(0.0_SP, -Fdata(ij), -Fdata(ji))
            k_ij = Fdata(ij) + d_ij
            k_ji = Fdata(ji) + d_ij

            ! Non-symmetric AFC w/o edge orientation
            Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)

            ! Update the global operator
            Fdata(ii) = Fdata(ii) - d_ij
            Fdata(jj) = Fdata(jj) - d_ij
            Fdata(ij) = Fdata(ij) + d_ij
            Fdata(ji) = Fdata(ji) + d_ij
          end do
          !$omp end do
        end do ! igroup

      end if

      !$omp end parallel

    end subroutine doOperatorEdgeAFCSP

  end subroutine afcsc_renderOperatorLED

end module afcstabscalar

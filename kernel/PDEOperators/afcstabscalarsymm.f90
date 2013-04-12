!##############################################################################
!# ****************************************************************************
!# <name> afcstabscalarsymm </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the algebraic
!# flux correction methodology proposed by Kuzmin, Moeller and Turek
!# in a series of publications for symmetric diffusion operators. As a
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
!# afcstabbase.
!#
!# The following routines are available:
!#
!#
!# 1.) afcsc_buildVectorSymm = afcsc_buildVectorSymmScalar /
!#                             afcsc_buildVectorSymmBlock
!#     -> Assembles the vector for AFC stabilisation of symmetric type
!#
!# 2.) afcsc_buildJacobianSymm = afcsc_buildJacobianSymmScalar /
!#                               afcsc_buildJacobianSymmBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of
!#         symmetric flux limiting.
!#
!# </purpose>
!##############################################################################

module afcstabscalarsymm

!$use omp_lib
  use afcstabbase
  use afcstabscalar
  use basicgeometry
  use collection
  use fsystem
  use genoutput
  use groupfembase
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use storage

  implicit none

  private

  public :: afcsc_buildVectorSymm
  public :: afcsc_buildJacobianSymm

  !*****************************************************************************

  interface afcsc_buildVectorSymm
    module procedure afcsc_buildVectorSymmScalar
    module procedure afcsc_buildVectorSymmBlock
  end interface

  interface afcsc_buildJacobianSymm
    module procedure afcsc_buildJacobianSymmScalar
    module procedure afcsc_buildJacobianSymmBlock
  end interface

contains

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
    call doLimitDP(p_IedgeListIdx, p_IedgeList,&
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

    subroutine doLimitDP(IedgeListIdx, IedgeList,&
        DcoefficientsAtEdge, Dx, dscale, NEDGE,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(:), intent(inout) :: Dflux,Dy

      ! local variables
      real(DP) :: d_ij,f_ij,s_ij,diff
      integer :: i,iedge,igroup,j


      ! Clear nodal vectors
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      !$omp parallel default(shared)&
      !$omp private(d_ij,diff,f_ij,i,j,s_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
      ! in parallel without the need to synchronise memory access
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

    end subroutine doLimitDP

  end subroutine afcsc_buildVectorSymmScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianSymmBlock(rx, dscale, hstep, bclear,&
      rafcstab, rjacobian, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete diffusion operator for a
    ! scalar convection equation. Note that this routine serves as a
    ! wrapper for block vectors. If there is only one block, then the
    ! corresponding scalar routine is called. Otherwise, an error is
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

      ! input parameters
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

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Drploc,Drmloc,Dfluxloc
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ild,iedge,i,j,k,l,iloc,nloc


      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ

        ! Assemble nodal coefficients P and Q for node k and all vertices
        ! surrounding node k. Note that it suffices to initialise only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0

        ! Initialise local counter
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

      ! input parameters
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

end module afcstabscalarsymm

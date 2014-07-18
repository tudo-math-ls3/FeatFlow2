!##############################################################################
!# ****************************************************************************
!# <name> afcstabscalargp </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the
!# general-purpose variant of the algebraic flux correction
!# methodology proposed by Kuzmin, Moeller and Turek in a series of
!# publications. As a starting point for scalar conservation laws, the
!# reader is referred to the book chapter
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
!# 1.) afcsc_buildVectorGP = afcsc_buildVectorGPScalar /
!#                           afcsc_buildVectorGPBlock
!#     -> Assembles the vector for AFC stabilisation of general-purpose type
!#
!# 2.) afcsc_buildJacobianGP = afcsc_buildJacLinearGPScalar /
!#                             afcsc_buildJacLinearGPBlock /
!#                             afcsc_buildJacobianGPScalar /
!#                             afcsc_buildJacobianGPBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of
!#         general purpose limiter. For the first two routines, the
!#         velocity is assumed to be linear which simplifies the
!#         evaluation of the Jacobian matrix significantly. For the
!#         second two routines, the velocity can be arbitrary.
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) minmodQP / minmodDP / minmodSP
!#     -> Computes minmod(a,b)
!#
!# </purpose>
!##############################################################################

module afcstabscalargp

!$ use omp_lib
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

  public :: afcsc_buildVectorGP
  public :: afcsc_buildJacobianGP
  
  interface afcsc_buildVectorGP
    module procedure afcsc_buildVectorGPScalar
    module procedure afcsc_buildVectorGPBlock
  end interface

  interface afcsc_buildJacobianGP
    module procedure afcsc_buildJacLinearGPScalar
    module procedure afcsc_buildJacLinearGPBlock
    module procedure afcsc_buildJacobianGPScalar
    module procedure afcsc_buildJacobianGPBlock
  end interface

contains

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorGPBlock(rafcstab, rmatrix,&
      rx, rx0, theta, dscale, bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-GP type. Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
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
    !     J. Comput. Phys. 219  (2006) 513-531.
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
    call doLimitDP(p_IedgeListIdx, p_IedgeList,&
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

    subroutine doLimitDP(IedgeListIdx, IedgeList,&
        DcoefficientsAtEdge, MC, Dx, Dx0, theta, dscale, NEDGE,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dflux0, Dy)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0
      real(DP), intent(in) :: theta,dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(:), intent(inout) :: Dflux,Dflux0,Dy

      ! local variables
      real(DP) :: d_ij,df_ij,f_ij,l_ij,l_ji,m_ij,p_ij,pf_ij,q_ij,q_ji
      real(DP) :: diff,diff0,diff1
      integer :: i,iedge,igroup,ij,j


      ! Clear nodal vectors
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)

      !$omp parallel default(shared)&
      !$omp private(d_ij,df_ij,diff,diff0,diff1,&
      !$omp         f_ij,i,ij,j,l_ij,l_ji,m_ij,p_ij,pf_ij,q_ij,q_ji)&
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

    end subroutine doLimitDP

  end subroutine afcsc_buildVectorGPScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearGPBlock(rmatrix, rx,&
      rx0, theta, tstep, hstep, bclear, rafcstab, rjacobian,&
      bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation. Note that the velocity is assumed
    ! to be linear. Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the
    ! corresponding scalar routine is called. Otherwise, an error is
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

    case default
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

      ! input/ parameters
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

      ! input/output parameters
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

      ! input parameters
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

  subroutine afcsc_buildJacobianGPBlock(rgroupFEMSet, rmatrix,&
      rx, rx0, fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, bextendedSparsity, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation. The velocity is assumed to be
    ! nonlinear/arbitrary. Note that this routine serves as a wrapper
    ! for block vectors. If there is only one block, then the
    ! corresponding scalar routine is called. Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
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

      call afcsc_buildJacobianGPScalar(rgroupFEMSet,&
          rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), fcb_calcMatrixSc_sim, theta, tstep, hstep,&
          bclear, rafcstab, rjacobian,bextendedSparsity, rcollection, rperfconfig)

    end if
  end subroutine afcsc_buildJacobianGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianGPScalar(rgroupFEMSet, rmatrix,&
      rx, rx0, fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, bextendedSparsity, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary.
    ! This routine will also work for linear velocities but then it is inefficient
    ! since the solution perturbation does not affect the velocity.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
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

!!$    ! How many dimensions do we have?
!!$    ndim = size(RcoeffMatrices,1)
!!$    select case(ndim)
!!$    case (NDIM1D)
!!$      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
!!$
!!$    case (NDIM2D)
!!$      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
!!$      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
!!$
!!$    case (NDIM3D)
!!$      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
!!$      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
!!$      call lsyssc_getbase_double(RcoeffMatrices(3), p_DcoeffZ)
!!$
!!$    case default
!!$      call output_line('Unsupported spatial dimension!',&
!!$          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianGPScalar')
!!$      call sys_halt()
!!$    end select

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

    case default
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

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      ! input/output parameters
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

      ! input parameters
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

      ! input/output parameters
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

      ! input parameters
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

      ! input/output parameters
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

      ! input parameters
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

      ! input/output parameters
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

      ! input parameters
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

!<function>

  elemental function minmodQP(qa,qb)

!<description>
    ! This function computes minmod(a,b) following the implementation
    ! suggested by N. Robidoux, M. Gong, J. Cupitt, A. Turcott, and
    ! K. Martinez in "CPU, SMP and GPU Implementations of Nohalo Level
    ! 1, a Fast Co-Convex Antialiasing Image Resampler", ACM 2009
!</description>

!<input>
    real(QP), intent(in) :: qa,qb
!</input>

!<result>
    real(QP) :: minmodQP
!</result>

!</function>

    ! local variables
    real(QP) :: qsa,qsb

    qsa = sign(1._QP,qa)
    qsb = sign(1._QP,qb)

    minmodQP = 0.5_QP * (qsa+qsb) * min(qa*qsa, qb*qsb)

  end function minmodQP

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

!</function>

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

!</function>

    ! local variables
    real(SP) :: fsa,fsb

    fsa = sign(1._SP,fa)
    fsb = sign(1._SP,fb)

    minmodSP = 0.5_SP * (fsa+fsb) * min(fa*fsa, fb*fsb)

  end function minmodSP

end module afcstabscalargp

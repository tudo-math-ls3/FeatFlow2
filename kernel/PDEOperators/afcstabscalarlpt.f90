!##############################################################################
!# ****************************************************************************
!# <name> afcstabscalarlpt </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the
!# linearity-preserving variant of the algebraic flux correction
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
!# 1.) afcsc_buildVectorLPT = afcsc_buildVecLPTScalar /
!#                            afcsc_buildVecLPTBlock
!#      -> Assembles the vector for stabilisation by means of
!#         linearity-preserving flux correction
!#
!# 2.) afcsc_buildFluxLPT = afcsc_buildFluxLPTScalar /
!#                          afcsc_buildFluxLPTBlock
!#     -> Assembles the raw antidiffusive flux for the
!#        linearity-preserving stabilisation
!#
!# The following auxiliary routines are available:
!#
!# 1.) minmodQP / minmodDP / minmodSP
!#     -> Computes minmod(a,b)
!#
!# </purpose>
!##############################################################################

module afcstabscalarlpt

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

  public :: afcsc_buildVectorLPT
  public :: afcsc_buildFluxLPT

  !*****************************************************************************

  interface afcsc_buildVectorLPT
    module procedure afcsc_buildVecLPTScalar
    module procedure afcsc_buildVecLPTBlock
  end interface

  interface afcsc_buildFluxLPT
    module procedure afcsc_buildFluxLPTScalar
    module procedure afcsc_buildFluxLPTBlock
  end interface

contains

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVecLPTBlock(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rmatrix, rperfconfig)

!<description>
    ! This subroutine assembles the vector resulting from the
    ! application of linearyity-preserving flux correction. Note that
    ! this routine serves as a wrapper for block vectors. If there is
    ! only one block, then the corresponding scalar routine is
    ! called. Otherwise, an error is thrown.
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

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: lumped mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

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
    if (rx%nblocks .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTBlock')
      call sys_halt()

    else

      call afcsc_buildVecLPTScalar(rafcstab, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1),&
          rmatrix, rperfconfig)

    end if

  end subroutine afcsc_buildVecLPTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVecLPTScalar(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rmatrix, rperfconfig)

!<description>
    ! This subroutine assembles the vector resulting from the
    ! application of linearyity-preserving flux correction.
    !
    ! This subroutine provides different algorithms:
    !
    ! Nonlinear/Linearised FEM-LPT
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Linearity-perserving flux correction for upwind-biased and
    ! symmetric antidiffusive fluxes.
    !
    ! The details of this method can be found in:
    !
    ! D. Kuzmin, Linearity-preserving flux correction and convergence
    ! acceleration for constrained Galerkin schemes, Ergebnisberichte
    ! Angew. Math. 430, Technische Universitaet Dortmund, 2011.
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

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: lumped mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

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
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy,p_Dq
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if

    ! Check if stabilisation is prepeared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
    ! The FEM-LPT algorithm is split into the following steps which
    ! can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Compute the antidiffusive increments (Pp, Pm)
    !
    ! 3) Compute the local solution bounds (Qp, Qm).
    !
    ! 4) Compute the nodal correction factors (Rp, Rm).
    !
    ! 5) Compute edgewise correction factors (Alpha).
    !
    ! 6) Apply the limited antidifusive fluxes
    !-------------------------------------------------------------------------

    if (iand(ioperationSpec, AFCSTAB_LPTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 2) Compute sums of antidiffusive increments
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Compute sums of antidiffusive increments based on the
      ! raw-antidiffusive fluxes ...
      select case(rafcstab%climitingType)
      case (AFCSTAB_LIMITING_SYMMETRIC)
        ! ... in a symmetric fashion
        call doADIncrementsSymmDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)

      case (AFCSTAB_LIMITING_UPWINDBIASED)
        ! ... in an upwind-biased fashion
        call doADIncrementsUpwDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)

      case default
        call output_line('Unsupported type of flux limiting!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end select

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 3) Compute local solution bounds
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Compute bounds
      call doBoundsDP(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEDGE, p_Dx, p_Dqp, p_Dqm)

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 4) Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Set additional pointer
      call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)

      ! Compute nodal correction factors
      call doLimitNodalConstrainedDP(rafcstab%NEQ,&
          p_Dq, p_Dx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 5) Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide nodal correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Compute edgewise correction factors ...
      select case(rafcstab%climitingType)
      case (AFCSTAB_LIMITING_SYMMETRIC)
        ! ... in a symmetric fashion
        call doLimitEdgewiseSymmDP(p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dpp, p_Dpm, p_Drp, p_Drm, p_Dalpha)

      case (AFCSTAB_LIMITING_UPWINDBIASED)
        ! ... in an upwind-biased fashion
        call doLimitEdgewiseUpwDP(p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dpp, p_Dpm, p_Drp, p_Drm, p_Dalpha)

      case default
        call output_line('Unsupported type of flux limiting!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end select

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 6) Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide edgewise correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_LPTALGO_SCALEBYMASS) .ne. 0) then
        if (present(rmatrix)) then
          call lsyssc_getbase_double(rmatrix, p_ML)

          call doCorrectScaleByMassDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
        else
          call output_line('Lumped mass matrix is not provided!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
          call sys_halt()
        end if
      else

        call doCorrectDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Dy)
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes in a symmetric fashion

    subroutine doADIncrementsSymmDP(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! The sums of positive/negative antidiffusive increments
      real(DP), dimension(:), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP) :: f_ij,fp_ij,fm_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij,fp_ij,fm_ij)

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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)

          ! Separate fluxes into positive/negative contributions
          fp_ij = max(0.0_DP,f_ij)
          fm_ij = min(0.0_DP,f_ij)

          ! Compute the sums of antidiffusive increments
          Dpp(i) = Dpp(i) + fp_ij   ! += max(0.0_DP, f_ij)
          Dpp(j) = Dpp(j) - fm_ij   ! += max(0.0_DP,-f_ij)
          Dpm(i) = Dpm(i) + fm_ij   ! += min(0.0_DP, f_ij)
          Dpm(j) = Dpm(j) - fp_ij   ! += min(0.0_DP,-f_ij)
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsSymmDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes in an upwind-biased fashion

    subroutine doADIncrementsUpwDP(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup

      !$omp parallel default(shared) private(i,f_ij)

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

          ! Get node number of the upwind nodeee
          i = IedgeList(1,iedge)

          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)

          ! Compute the sums of antidiffusive increments
          Dpp(i) = Dpp(i) + max(0.0_DP,f_ij)
          Dpm(i) = Dpm(i) + min(0.0_DP,f_ij)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsUpwDP

    !**************************************************************
    ! Assemble the local bounds from the predicted solution

    subroutine doBoundsDP(IedgeListIdx, IedgeList, NEDGE, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dqp,Dqm

      ! local variables
      integer :: i,iedge,igroup,j


      !$omp parallel default(shared) private(i,j)

      ! Initialise Q`s by solution
      !$omp sections
      !$omp section
      call lalg_copyVector(Dx, Dqp)
      !$omp section
      call lalg_copyVector(Dx, Dqm)
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Compute local upper and lower bounds
          Dqp(i) = max(Dqp(i), Dx(j))
          Dqm(i) = min(Dqm(i), Dx(j))
          Dqp(j) = max(Dqp(j), Dx(i))
          Dqm(j) = min(Dqm(j), Dx(i))
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDP

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrainedDP(NEQ, Dq, Dx,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dq,Dx
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      integer, intent(in) :: NEQ

      ! output parameters
      real(DP), dimension(:), intent(out) :: Drp,Drm

      ! local variables
      integer :: ieq

      !$omp parallel sections default(shared) private(ieq)

      !$omp section

      !$omp parallel do default(shared) private(ieq)
      do ieq = 1, NEQ
        Drp(ieq) = min(1.0_DP, (Dqp(ieq)-Dx(ieq)+AFCSTAB_EPSABS) /&
                               (Dpp(ieq)+AFCSTAB_EPSABS) * Dq(ieq))
      end do
      !$omp end parallel do

      !$omp section

      !$omp parallel do default(shared) private(ieq)
      do ieq = 1, NEQ
        Drm(ieq) = min(1.0_DP, (Dqm(ieq)-Dx(ieq)-AFCSTAB_EPSABS) /&
                               (Dpm(ieq)-AFCSTAB_EPSABS) * Dq(ieq))
      end do
      !$omp end parallel do

      !$omp end parallel sections

    end subroutine doLimitNodalConstrainedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes in
    ! a symmetric fashion, e.g., consider the nodal correction factors
    ! at both endpoints of the edge

    subroutine doLimitEdgewiseSymmDP(IedgeList,&
        NEDGE, Dflux, Dpp, Dpm, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dpp,Dpm
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factors from previous
      !           multiplicative correction steps
      ! On exit: the edge-wise correction factors resulting from
      !          the nodal correction factors Rp and Rm
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP) :: f_ij,r_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,f_ij,r_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)

        ! Compute nodal correction factors
        if (f_ij .gt. AFCSTAB_EPSABS) then
          r_ij = min(Drp(i),Drm(j))
        elseif (f_ij .lt. -AFCSTAB_EPSABS) then
          r_ij = min(Drp(j),Drm(i))
        else
          r_ij = 1.0_DP
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseSymmDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes in
    ! an upwind-biased fashion, e.g., consider the nodal correction
    ! factors only at the endpoint of the edge located upwind

    subroutine doLimitEdgewiseUpwDP(IedgeList,&
        NEDGE, Dflux, Dpp, Dpm, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dpp,Dpm
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factors from previous
      !           multiplicative correction steps
      ! On exit: the edge-wise correction factors resulting from
      !          the nodal correction factors Rp and Rm
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP) :: f_ij,r_ij
      integer :: iedge,i

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,f_ij,r_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node number of the upwind node
        i = IedgeList(1,iedge)

        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)

        ! Compute nodal correction factors
        if (f_ij .gt. AFCSTAB_EPSABS) then
          r_ij = Drp(i)
        elseif (f_ij .lt. -AFCSTAB_EPSABS) then
          r_ij = Drm(i)
        else
          r_ij = 1.0_DP
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseUpwDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectDP(IedgeListIdx, IedgeList,&
        NEDGE, dscale, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dy

      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij)&
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)

          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij
          Dy(j) = Dy(j) - f_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectScaleByMassDP(IedgeListIdx, IedgeList,&
        NEDGE, dscale, ML, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: ML,Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dy

      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij)&
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)

          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij/ML(i)
          Dy(j) = Dy(j) - f_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectScaleByMassDP

  end subroutine afcsc_buildVecLPTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildFluxLPTBlock(rafcstab, rx, dscale, bclear,&
      ioperationSpec, rmatrix, rperfconfig)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! linearity-preserving flux correction, whereby the coefficients
    ! are determined from the off-diagonal entries of the matrix. Note
    ! that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: coefficient matrix which is used instead of the
    ! coefficients at edges provided by the stabilisation structure
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    if (rx%nblocks .ne. 1) then
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTBlock')
      call sys_halt()
    end if

    ! Call subroutine for scalar vectors
    call afcsc_buildFluxLPTScalar(rafcstab, rx%RvectorBlock(1),&
        dscale, bclear, ioperationSpec, rmatrix, rperfconfig)

  end subroutine afcsc_buildFluxLPTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildFluxLPTScalar(rafcstab, rx, dscale, bclear,&
      ioperationSpec, rmatrix, rperfconfig)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! linearity-preserving flux correction, whereby the coefficients
    ! are determined from the off-diagonal entries of the matrix.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: coefficient matrix which is used instead of the
    ! coefficients at edges provided by the stabilisation structure
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx
    real(DP), dimension(:), pointer :: p_Dflux,p_Dq
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:,:), pointer :: p_IedgeList

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
        call sys_halt()
      end if
    end if


    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case(AFCSTAB_NLINLPT_MASS,&
         AFCSTAB_LINLPT_MASS)

      if (present(rmatrix)) then

        ! Set pointer
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX) .ne. 0) then

          ! Set pointers
          call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Assemble raw antidiffusive fluxes for mass antidiffusion
          call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
              p_Dmatrix, p_Dx, dscale, bclear, p_Dflux)

          ! Set specifiers for raw antidiffusive fluxes
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        end if

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX_BOUNDS) .ne. 0) then

          ! Check if stabilisation provides edgewise bounds
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS) .ne. 0) then

            ! Set pointers
            call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
            call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
            call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
            call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)

            ! Assemble bounds based on matrix coefficients
            call doBoundsByMatrixDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dmatrix, p_DboundsAtEdge, p_Dq)

            ! Set specifiers for nodal bounds
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)

          else

            call output_line('Stabilisation does not provide bounds at edges!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
        call sys_halt()
      end if

      !-------------------------------------------------------------------------

    case(AFCSTAB_NLINLPT_UPWINDBIASED)

      ! Check if stabilisation provides edge data
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .ne. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .ne. 0)) then

        ! Set pointer
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX) .ne. 0) then

          ! Set pointers
          call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Assemble upwind-biased raw antidiffusive fluxes
          call doFluxesUpwindBiasedDP(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale, bclear, p_Dflux)

          ! Set specifiers for raw antidiffusive fluxes
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        end if

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX_BOUNDS) .ne. 0) then

          ! Check if stabilisation provides edgewise bounds
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS) .ne. 0) then

            ! Set pointers
            call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
            call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
            call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
            call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)

            ! Assemble bounds based on edge coefficients
            call doBoundsByCoeffDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_DcoefficientsAtEdge, p_DboundsAtEdge, p_Dq)

            ! Set specifiers for nodal bounds
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)

          else

            call output_line('Stabilisation does not provide bounds at edges!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Stabilisation does not provide (oriented) edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
        call sys_halt()
      end if

      !-------------------------------------------------------------------------

    case(AFCSTAB_NLINLPT_SYMMETRIC,&
         AFCSTAB_LINLPT_SYMMETRIC,&
         AFCSTAB_LINLPT_UPWINDBIASED)

      ! Remark: in the linearised version, all fluxes must be limited
      ! in a symmetric fashion since the flux into the downwind node
      ! is no longer balanced by the nonoscillatory part of the
      ! Galerkin operator; cf. comment by D. Kuzmin, 2011

      ! Check if stabilisation provides edge data
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .ne. 0) then

        ! Set pointer
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX) .ne. 0) then

          ! Set pointers
          call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Assemble upwind-biased raw antidiffusive fluxes
          call doFluxesSymmetricDP(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale, bclear, p_Dflux)

          ! Set specifiers for raw antidiffusive fluxes
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        end if

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX_BOUNDS) .ne. 0) then

          ! Check if stabilisation provides edgewise bounds
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS) .ne. 0) then

            ! Set pointers
            call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
            call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
            call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
            call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)

            ! Assemble bounds based on edge coefficients
            call doBoundsByCoeffDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_DcoefficientsAtEdge, p_DboundsAtEdge, p_Dq)

            ! Set specifiers for nodal bounds
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)

          else

            call output_line('Stabilisation does not provide bounds at edges!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Stabilisation does not provide edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
        call sys_halt()
      end if

    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble symmetric raw antidiffusive fluxes using the
    ! coefficients supplied by the CSR-matrix Dmatrix

    subroutine doFluxesByMatrixDP(IedgeList, NEDGE,&
        Dmatrix, Dx, dscale, bclear, Dflux)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dmatrix, Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,i,j,ij

      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Store raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(i)-Dx(j))
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

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(i)-Dx(j))
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

            ! Store raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(j)-Dx(i))
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

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(j)-Dx(i))
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

            ! Store raw antidiffusive flux
            Dflux(iedge) = dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
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

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByMatrixDP

    !**************************************************************
    ! Assemble upwind-biased raw antidiffusive fluxes using the
    ! coefficients supplied by the edge-by-edge array DcoefficientsAtEdge

    subroutine doFluxesUpwindBiasedDP(IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux

      ! local variables
      real(DP) :: f_ij,d_ij,l_ji,g_ij
      integer :: iedge,i,j

      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))

            ! Store raw antidiffusive flux f_ij := minmod(f_ij,g_ij)
            Dflux(iedge) = minmodDP(f_ij, g_ij)
          end do
          !$omp end parallel do

        else
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))

            ! Update raw antidiffusive flux by f_ij := minmod(f_ij,g_ij)
            Dflux(iedge) = Dflux(iedge) + minmodDP(f_ij, g_ij)
          end do
          !$omp end parallel do
        end if

      elseif (dscale .eq. -1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))

            ! Store raw antidiffusive flux f_ij := minmod(f_ij,g_ij)
            Dflux(iedge) = -minmodDP(f_ij, g_ij)
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))

            ! Update raw antidiffusive flux by f_ij := minmod(f_ij,g_ij)
            Dflux(iedge) = Dflux(iedge) - minmodDP(f_ij, g_ij)
          end do
          !$omp end parallel do
        end if

      else

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))

            ! Store raw antidiffusive flux f_ij := minmod(f_ij,g_ij)
            Dflux(iedge) = dscale*minmodDP(f_ij, g_ij)
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))

            ! Update raw antidiffusive flux by f_ij := minmod(f_ij,g_ij)
            Dflux(iedge) = Dflux(iedge) + dscale*minmodDP(f_ij, g_ij)
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesUpwindBiasedDP

    !**************************************************************
    ! Assemble symmetric raw antidiffusive fluxes using the
    ! coefficients supplied by the edge-by-edge array DcoefficientsAtEdge

    subroutine doFluxesSymmetricDP(IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux

            ! local variables
      integer :: iedge,i,j

      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Store raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
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

            ! Store raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
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

            ! Store raw antidiffusive flux
            Dflux(iedge) = dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesSymmetricDP

    !**************************************************************
    ! Assemble nodal bounds using the coefficients supplied by the
    ! CSR-matrix Dmatrix and the precomputed bounds at edges

    subroutine doBoundsByMatrixDP(IedgeListIdx, IedgeList,&
        NEDGE, Dmatrix, DboundsAtEdge, Dq)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DboundsAtEdge
      real(DP), dimension(:), intent(in) :: Dmatrix
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dq

      ! local variables
      integer :: i,iedge,igroup,ij,j,ji

      ! Clear bounds`s
      call lalg_clearVector(Dq)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      !$omp parallel default(shared) private(iedge,i,j,ij,ji)
      do igroup = 1, size(IedgeListIdx)-1

        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers and matrix positions
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Update bounds
          Dq(i) = Dq(i) + DboundsAtEdge(1,iedge)*Dmatrix(ij)
          Dq(j) = Dq(j) + DboundsAtEdge(2,iedge)*Dmatrix(ji)
        end do
        !$omp end do
      end do
      !$omp end parallel

    end subroutine doBoundsByMatrixDP

    !**************************************************************
    ! Assemble nodal bounds using the coefficients supplied by the
    !  edge-by-edge array DcoefficientsAtEdge

    subroutine doBoundsByCoeffDP(IedgeListIdx, IedgeList,&
        NEDGE, DcoefficientsAtEdge, DboundsAtEdge, Dq)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:,:), intent(in) :: DboundsAtEdge
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dq

      ! local variables
      integer :: i,iedge,igroup,j

      ! Clear bounds`s
      call lalg_clearVector(Dq)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      !$omp parallel default(shared) private(iedge,i,j)
      do igroup = 1, size(IedgeListIdx)-1

        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Update bounds
          Dq(i) = Dq(i) + DboundsAtEdge(1,iedge)*DcoefficientsAtEdge(1,iedge)
          Dq(j) = Dq(j) + DboundsAtEdge(2,iedge)*DcoefficientsAtEdge(1,iedge)
        end do
        !$omp end do
      end do
      !$omp end parallel

    end subroutine doBoundsByCoeffDP

  end subroutine afcsc_buildFluxLPTScalar

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

end module afcstabscalarlpt

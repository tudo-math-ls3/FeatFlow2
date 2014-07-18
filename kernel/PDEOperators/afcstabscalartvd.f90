!##############################################################################
!# ****************************************************************************
!# <name> afcstabscalartvd </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the TVD-type
!# algebraic flux correction methodology proposed by Kuzmin, Moeller
!# and Turek in a series of publications. As a starting point for
!# scalar conservation laws, the reader is referred to the book chapter
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
!# 1.) afcsc_buildVectorTVD = afcsc_buildVectorTVDScalar /
!#                            afcsc_buildVectorTVDBlock
!#     -> Assembles the vector for AFC stabilisation of TVD type
!#
!# 2.) afcsc_buildJacobianTVD = afcsc_buildJacLinearTVDScalar /
!#                              afcsc_buildJacLinearTVDBlock /
!#                              afcsc_buildJacobianTVDScalar /
!#                              afcsc_buildJacobianTVDBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of TVD
!#         type; For the first two routines, the velocity is assumed
!#         to be linear which simplifies the evaluation of the
!#         Jacobian matrix significantly. For the second two
!#         routines, the velocity can be arbitrary.
!#
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) minmodQP / minmodDP / minmodSP
!#     -> Computes minmod(a,b)
!#
!# </purpose>
!##############################################################################

module afcstabscalartvd

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

  public :: afcsc_buildVectorTVD
  public :: afcsc_buildJacobianTVD

  interface afcsc_buildVectorTVD
    module procedure afcsc_buildVectorTVDScalar
    module procedure afcsc_buildVectorTVDBlock
  end interface

  interface afcsc_buildJacobianTVD
    module procedure afcsc_buildJacLinearTVDScalar
    module procedure afcsc_buildJacLinearTVDBlock
    module procedure afcsc_buildJacobianTVDScalar
    module procedure afcsc_buildJacobianTVDBlock
  end interface

contains

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorTVDBlock(rafcstab, rx, dscale, bclear,&
      ioperationSpec, ry, fcb_calcADFluxes, fcb_calcADIncrements,&
      fcb_calcBounds, fcb_limitNodal, fcb_limitEdgewise,&
      fcb_calcCorrection, rcollection, rperfconfig)


!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-TVD type. Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
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
    ! combination of different AFCSTAB_TVDALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: callback functions to overwrite the standard operations
    include 'intf_calcADFluxes.inc'
    optional :: fcb_calcADFluxes

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks   .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDBlock')
      call sys_halt()

    else

      call afcsc_buildVectorTVDScalar(rafcstab, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1),&
          fcb_calcADFluxes, fcb_calcADIncrements, fcb_calcBounds,&
          fcb_limitNodal, fcb_limitEdgewise, fcb_calcCorrection,&
          rcollection, rperfconfig)

    end if

  end subroutine afcsc_buildVectorTVDBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorTVDScalar(rafcstab, rx, dscale, bclear,&
      ioperationSpec, ry, fcb_calcADFluxes, fcb_calcADIncrements,&
      fcb_calcBounds, fcb_limitNodal, fcb_limitEdgewise,&
      fcb_calcCorrection, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-TVD type.
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

    ! OPTIONAL: callback functions to overwrite the standard operations
    include 'intf_calcADFluxes.inc'
    optional :: fcb_calcADFluxes

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm
    real(DP), dimension(:), pointer :: p_Dqp,p_Dqm
    real(DP), dimension(:), pointer :: p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dx,p_Dy,p_Dalpha,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    logical :: bquickPrepare,bquickCorrect

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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
      call sys_halt()
    end if

    ! Clear destination vector?
    if (bclear) call lsyssc_clearVector(ry)

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
    call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

    !---------------------------------------------------------------------------
    ! The FEM-FCT algorithm is split into the following steps which
    ! can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Compute the raw antidiffusive fluxes (Flux).
    !
    ! 3) Compute the antidiffusive increments (Pp, Pm).
    !
    ! 4) Compute the local solution bounds (Qp, Qm).
    !
    ! 5) Compute the nodal correction factors (Rp, Rm).
    !
    ! 6) Compute edgewise correction factors (Alpha).
    !
    ! 7) Apply the limited antiddifusive fluxes.
    !
    ! Note that steps 2)-4) can be performed simultaneously within a single
    ! loop over the edges to increase performance. This also applies to 6)-7)
    !-------------------------------------------------------------------------

    if (iand(ioperationSpec, AFCSTAB_TVDALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if


    if ((iand(ioperationSpec, AFCSTAB_TVDALGO_ADFLUXES)     .ne. 0) .and.&
        (iand(ioperationSpec, AFCSTAB_TVDALGO_ADINCREMENTS) .ne. 0) .and.&
        (iand(ioperationSpec, AFCSTAB_TVDALGO_BOUNDS)       .ne. 0)) then
      !-------------------------------------------------------------------------
      ! 2)-4) Compute raw antidiffusive fluxes, sums of antidiffusive
      !       increments, and local solution bounds simultaneously
      !-------------------------------------------------------------------------

      ! Check if stabilisation is prepared
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0)) then
        call output_line('Stabilisation does not provide required structures!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
        call sys_halt()
      end if

      ! Check that no user-defined callback routinee are given
      if (.not.present(fcb_calcADFluxes)     .and.&
          .not.present(fcb_calcADIncrements) .and.&
          .not.present(fcb_calcBounds)) then

        ! Use quick prepare routine
        call doQuickPrepareDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_DcoefficientsAtEdge, p_Dx, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Dflux)

        bquickPrepare = .true.

        ! Set specifier
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
      else
        bquickPrepare = .false.
      end if
    end if

    ! Perform 2)-4) step-by-step?
    if (.not.bquickPrepare) then
      if (iand(ioperationSpec, AFCSTAB_TVDALGO_ADFLUXES) .ne. 0) then
        !-----------------------------------------------------------------------
        ! 2) Compute raw antidiffusive fluxes
        !-----------------------------------------------------------------------

        ! Check if stabilisation is prepared
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0)) then
          call output_line('Stabilisation does not provide required structures!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        if (present(fcb_calcADFluxes)) then
          ! User-supplied callback routine
          call fcb_calcADFluxes(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ,&
              dscale, p_Dx, p_Dflux, rcollection)
        else
          ! Standard routine
          call doADFluxesDP(p_IedgeListIdx, p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale, p_Dflux)
        end if

        ! Set specifier
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
      end if


      if (iand(ioperationSpec, AFCSTAB_TVDALGO_ADINCREMENTS) .ne. 0) then
        !-----------------------------------------------------------------------
        ! 3) Compute sums of antidiffusive increments
        !-----------------------------------------------------------------------

        ! Check if stabilisation is prepared
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide required structures!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        if (present(fcb_calcADIncrements)) then
          ! User-supplied callback routine
          call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ, &
              p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm, rcollection=rcollection)
        else
          ! Standard routine
          call doADIncrementsDP(p_IedgeListIdx, p_IedgeList, rafcstab%NEDGE,&
              p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        end if

        ! Set specifier
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
      end if


      if (iand(ioperationSpec, AFCSTAB_TVDALGO_BOUNDS) .ne. 0) then
        !-----------------------------------------------------------------------
        ! 4) Compute local solution bounds
        !-----------------------------------------------------------------------

        ! Check if stabilisation is prepared
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide required structures!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        if (present(fcb_calcBounds)) then
          ! User-supplied callback routine
          call fcb_calcBounds(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ,&
              p_Dx, p_Dqp, p_Dqm, rcollection=rcollection)
        else
          ! Standard routine
          call doBoundsDP(p_IedgeListIdx, p_IedgeList, rafcstab%NEDGE,&
              p_Dflux, p_Dalpha, p_Dqp, p_Dqm)
        end if

        ! Set specifier
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
      end if
    end if


    if (iand(ioperationSpec, AFCSTAB_TVDALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 5) Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
        call sys_halt()
      end if

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, 1, 1.0_DP,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection=rcollection)
      else
        ! Apply the standard nodal limiter
        !$omp parallel sections
        !$omp section
        p_Drp = afcstab_limit(p_Dpp, p_Dqp, 0.0_DP, 1.0_DP)

        !$omp section
        p_Drm = afcstab_limit(p_Dpm, p_Dqm, 0.0_DP, 1.0_DP)
        !$omp end parallel sections
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


     if ((iand(ioperationSpec, AFCSTAB_TVDALGO_LIMITEDGE) .ne. 0) .and.&
         (iand(ioperationSpec, AFCSTAB_TVDALGO_CORRECT)   .ne. 0)) then
       !------------------------------------------------------------------------
       ! 6)-7) Compute edgewise correction factors and apply correction
       !------------------------------------------------------------------------

       ! Check if stabilisation is prepared
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide required structures!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        ! Check that no user-defined callback rouintes are given
        if (.not.(present(fcb_limitEdgewise)) .and.&
            .not.(present(fcb_calcCorrection))) then

          ! Check if stabilisation provides raw antidiffusive fluxes
          ! and nodal correction factors
          if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)    .eq. 0) .or.&
              (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0)) then
            call output_line('Stabilisation does not provide antidiffusive fluxes '//&
                'and/or nodal correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
            call sys_halt()
          end if

          ! Use quick correction routin
          call doQuickCorrectDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, p_Dflux, p_Drp, p_Drm, p_Dalpha, p_Dy)

          bquickCorrect = .true.

          ! Set specifier
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
        else
          bquickCorrect = .false.
        end if
      end if


      ! Perform 6)-7) step-by-step
      if (.not.bquickCorrect) then
        if (iand(ioperationSpec, AFCSTAB_TVDALGO_LIMITEDGE) .ne. 0) then
          !---------------------------------------------------------------------
          ! 6) Compute edgewise correction factors
          !---------------------------------------------------------------------

          ! Check if stabilisation is prepared
          if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
              (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
            call output_line('Stabilisation does not provide required structures!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
            call sys_halt()
          end if

          ! Check if stabilisation provides raw antidiffusive fluxes
          ! and nodal correction factors
          if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)    .eq. 0) .or.&
              (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0)) then
            call output_line('Stabilisation does not provide antidiffusive fluxes '//&
                'and/or nodal correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
            call sys_halt()
          end if

          if (present(fcb_limitEdgewise)) then
            ! User-supplied callback routine
            call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ, &
                p_Dx, p_Dflux, p_Dalpha, p_Drp, p_Drm, rcollection=rcollection)
          else
            ! Standard routine
            call doLimitEdgewiseDP(p_IedgeList,&
                rafcstab%NEDGE, p_Dflux, p_Drp, p_Drm, p_Dalpha)
          end if

          ! Set specifier
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
        end if


        if (iand(ioperationSpec, AFCSTAB_TVDALGO_CORRECT) .ne. 0) then
          !---------------------------------------------------------------------
          ! 7) Apply limited antidiffusive fluxes
          !---------------------------------------------------------------------

          ! Check if stabilisation is prepared
          if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
              (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
            call output_line('Stabilisation does not provide required structures!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
            call sys_halt()
          end if

          ! Check if stabilisation provides raw antidiffusive fluxes
          ! and edgewise correction factors
          if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)    .eq. 0) .or.&
              (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0)) then
            call output_line('Stabilisation does not provide antidiffusive fluxes '//&
                'and/or edgewise correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
            call sys_halt()
          end if

          if (present(fcb_calcCorrection)) then
            ! User-supplied callback routine
            call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, 1, 1, rafcstab%NEQ, dscale,&
                p_Dx, p_Dalpha, p_Dflux, p_Dy, rcollection=rcollection)
          else
            ! Standard routine
            call doCorrectDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dalpha, p_Dflux, p_Dy)
          end if
        end if
      end if

  contains

    ! Here, the working routine follows

    !**************************************************************
    ! Assemble the raw antidiffusive fluxes, the sums of antidiffusive
    ! increments and the local bounds simultaneously

    subroutine doQuickPrepareDP(IedgeListIdx, IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, Dpp, Dpm, Dqp, Dqm, Dflux)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dpp,Dpm,Dqp,Dqm,Dflux

      ! local variables
      real(DP) :: d_ij,diff,f_ij,fm_ij,fp_ij,l_ji
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(d_ij,diff,f_ij,fm_ij,fp_ij,i,iedge,j,l_ji)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear P`s and Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

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
          l_ji = DcoefficientsAtEdge(3,iedge)

          ! Determine solution difference
          diff = Dx(i)-Dx(j)

          ! Prelimit the antidiffusive flux
          ! F`_IJ=MIN(-P_IJ,L_JI)(DX_I-DX_J)
          f_ij = dscale*min(d_ij,l_ji)*diff

          ! And store it
          Dflux(iedge) = f_ij

          ! Separate fluxes into positive/negative contributions
          fp_ij = max(0.0_DP,f_ij)
          fm_ij = min(0.0_DP,f_ij)

          ! Assemble P`s accordingly
          Dpp(i) = Dpp(i) + fp_ij   ! += max(0.0_DP, f_ij)
          Dpm(i) = Dpm(i) + fm_ij   ! += min(0.0_DP, f_ij)

          ! Assemble Q`s
          Dqp(i) = Dqp(i) - fm_ij   ! += max(0.0_DP,-f_ij)
          Dqp(j) = Dqp(j) + fp_ij   ! += max(0.0_DP, f_ij)
          Dqm(i) = Dqm(i) - fp_ij   ! += min(0.0_DP,-f_ij)
          Dqm(j) = Dqm(j) + fm_ij   ! += min(0.0_DP, f_ij)
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doQuickPrepareDP

    !**************************************************************
    ! Assemble the raw antidiffusive fluxes

    subroutine doADFluxesDP(IedgeListIdx, IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, Dflux)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dflux

      ! local variables
      real(DP) :: d_ij,diff,f_ij,l_ji
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(d_ij,diff,f_ij,i,iedge,j,l_ji)&
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
          l_ji = DcoefficientsAtEdge(3,iedge)

          ! Determine solution difference
          diff = Dx(i)-Dx(j)

          ! Prelimit the antidiffusive flux
          ! F`_IJ=MIN(-P_IJ,L_JI)(DX_I-DX_J)
          f_ij = dscale*min(d_ij,l_ji)*diff

          ! And store it
          Dflux(iedge) = f_ij
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doADFluxesDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes

    subroutine doADIncrementsDP(IedgeListIdx, IedgeList, NEDGE,&
        Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux,Dalpha
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup

      !$omp parallel default(shared) private(i,iedge,f_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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

        ! Loop over the edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Determine index
          i = IedgeList(1,iedge)

          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)

          ! Assemble P`s accordingly
          Dpp(i) = Dpp(i) + max(0.0_DP, f_ij)
          Dpm(i) = Dpm(i) + min(0.0_DP, f_ij)
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes

    subroutine doBoundsDP(IedgeListIdx, IedgeList, NEDGE,&
        Dflux, Dalpha, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux,Dalpha
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dqp,Dqm

      ! local variables
      real(DP) :: f_ij,fm_ij,fp_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(f_ij,fm_ij,fp_ij,i,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

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

          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)

          ! Separate fluxes into positive/negative contributions
          fp_ij = max(0.0_DP,f_ij)
          fm_ij = min(0.0_DP,f_ij)

          ! Assemble Q`s
          Dqp(i) = Dqp(i) - fm_ij   ! += max(0.0_DP,-f_ij)
          Dqp(j) = Dqp(j) + fp_ij   ! += max(0.0_DP, f_ij)
          Dqm(i) = Dqm(i) - fp_ij   ! += min(0.0_DP,-f_ij)
          Dqm(j) = Dqm(j) + fm_ij   ! += min(0.0_DP, f_ij)
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDP

    !**************************************************************
    ! Perform edgewise limiting and apply limited antidiffusive fluxes

    subroutine doQuickCorrectDP(IedgeListIdx, IedgeList, NEDGE,&
        Dflux, Drp, Drm, Dalpha, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux,Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dy

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dalpha

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

        ! Loop over the edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Get precomputed raw antidiffusive flux
          f_ij = Dflux(iedge)

          ! Determine upwind correction factor
          if (f_ij .gt. 0.0_DP) then
            Dalpha(iedge) = Drp(i)
          else
            Dalpha(iedge) = Drm(i)
          end if

          ! Limit raw antidiffusive fluxe
          f_ij = Dalpha(iedge) * f_ij

          ! Update the vector
          Dy(i) = Dy(i)+f_ij
          Dy(j) = Dy(j)-f_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doQuickCorrectDP

    !**************************************************************
    ! Perform edgewise limiting

    subroutine doLimitEdgewiseDP(IedgeList, NEDGE,&
        Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dflux,Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE

      ! output parameters
      real(DP), dimension(:), intent(out) :: Dalpha

      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,f_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Determine index
        i = IedgeList(1,iedge)

        ! Get precomputed raw antidiffusive flux
        f_ij = Dflux(iedge)

        ! Determine upwind correction factor
        if (f_ij .gt. 0.0_DP) then
          Dalpha(iedge) = Drp(i)
        else
          Dalpha(iedge) = Drm(i)
        end if
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectDP(IedgeListIdx, IedgeList, NEDGE,&
        Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha,Dflux
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
          f_ij = Dalpha(iedge) * Dflux(iedge)

          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij
          Dy(j) = Dy(j) - f_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectDP

  end subroutine afcsc_buildVectorTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearTVDBlock(rx, tstep, hstep,&
      bclear, rafcstab, rjacobian, bextendedSparsity)

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
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearTVDBlock')
      call sys_halt()

    else

      call afcsc_buildJacLinearTVDScalar(rx%RvectorBlock(1), tstep,&
          hstep, bclear, rafcstab, rjacobian, bextendedSparsity)

    end if
  end subroutine afcsc_buildJacLinearTVDBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearTVDScalar(rx, tstep, hstep,&
      bclear, rafcstab, rjacobian, bextendedSparsity)

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
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm
    real(DP), dimension(:), pointer :: p_Jac,p_Dx,p_Dflux
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearTVDScalar')
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

      ! Create diagonal separator and increase it by one
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)

      call doJacobianMat79_TVD(&
          p_IsuperdiagEdgesIdx, p_IedgeList,&
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
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)

      call doJacobianMat79_TVD(&
          p_IsuperdiagEdgesIdx, p_IedgeList,&
          p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
          p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
          p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
          tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
          rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)

      ! Free storage
      call storage_free(h_Ksep)

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearTVDScalar')
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
    ! Assemble the Jacobian matrix for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
        Dx, Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
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
          call updateJacobianMat79_TVD(&
              IedgeList, DcoefficientsAtEdge, Dx,&
              Dpp, Dpm, Dqp, Dqm, tstep, hstep, iedge, iloc, k,&
              Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1

          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              IedgeList, DcoefficientsAtEdge, Dx,&
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
                IedgeList, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1

            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux,&
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
    subroutine updateJacobianMat79_TVD(IedgeList,&
        DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        iedge, iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
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
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)

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
    subroutine assembleJacobianMat79_TVD(IedgeList, Kdiagonal,&
        Dflux, Kloc, Drploc, Drmloc, Dfluxloc, hstep, iedge, iloc, k, l,&
        bisExtended, Ksep, Jac)

      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IedgeList
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
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)
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

  end subroutine afcsc_buildJacLinearTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianTVDBlock(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, tstep, hstep, bclear, rafcstab,&
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

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDBlock')
      call sys_halt()

    else

      call afcsc_buildJacobianTVDScalar(rgroupFEMSet,&
          rx%RvectorBlock(1), fcb_calcMatrixSc_sim,&
          tstep, hstep, bclear, rafcstab, rjacobian,&
          bextendedSparsity, rcollection, rperfconfig)

    end if
  end subroutine afcsc_buildJacobianTVDBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianTVDScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, tstep, hstep, bclear, rafcstab,&
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
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Dflux
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Jac,p_Dx
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDScalar')
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
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)

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
!!$          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDScalar')
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
        call doJacobianMat79_TVD_1D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dpp, p_Dpm,&
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
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian,   p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)

      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_TVD_1D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dpp, p_Dpm,&
            p_Dqp, p_Dqm, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      end select

      ! Free storage
      call storage_free(h_Ksep)

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDScalar')
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
    ! Assemble the Jacobian matrix for FEM-TVD in 1D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD_1D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
        DcoeffX, Dx, Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
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
      real(DP), dimension(NDIM1D) :: c_ij, c_ji
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = DcoeffX(ij)
          c_ji = DcoeffX(ji)

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
                IedgeList, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1

            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux,&
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
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, DcoeffY, Dx, Dflux,&
        Dpp, Dpm, Dqp, Dqm, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
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
      real(DP), dimension(NDIM2D) :: c_ij, c_ji
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji)/)

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
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1

            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
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
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, Dx,&
        Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
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
      real(DP), dimension(NDIM3D) :: c_ij, c_ji
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij),DcoeffZ(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji),DcoeffZ(ji)/)

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
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1

            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
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

      ! input parameters
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

!!$        ! Compute perturbed coefficients k_ij and k_ji
!!$        call fcb_calcMatrix(Dx(i)+dsign*hstep_ik, Dx(j)+dsign*hstep_jk,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
    subroutine assembleJacobianMat79_TVD(IedgeList, Kdiagonal,&
        Dflux, Kloc, Drploc, Drmloc, Dfluxloc, hstep, iedge, iloc, k, l,&
        bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IedgeList,Kloc
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

  end subroutine afcsc_buildJacobianTVDScalar
  
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

end module afcstabscalartvd

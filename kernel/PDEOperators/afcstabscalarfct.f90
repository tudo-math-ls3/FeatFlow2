!#############################################################################
!# ****************************************************************************
!# <name> afcstabscalarfct </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the FCT-type
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
!# 1.) afcsc_buildVectorFCT = afcsc_buildVectorFCTScalar /
!#                            afcsc_buildVectorFCTBlock /
!#     -> Assembles the vector for AFC stabilisation of FCT-type
!#
!# 2.) afcsc_buildFluxFCT = afcsc_buildFluxFCTScalar /
!#                          afcsc_buildFluxFCTBlock
!#     -> Assembles the raw antidiffusive flux for AFC stabilisation of FCT-type
!#
!# 3.) afcsc_buildJacobianFCT = afcsc_buildJacLinearFCTScalar /
!#                              afcsc_buildJacLinearFCTBlock /
!#                              afcsc_buildJacobianFCTScalar /
!#                              afcsc_buildJacobianFCTBlock
!#     -> Assembles the Jacobian matrix for the stabilisation part of FCT
!#        type; For the first two routines, the velocity is assumed
!#        to be linear which simplifies the evaluation of the
!#        Jacobian matrix significantly. For the second two
!#        routines, the velocity can be arbitrary.
!#
!# </purpose>
!##############################################################################

module afcstabscalarfct

#include "kernel/feat2macros.h"

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

  public :: afcsc_buildVectorFCT
  public :: afcsc_buildFluxFCT
  public :: afcsc_buildJacobianFCT

  !*****************************************************************************

  interface afcsc_buildVectorFCT
    module procedure afcsc_buildVectorFCTScalar
    module procedure afcsc_buildVectorFCTBlock
  end interface

    interface afcsc_buildFluxFCT
    module procedure afcsc_buildFluxFCTScalar
    module procedure afcsc_buildFluxFCTBlock
  end interface

  interface afcsc_buildJacobianFCT
    module procedure afcsc_buildJacLinearFCTScalar
    module procedure afcsc_buildJacLinearFCTBlock
    module procedure afcsc_buildJacobianFCTScalar
    module procedure afcsc_buildJacobianFCTBlock
  end interface

contains

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorFCTBlock(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry,&
      fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
      fcb_limitEdgewise, fcb_calcCorrection, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FCT-type. Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: callback functions to overwrite the standard operations
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
    if (rx%nblocks .ne. 1 .or. ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTBlock')
      call sys_halt()

    else

      call afcsc_buildVectorFCTScalar(rafcstab, rmatrix, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1),&
          fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
          fcb_limitEdgewise, fcb_calcCorrection, rcollection, rperfconfig)

    end if

  end subroutine afcsc_buildVectorFCTBlock

  ! *****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorFCTScalar(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry,&
      fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
      fcb_limitEdgewise, fcb_calcCorrection, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FCT-type. The idea of flux corrected transport can be
    ! traced back to the early SHASTA algorithm by Boris and Bock in
    ! the early 1970s. Zalesak suggested a fully multi-dimensional
    ! generalisation of this approach and paved the way for a large
    ! family of FCT algorithms.
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
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: callback functions to overwrite the standard operations
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
    real(QP), dimension(:), pointer :: p_Qml,p_Qx,p_Qy
    real(DP), dimension(:), pointer :: p_Dml,p_Dx,p_Dy
    real(SP), dimension(:), pointer :: p_Fml,p_Fx,p_Fy
    real(QP), dimension(:), pointer :: p_Qpp,p_Qpm,p_Qqp,p_Qqm,p_Qrp,p_Qrm
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(SP), dimension(:), pointer :: p_Fpp,p_Fpm,p_Fqp,p_Fqm,p_Frp,p_Frm
    real(QP), dimension(:), pointer :: p_Qalpha,p_Qflux,p_QfluxPrel
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_DfluxPrel
    real(SP), dimension(:), pointer :: p_Falpha,p_Fflux,p_FfluxPrel
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
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set integer pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)

    ! Set data pointers according to precision
    select case(rafcstab%cdataType)
    case (ST_QUAD)
       call lsyssc_getbase_quad(rafcstab%p_rvectorPp, p_Qpp)
       call lsyssc_getbase_quad(rafcstab%p_rvectorPm, p_Qpm)
       call lsyssc_getbase_quad(rafcstab%p_rvectorQp, p_Qqp)
       call lsyssc_getbase_quad(rafcstab%p_rvectorQm, p_Qqm)
       call lsyssc_getbase_quad(rafcstab%p_rvectorRp, p_Qrp)
       call lsyssc_getbase_quad(rafcstab%p_rvectorRm, p_Qrm)
       call lsyssc_getbase_quad(rafcstab%p_rvectorAlpha, p_Qalpha)
       call lsyssc_getbase_quad(rafcstab%p_rvectorFlux, p_Qflux)

    case (ST_DOUBLE)
       call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
       call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
       call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
       call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
       call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
       call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
       call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
       call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)

    case (ST_SINGLE)
       call lsyssc_getbase_single(rafcstab%p_rvectorPp, p_Fpp)
       call lsyssc_getbase_single(rafcstab%p_rvectorPm, p_Fpm)
       call lsyssc_getbase_single(rafcstab%p_rvectorQp, p_Fqp)
       call lsyssc_getbase_single(rafcstab%p_rvectorQm, p_Fqm)
       call lsyssc_getbase_single(rafcstab%p_rvectorRp, p_Frp)
       call lsyssc_getbase_single(rafcstab%p_rvectorRm, p_Frm)
       call lsyssc_getbase_single(rafcstab%p_rvectorAlpha, p_Falpha)
       call lsyssc_getbase_single(rafcstab%p_rvectorFlux, p_Fflux)
       
    case default
       call output_line('Unsupported data type',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
       call sys_halt()
    end select

    ! Set pointers according to precision
    select case(rx%cdataType)
    case (ST_QUAD)
       call lsyssc_getbase_quad(rx, p_Qx)
    case (ST_DOUBLE)
       call lsyssc_getbase_double(rx, p_Dx)
    case (ST_SINGLE)
       call lsyssc_getbase_single(rx, p_Fx)
    case default
       call output_line('Unsupported data type',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
       call sys_halt()
    end select
    
    select case(ry%cdataType)
    case (ST_QUAD)
       call lsyssc_getbase_quad(ry, p_Qy)
    case (ST_DOUBLE)
       call lsyssc_getbase_double(ry, p_Dy)
    case (ST_SINGLE)
       call lsyssc_getbase_single(ry, p_Fy)
    case default
       call output_line('Unsupported data type',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
       call sys_halt()
    end select

    !---------------------------------------------------------------------------
    ! The FEM-FCT algorithm is split into the following steps which
    ! can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Prelimit the antidiffusive fluxes (alpha).
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
    !-------------------------------------------------------------------------

    if (iand(ioperationSpec, AFCSTAB_FCTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------

      ! Initialise alpha by unity
      select case(rafcstab%cdataType)
      case(ST_QUAD)
         call lalg_setVector(p_Qalpha, 1.0_QP)
      case(ST_DOUBLE)
         call lalg_setVector(p_Dalpha, 1.0_DP)
      case(ST_SINGLE)
         call lalg_setVector(p_Falpha, 1.0_SP)
      end select
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_PRELIMIT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 2) Prelimit the raw antidiffusive fluxes (if required)
      !-------------------------------------------------------------------------
      if (rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) then

        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
          call sys_halt()
        end if


        if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
          ! Perform standard prelimiting
          select case(rafcstab%cdataType)
          case (ST_QUAD)
            call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
            call doStdPrelimitQP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Qflux, p_QfluxPrel, p_Qalpha)
          case (ST_DOUBLE)
            call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
            call doStdPrelimitDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dflux, p_DfluxPrel, p_Dalpha)
          case (ST_SINGLE)
            call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
            call doStdPrelimitSP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Fflux, p_FfluxPrel, p_Falpha)
          end select
        elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
          ! Perform minmod prelimiting
          select case(rafcstab%cdataType)
          case (ST_QUAD)
             call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
             call doMinModPrelimitQP(p_IedgeListIdx, p_IedgeList,&
                  rafcstab%NEDGE, p_Qflux, p_QfluxPrel, p_Qalpha)
          case (ST_DOUBLE)
             call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
             call doMinModPrelimitDP(p_IedgeListIdx, p_IedgeList,&
                  rafcstab%NEDGE, p_Dflux, p_DfluxPrel, p_Dalpha)
          case (ST_SINGLE)
             call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
             call doMinModPrelimitSP(p_IedgeListIdx, p_IedgeList,&
                  rafcstab%NEDGE, p_Fflux, p_FfluxPrel, p_Falpha)
          end select
        else
          call output_line('Invalid type of prelimiting!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
          call sys_halt()
        end if
      end if
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 3) Compute sums of antidiffusive increments
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Set additional pointer
          select case(rafcstab%cdataType)
          case (ST_QUAD)
             call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
          case (ST_DOUBLE)
             call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
          case (ST_SINGLE)
             call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
          end select

        ! Compute sums of antidiffusive increments
        ! based on the prelimiting fluxes
        if (present(fcb_calcADIncrements)) then
          ! User-defined callback routine
          call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ,&
              p_Dx, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm, rcollection=rcollection)
        else
          ! Standard routine
          select case(rafcstab%cdataType)
          case (ST_QUAD)
             call doADIncrementsQP(p_IedgeListIdx, p_IedgeList,&
                  rafcstab%NEDGE, p_QfluxPrel, p_Qalpha, p_Qpp, p_Qpm)
          case (ST_DOUBLE)
             call doADIncrementsDP(p_IedgeListIdx, p_IedgeList,&
                  rafcstab%NEDGE, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          case (ST_SINGLE)
             call doADIncrementsSP(p_IedgeListIdx, p_IedgeList,&
                  rafcstab%NEDGE, p_FfluxPrel, p_Falpha, p_Fpp, p_Fpm)
          end select
        end if

      else

        ! Compute sums of antidiffusive increments
        ! based on the raw-antidiffusive fluxes
        if (present(fcb_calcADIncrements)) then
           ! User-defined callback routine
           call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ, &
                p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm, rcollection=rcollection)
        else
           ! Standard routine
           select case(rafcstab%cdataType)
           case (ST_QUAD)
              call doADIncrementsQP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, p_Qflux, p_Qalpha, p_Qpp, p_Qpm)
           case (ST_DOUBLE)
              call doADIncrementsDP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
           case (ST_SINGLE)
              call doADIncrementsSP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, p_Fflux, p_Falpha, p_Fpp, p_Fpm)
          end select
        end if
      end if

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 4) Compute local solution bounds
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Compute bounds
      if (present(fcb_calcBounds)) then
        ! User-supplied callback routine
        call fcb_calcBounds(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ,&
            p_Dx, p_Dqp, p_Dqm, rcollection=rcollection)
      else
        ! Standard routine
        select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
        case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
           call doBoundsQPQP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Qx, p_Qqp, p_Qqm)
        case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
           call doBoundsDPQP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dx, p_Qqp, p_Qqm)
        case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
           call doBoundsSPQP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Fx, p_Qqp, p_Qqm)

        case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
           call doBoundsQPDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Qx, p_Dqp, p_Dqm)
        case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
           call doBoundsDPDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dx, p_Dqp, p_Dqm)
        case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
           call doBoundsSPDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Fx, p_Dqp, p_Dqm)

        case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
           call doBoundsQPSP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Qx, p_Fqp, p_Fqm)
        case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
           call doBoundsDPSP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dx, p_Fqp, p_Fqm)
        case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
           call doBoundsSPSP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Fx, p_Fqp, p_Fqm)
        end select
      end if

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 5) Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Set additional pointers
      select case(rx%cdataType)
      case(ST_QUAD)
         call lsyssc_getbase_quad(rmatrix, p_Qml)
      case(ST_DOUBLE)
         call lsyssc_getbase_double(rmatrix, p_Dml)
      case(ST_SINGLE)
         call lsyssc_getbase_single(rmatrix, p_Fml)
      end select

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, 1, dscale,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dml, rcollection)
      elseif (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
         ! Standard routine without constraints
         select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
         case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
            call doLimitNodalQPQP(rafcstab%NEQ, dscale,&
                 p_Qml, p_Qx, p_Qpp, p_Qpm, p_Qqp, p_Qqm, p_Qrp, p_Qrm)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
            call doLimitNodalDPQP(rafcstab%NEQ, dscale,&
                 p_Dml, p_Dx, p_Qpp, p_Qpm, p_Qqp, p_Qqm, p_Qrp, p_Qrm)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
            call doLimitNodalSPQP(rafcstab%NEQ, dscale,&
                 p_Fml, p_Fx, p_Qpp, p_Qpm, p_Qqp, p_Qqm, p_Qrp, p_Qrm)
            
         case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
            call doLimitNodalQPDP(rafcstab%NEQ, dscale,&
                 p_Qml, p_Qx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
            call doLimitNodalDPDP(rafcstab%NEQ, dscale,&
                 p_Dml, p_Dx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
            call doLimitNodalSPDP(rafcstab%NEQ, dscale,&
                 p_Fml, p_Fx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
            
         case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
            call doLimitNodalQPSP(rafcstab%NEQ, dscale,&
                 p_Qml, p_Qx, p_Fpp, p_Fpm, p_Fqp, p_Fqm, p_Frp, p_Frm)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
            call doLimitNodalDPSP(rafcstab%NEQ, dscale,&
                 p_Dml, p_Dx, p_Fpp, p_Fpm, p_Fqp, p_Fqm, p_Frp, p_Frm)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
            call doLimitNodalSPSP(rafcstab%NEQ, dscale,&
                 p_Fml, p_Fx, p_Fpp, p_Fpm, p_Fqp, p_Fqm, p_Frp, p_Frm)
         end select

      else
         ! Standard routine with constraints
         select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
         case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
            call doLimitNodalConstrainedQPQP(rafcstab%NEQ, dscale,&
                 p_Qml, p_Qx, p_Qpp, p_Qpm, p_Qqp, p_Qqm, p_Qrp, p_Qrm)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
            call doLimitNodalConstrainedDPQP(rafcstab%NEQ, dscale,&
                 p_Dml, p_Dx, p_Qpp, p_Qpm, p_Qqp, p_Qqm, p_Qrp, p_Qrm)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
            call doLimitNodalConstrainedSPQP(rafcstab%NEQ, dscale,&
                 p_Fml, p_Fx, p_Qpp, p_Qpm, p_Qqp, p_Qqm, p_Qrp, p_Qrm)
            
         case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
            call doLimitNodalConstrainedQPDP(rafcstab%NEQ, dscale,&
                 p_Qml, p_Qx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
            call doLimitNodalConstrainedDPDP(rafcstab%NEQ, dscale,&
                 p_Dml, p_Dx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
            call doLimitNodalConstrainedSPDP(rafcstab%NEQ, dscale,&
                 p_Fml, p_Fx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
            
         case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
            call doLimitNodalConstrainedQPSP(rafcstab%NEQ, dscale,&
                 p_Qml, p_Qx, p_Fpp, p_Fpm, p_Fqp, p_Fqm, p_Frp, p_Frm)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
            call doLimitNodalConstrainedDPSP(rafcstab%NEQ, dscale,&
                 p_Dml, p_Dx, p_Fpp, p_Fpm, p_Fqp, p_Fqm, p_Frp, p_Frm)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
            call doLimitNodalConstrainedSPSP(rafcstab%NEQ, dscale,&
                 p_Fml, p_Fx, p_Fpp, p_Fpm, p_Fqp, p_Fqm, p_Frp, p_Frm)
         end select
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 6) Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      ! and nodal correction factors
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)    .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0)) then
        call output_line('Stabilisation does not provide antidiffusive fluxes '//&
            'and/or nodal correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Compute edgewise correction factors
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
         select case(rafcstab%cdataType)
         case(ST_QUAD)
            call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
         case(ST_DOUBLE)
            call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
         case(ST_SINGLE)
            call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
         end select

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ, &
              p_Dx, p_Dflux, p_Dalpha, p_Drp, p_Drm, DfluxConstr=p_DfluxPrel,&
              rcollection=rcollection)
        else
          ! Standard routine
           select case(rafcstab%cdataType)
           case(ST_QUAD)
              call doLimitEdgewiseConstrainedQP(p_IedgeList,&
                   rafcstab%NEDGE, p_QfluxPrel, p_Qflux, p_Qrp, p_Qrm, p_Qalpha)
           case(ST_DOUBLE)
              call doLimitEdgewiseConstrainedDP(p_IedgeList,&
                   rafcstab%NEDGE, p_DfluxPrel, p_Dflux, p_Drp, p_Drm, p_Dalpha)
           case(ST_SINGLE)
              call doLimitEdgewiseConstrainedSP(p_IedgeList,&
                   rafcstab%NEDGE, p_FfluxPrel, p_Fflux, p_Frp, p_Frm, p_Falpha)
           end select
        end if

      else

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, 1, rafcstab%NEQ,&
              p_Dx, p_Dflux, p_Dalpha, p_Drp, p_Drm, rcollection=rcollection)
        else
          ! Standard routine
           select case(rafcstab%cdataType)
           case(ST_QUAD)
              call doLimitEdgewiseQP(p_IedgeList,&
                   rafcstab%NEDGE, p_Qflux, p_Qrp, p_Qrm, p_Qalpha)
           case(ST_DOUBLE)
              call doLimitEdgewiseDP(p_IedgeList,&
                   rafcstab%NEDGE, p_Dflux, p_Drp, p_Drm, p_Dalpha)
           case(ST_SINGLE)
              call doLimitEdgewiseSP(p_IedgeList,&
                   rafcstab%NEDGE, p_Fflux, p_Frp, p_Frm, p_Falpha)
           end select
        end if
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 7) Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide edgewise correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then

        ! Set pointer
         select case(rx%cdataType)
         case(ST_QUAD)
            call lsyssc_getbase_quad(rmatrix, p_Qml)
         case(ST_DOUBLE)
            call lsyssc_getbase_double(rmatrix, p_Dml)
         case(ST_SINGLE)
            call lsyssc_getbase_single(rmatrix, p_Fml)
         end select

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, rafcstab%NEQ, dscale,&
              p_Dx, p_Dalpha, p_Dflux, p_Dy, p_Dml, rcollection)
        else
          ! Standard routine
         select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
         case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
          call doCorrectScaleByMassQPQP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Qml, p_Qalpha, p_Qflux, p_Qy)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
          call doCorrectScaleByMassDPQP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Dml, p_Qalpha, p_Qflux, p_Dy)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
          call doCorrectScaleByMassSPQP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Fml, p_Qalpha, p_Qflux, p_Fy)
            
         case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
          call doCorrectScaleByMassQPDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Qml, p_Dalpha, p_Dflux, p_Qy)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
          call doCorrectScaleByMassDPDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Dml, p_Dalpha, p_Dflux, p_Dy)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
          call doCorrectScaleByMassSPDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Fml, p_Dalpha, p_Dflux, p_Fy)
            
         case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
          call doCorrectScaleByMassQPSP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Qml, p_Falpha, p_Fflux, p_Qy)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
          call doCorrectScaleByMassDPSP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Dml, p_Falpha, p_Fflux, p_Dy)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
          call doCorrectScaleByMassSPSP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_Fml, p_Falpha, p_Fflux, p_Fy)
         end select

        end if

      else

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, 1, 1, rafcstab%NEQ, dscale,&
              p_Dx, p_Dalpha, p_Dflux, p_Dy, rcollection=rcollection)
        else
          ! Standard routine
           select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
           case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
              call doCorrectQPQP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Qalpha, p_Qflux, p_Qy)
           case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
              call doCorrectDPQP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Qalpha, p_Qflux, p_Dy)
           case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
              call doCorrectSPQP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Qalpha, p_Qflux, p_Fy)

           case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
              call doCorrectQPDP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Qy)
           case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
              call doCorrectDPDP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Dy)
           case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
              call doCorrectSPDP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Fy)

           case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
              call doCorrectQPSP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Falpha, p_Fflux, p_Qy)
           case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
              call doCorrectDPSP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Falpha, p_Fflux, p_Dy)
           case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
              call doCorrectSPSP(p_IedgeListIdx, p_IedgeList,&
                   rafcstab%NEDGE, dscale, p_Falpha, p_Fflux, p_Fy)
           end select
        end if
      end if
    end if

  contains

    ! Here, the working routines follow


    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes the standard way, as
    ! suggested by Boris and Book in their first FCT algorithm

#define TemplateType_AFC QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doStdPrelimit.h"
#undef TemplateType_AFC

#define TemplateType_AFC DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doStdPrelimit.h"
#undef TemplateType_AFC

#define TemplateType_AFC SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doStdPrelimit.h"
#undef TemplateType_AFC



    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter

#define TemplateType_AFC QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doMinModPrelimit.h"
#undef TemplateType_AFC

#define TemplateType_AFC DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doMinModPrelimit.h"
#undef TemplateType_AFC

#define TemplateType_AFC SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doMinModPrelimit.h"
#undef TemplateType_AFC



    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without prelimiting

#define TemplateType_AFC QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doADIncrements.h"
#undef TemplateType_AFC

#define TemplateType_AFC DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doADIncrements.h"
#undef TemplateType_AFC

#define TemplateType_AFC SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doADIncrements.h"
#undef TemplateType_AFC



    !**************************************************************
    ! Assemble the local bounds from the predicted solution

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doBounds.h"
#undef TemplateType_Vector
#undef TemplateType_AFC



    !**************************************************************
    ! Compute the nodal correction factors without constraints

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodal.h"
#undef TemplateType_Vector
#undef TemplateType_AFC



    !**************************************************************
    ! Compute nodal correction factors with constraints

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitNodalConstrained.h"
#undef TemplateType_Vector
#undef TemplateType_AFC



    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

#define TemplateType_AFC    QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitEdgewise.h"
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitEdgewise.h"
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitEdgewise.h"
#undef TemplateType_AFC



    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

#define TemplateType_AFC    QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitEdgewiseConstrained.h"
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitEdgewiseConstrained.h"
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doLimitEdgewiseConstrained.h"
#undef TemplateType_AFC



    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrect.h"
#undef TemplateType_Vector
#undef TemplateType_AFC



    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildVectorFCTScalar_doCorrectScaleByMass.h"
#undef TemplateType_Vector
#undef TemplateType_AFC


  end subroutine afcsc_buildVectorFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildFluxFCTBlock(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, fcb_calcFluxFCTSc_sim,&
      rgroupFEMSet, rmatrix, rxTimeDeriv, rxPredictor, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! algebraic flux correction of FCT-type with or without the
    ! contributions of the consistent mass matrix. Note that this
    ! routine serves as a wrapper for block vectors. If there is only
    ! one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Switch for flux assembly
    ! TRUE  : fluxes are not modified externally so that
    !         quicker assembly procedures may be feasible
    ! FALSE : fluxes are truely assembled even if this
    !         leads to an expensive addition of zeros
    logical, intent(in) :: bquickAssembly

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: callback functions to compute antidiffusive fluxes
    include 'intf_calcFluxFCTSc_sim.inc'
    optional :: fcb_calcFluxFCTSc_sim

    ! OPTIONAL: group finite element set
    type(t_groupFEMSet), intent(in), optional :: rgroupFEMSet

    ! OPTIONAL: Consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorBlock), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorBlock), intent(in), optional :: rxPredictor

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    integer :: nblocks

    ! Check if block vector(s) contains exactly one block
    nblocks = rx%nblocks
    if (present(rxTimeDeriv)) nblocks = max(nblocks, rxTimeDeriv%nblocks)
    if (present(rxPredictor)) nblocks = max(nblocks, rxPredictor%nblocks)

    if (nblocks .ne. 1) then
      call output_line('Vector(s) must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Call subroutine for scalar vectors
    if (present(rxTimeDeriv)) then
      if (present(rxPredictor)) then
        ! ... both approximate time derivative and predictor are present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rxTimeDeriv%RvectorBlock(1), rxPredictor%RvectorBlock(1), &
            rcollection, rperfconfig)
      else
        ! ... only the approximate time derivative is present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rxTimeDeriv=rxTimeDeriv%RvectorBlock(1),&
            rcollection=rcollection, rperfconfig=rperfconfig)
      end if
    else
      if (present(rxPredictor)) then
        ! ... only the predictor is present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rxPredictor=rxPredictor%RvectorBlock(1),&
            rcollection=rcollection, rperfconfig=rperfconfig)
      else
        ! ... neither the approximate time derivative nor the predictor is present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rcollection=rcollection, rperfconfig=rperfconfig)
      end if
    end if

  end subroutine afcsc_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>


  subroutine afcsc_buildFluxFCTScalar(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, fcb_calcFluxFCTSc_sim,&
      rgroupFEMSet, rmatrix, rxTimeDeriv, rxPredictor, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! algebraic flux correction of FCT-type with or without the
    ! contribution of the consistent mass matrix.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Switch for flux assembly
    ! TRUE  : fluxes are not modified externally so that
    !         quicker assembly procedures may be feasible
    ! FALSE : fluxes are truely assembled even if this
    !         leads to an expensive addition of zeros
    logical, intent(in) :: bquickAssembly

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: callback functions to compute antidiffusive fluxes
    include 'intf_calcFluxFCTSc_sim.inc'
    optional :: fcb_calcFluxFCTSc_sim

    ! OPTIONAL: group finite element set
    type(t_groupFEMSet), intent(in), optional :: rgroupFEMSet

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorScalar), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorScalar), intent(in), optional :: rxPredictor

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(QP), dimension(:), pointer :: p_Qmatrix,p_Qx
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx
    real(SP), dimension(:), pointer :: p_Fmatrix,p_Fx
    real(QP), dimension(:), pointer :: p_QxTimeDeriv, p_QxPredictor
    real(DP), dimension(:), pointer :: p_DxTimeDeriv, p_DxPredictor
    real(SP), dimension(:), pointer :: p_FxTimeDeriv, p_FxPredictor
    real(QP), dimension(:), pointer :: p_Qflux0,p_Qflux,p_QfluxPrel,p_Qalpha
    real(DP), dimension(:), pointer :: p_Dflux0,p_Dflux,p_DfluxPrel,p_Dalpha
    real(SP), dimension(:), pointer :: p_Fflux0,p_Fflux,p_FfluxPrel,p_Falpha
    real(QP), dimension(:,:), pointer :: p_Qcoefficients
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    real(SP), dimension(:,:), pointer :: p_Fcoefficients
    real(QP), dimension(:,:,:), pointer :: p_QcoeffsAtEdge
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(SP), dimension(:,:,:), pointer :: p_FcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    logical :: buseCallback

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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
        call sys_halt()
      end if
    end if

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)

    select case(rx%cdataType)
    case (ST_QUAD)
       call lsyssc_getbase_quad(rx, p_Qx)
    case (ST_DOUBLE)
       call lsyssc_getbase_double(rx, p_Dx)
    case (ST_SINGLE)
       call lsyssc_getbase_single(rx, p_Fx)
    case default
       call output_line('Unsupported data type',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
       call sys_halt()
    end select

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)

      ! Set pointers
      select case(rafcstab%cdataType)
      case(ST_QUAD)
         call lsyssc_getbase_quad(rafcstab%p_rvectorFlux, p_Qflux)
         call lsyssc_getbase_quad(rafcstab%p_rvectorFlux0, p_Qflux0)
      case(ST_DOUBLE)
         call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
         call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
      case(ST_SINGLE)
         call lsyssc_getbase_single(rafcstab%p_rvectorFlux, p_Fflux)
         call lsyssc_getbase_single(rafcstab%p_rvectorFlux0, p_Fflux0)
      case default
         call output_line('Unsupported data type',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
         call sys_halt()
      end select

      ! Use callback routine?
      if (present(fcb_calcFluxFCTSc_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        select case(rgroupFEMSet%cdataType)
        case(ST_QUAD)
           call gfem_getbase_QcoeffsAtEdge(rgroupFEMSet, p_QcoeffsAtEdge)
        case(ST_DOUBLE)
           call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        case(ST_SINGLE)
           call gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)
        case default
           call output_line('Unsupported data type',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
           call sys_halt()
        end select
        buseCallback = .true.
      else

        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        select case(rafcstab%cdataType)
        case(ST_QUAD)
           call afcstab_getbase_QcoeffsAtEdge(rafcstab, p_Qcoefficients)
        case(ST_DOUBLE)
           call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        case(ST_SINGLE)
           call afcstab_getbase_FcoeffsAtEdge(rafcstab, p_Fcoefficients)
        end select
        buseCallback = .false.
      end if

      !-------------------------------------------------------------------------
      ! Classical, iterative and semi-implicit nonlinear FEM-FCT algorithm
      !
      ! The raw antidiffusive fluxes for all algorithms can be
      ! assembled essentially in the same way.
      !
      ! $$ f_{ij}^n = -m_{ij}*(u_i^n-u_j^n)+(1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
      ! $$ f_{ij}^m = f_{ij}^n + m_{ij}*(u_i^m-u_j^m)+\theta\Delta t d_{ij}^m(u_i^m-u_j^m) $$
      !
      ! The only difference is that the amount of rejected  antidiffusion
      ! is subtracted from the initial fluxes in subsequent iterations if
      ! the iterative FEM-FCT algorithm is applied.
      ! Moreover the initial flux without mass contribution is stored
      ! separately for the semi-implicit FEM-FCT algorithm since it is
      ! used to constrain the raw antidiffusive fluxes in each iteration.
      !-------------------------------------------------------------------------

      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_EXPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble explicit part of raw-antidiffive fluxes
        !-----------------------------------------------------------------------

        if (theta .ne. 1.0_DP) then
          ! Assemble the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = (1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE,&
                p_DcoeffsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          else
             
             select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
             case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
                call doFluxesByCoeffsQPQP(p_IedgeList, rafcstab%NEDGE,&
                     p_Qcoefficients, p_Qx, dscale*(1.0_DP-theta),&
                     bclear, p_Qflux0)
             case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
                call doFluxesByCoeffsDPQP(p_IedgeList, rafcstab%NEDGE,&
                     p_Qcoefficients, p_Dx, dscale*(1.0_DP-theta),&
                     bclear, p_Qflux0)
             case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
                call doFluxesByCoeffsSPQP(p_IedgeList, rafcstab%NEDGE,&
                     p_Qcoefficients, p_Fx, dscale*(1.0_DP-theta),&
                     bclear, p_Qflux0)

             case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
                call doFluxesByCoeffsQPDP(p_IedgeList, rafcstab%NEDGE,&
                     p_Dcoefficients, p_Qx, dscale*(1.0_DP-theta),&
                     bclear, p_Dflux0)
             case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
                call doFluxesByCoeffsDPDP(p_IedgeList, rafcstab%NEDGE,&
                     p_Dcoefficients, p_Dx, dscale*(1.0_DP-theta),&
                     bclear, p_Dflux0)
             case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
                call doFluxesByCoeffsSPDP(p_IedgeList, rafcstab%NEDGE,&
                     p_Dcoefficients, p_Fx, dscale*(1.0_DP-theta),&
                     bclear, p_Dflux0)

             case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
                call doFluxesByCoeffsQPSP(p_IedgeList, rafcstab%NEDGE,&
                     p_Fcoefficients, p_Qx, dscale*(1.0_DP-theta),&
                     bclear, p_Fflux0)
             case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
                call doFluxesByCoeffsDPSP(p_IedgeList, rafcstab%NEDGE,&
                     p_Fcoefficients, p_Dx, dscale*(1.0_DP-theta),&
                     bclear, p_Fflux0)
             case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
                call doFluxesByCoeffsSPSP(p_IedgeList, rafcstab%NEDGE,&
                     p_Fcoefficients, p_Fx, dscale*(1.0_DP-theta),&
                     bclear, p_Fflux0)
             end select
          end if
        elseif (.not.bquickAssembly .and. bclear) then
          ! Clear the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = 0 $$
          select case(rafcstab%cdataType)
          case(ST_QUAD)    
             call lalg_clearVector(p_Qflux0, rafcstab%NEDGE)
          case(ST_DOUBLE)    
             call lalg_clearVector(p_Dflux0, rafcstab%NEDGE)
          case(ST_SINGLE)    
             call lalg_clearVector(p_Fflux0, rafcstab%NEDGE)
          end select
          ! if bquickAssembly = TRUE then this step can be skipped
        end if

        !-----------------------------------------------------------------------

        ! Check for special treatment
        if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

          ! Set pointers
          select case(rafcstab%cdataType)
          case(ST_QUAD)
             call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
          case(ST_DOUBLE)
             call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
          case(ST_SINGLE)
             call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
          end select

          ! We have to store the raw-antidiffusive fluxes based on the
          ! initial solution without contribution of the consistent
          ! mass matrix and without scaling by the implicitness parameter
          ! $$ f_{ij} = \Delta t d_{ij}^n(u_i^n-u_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE,&
                p_DcoeffsAtEdge, p_Dx, dscale, .true., p_DfluxPrel)
          else
             select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
             case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
                call doFluxesByCoeffsQPQP(p_IedgeList, rafcstab%NEDGE,&
                     p_Qcoefficients, p_Qx, dscale, .true., p_QfluxPrel)
             case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
                call doFluxesByCoeffsDPQP(p_IedgeList, rafcstab%NEDGE,&
                     p_Qcoefficients, p_Dx, dscale, .true., p_QfluxPrel)
             case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
                call doFluxesByCoeffsSPQP(p_IedgeList, rafcstab%NEDGE,&
                     p_Qcoefficients, p_Fx, dscale, .true., p_QfluxPrel)
                
             case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
                call doFluxesByCoeffsQPDP(p_IedgeList, rafcstab%NEDGE,&
                     p_Dcoefficients, p_Qx, dscale, .true., p_DfluxPrel)
             case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
                call doFluxesByCoeffsDPDP(p_IedgeList, rafcstab%NEDGE,&
                     p_Dcoefficients, p_Dx, dscale, .true., p_DfluxPrel)
             case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
                call doFluxesByCoeffsSPDP(p_IedgeList, rafcstab%NEDGE,&
                     p_Dcoefficients, p_Fx, dscale, .true., p_DfluxPrel)
                
             case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
                call doFluxesByCoeffsQPSP(p_IedgeList, rafcstab%NEDGE,&
                     p_Fcoefficients, p_Qx, dscale, .true., p_FfluxPrel)
             case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
                call doFluxesByCoeffsDPSP(p_IedgeList, rafcstab%NEDGE,&
                     p_Fcoefficients, p_Dx, dscale, .true., p_FfluxPrel)
             case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
                call doFluxesByCoeffsSPSP(p_IedgeList, rafcstab%NEDGE,&
                     p_Fcoefficients, p_Fx, dscale, .true., p_FfluxPrel)
             end select
          end if

        elseif (rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) then

          ! We have to assemble the raw-antidiffusive fluxes for
          ! prelimiting separately based on the low-order predictor
          if (present(rxPredictor)) then

            ! Set pointers
            select case(rxPredictor%cdataType)
            case(ST_QUAD)
               call lsyssc_getbase_quad(rxPredictor, p_QxPredictor)
            case(ST_DOUBLE)
               call lsyssc_getbase_double(rxPredictor, p_DxPredictor)
            case(ST_SINGLE)
               call lsyssc_getbase_single(rxPredictor, p_FxPredictor)
            case default
               call output_line('Unsupported data type',&
                    OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
               call sys_halt()
            end select

            select case(rafcstab%cdataType)
            case(ST_QUAD)
               call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
            case(ST_DOUBLE)
               call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
            case(ST_SINGLE)
               call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
            end select
            
            if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
               select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
               case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
                  call doDifferencesQPQP(p_IedgeList, rafcstab%NEDGE,&
                       p_QxPredictor, p_QfluxPrel)
               case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
                  call doDifferencesDPQP(p_IedgeList, rafcstab%NEDGE,&
                       p_DxPredictor, p_QfluxPrel)
               case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
                  call doDifferencesSPQP(p_IedgeList, rafcstab%NEDGE,&
                       p_FxPredictor, p_QfluxPrel)

               case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
                  call doDifferencesQPDP(p_IedgeList, rafcstab%NEDGE,&
                       p_QxPredictor, p_DfluxPrel)
               case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
                  call doDifferencesDPDP(p_IedgeList, rafcstab%NEDGE,&
                       p_DxPredictor, p_DfluxPrel)
               case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
                  call doDifferencesSPDP(p_IedgeList, rafcstab%NEDGE,&
                       p_FxPredictor, p_DfluxPrel)

               case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
                  call doDifferencesQPSP(p_IedgeList, rafcstab%NEDGE,&
                       p_QxPredictor, p_FfluxPrel)
               case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
                  call doDifferencesDPSP(p_IedgeList, rafcstab%NEDGE,&
                       p_DxPredictor, p_FfluxPrel)
               case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
                  call doDifferencesSPSP(p_IedgeList, rafcstab%NEDGE,&
                       p_FxPredictor, p_FfluxPrel)
               end select
            elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
              ! Compute fluxes for minmod prelimiting
              if (buseCallback) then
                call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE,&
                    p_DcoeffsAtEdge, p_DxPredictor, dscale, .true., p_DfluxPrel)
              else
                 select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
                 case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
                    call doFluxesByCoeffsQPQP(p_IedgeList, rafcstab%NEDGE,&
                         p_Qcoefficients, p_QxPredictor, dscale, .true., p_QfluxPrel)
                 case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
                    call doFluxesByCoeffsDPQP(p_IedgeList, rafcstab%NEDGE,&
                         p_Qcoefficients, p_DxPredictor, dscale, .true., p_QfluxPrel)
                 case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
                    call doFluxesByCoeffsSPQP(p_IedgeList, rafcstab%NEDGE,&
                         p_Qcoefficients, p_FxPredictor, dscale, .true., p_QfluxPrel)

                 case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
                    call doFluxesByCoeffsQPDP(p_IedgeList, rafcstab%NEDGE,&
                         p_Dcoefficients, p_QxPredictor, dscale, .true., p_DfluxPrel)
                 case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
                    call doFluxesByCoeffsDPDP(p_IedgeList, rafcstab%NEDGE,&
                         p_Dcoefficients, p_DxPredictor, dscale, .true., p_DfluxPrel)
                 case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
                    call doFluxesByCoeffsSPDP(p_IedgeList, rafcstab%NEDGE,&
                         p_Dcoefficients, p_FxPredictor, dscale, .true., p_DfluxPrel)

                 case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
                    call doFluxesByCoeffsQPSP(p_IedgeList, rafcstab%NEDGE,&
                         p_Fcoefficients, p_QxPredictor, dscale, .true., p_FfluxPrel)
                 case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
                    call doFluxesByCoeffsDPSP(p_IedgeList, rafcstab%NEDGE,&
                         p_Fcoefficients, p_DxPredictor, dscale, .true., p_FfluxPrel)
                 case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
                    call doFluxesByCoeffsSPSP(p_IedgeList, rafcstab%NEDGE,&
                         p_Fcoefficients, p_FxPredictor, dscale, .true., p_FfluxPrel)
                 end select
              end if
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
            call sys_halt()
          end if
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then

          ! Set pointers
          select case(rmatrix%cdataType)
          case(ST_QUAD)
             call lsyssc_getbase_quad(rmatrix, p_Qmatrix)
          case(ST_DOUBLE)
             call lsyssc_getbase_double(rmatrix, p_Dmatrix)
          case(ST_SINGLE)
             call lsyssc_getbase_single(rmatrix, p_Fmatrix)
          case default
             call output_line('Unsupported data type',&
                  OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
             call sys_halt()
          end select

          ! Assemble the explicit part of the mass-antidiffusive fluxes
          ! $$ f_{ij}^n := f_{ij}^n - m_{ij}(u_i^n-u_j^n) $$
          select case(FEAT2_PP_ID3(rmatrix%cdataType,rx%cdataType,rafcstab%cdataType,10))
             ! AFC = QUAD
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixQPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, -dscale/tstep, .false., p_Qflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixDPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, -dscale/tstep, .false., p_Qflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixSPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, -dscale/tstep, .false., p_Qflux0)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixQPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, -dscale/tstep, .false., p_Qflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixDPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, -dscale/tstep, .false., p_Qflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixSPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, -dscale/tstep, .false., p_Qflux0)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixQPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, -dscale/tstep, .false., p_Qflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixDPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, -dscale/tstep, .false., p_Qflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixSPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, -dscale/tstep, .false., p_Qflux0)


             ! AFC = DOUBLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixQPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, -dscale/tstep, .false., p_Dflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixDPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, -dscale/tstep, .false., p_Dflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixSPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, -dscale/tstep, .false., p_Dflux0)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, -dscale/tstep, .false., p_Dflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, -dscale/tstep, .false., p_Dflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, -dscale/tstep, .false., p_Dflux0)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, -dscale/tstep, .false., p_Dflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, -dscale/tstep, .false., p_Dflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, -dscale/tstep, .false., p_Dflux0)


             ! AFC = SINGLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixQPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, -dscale/tstep, .false., p_Fflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixDPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, -dscale/tstep, .false., p_Fflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixSPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, -dscale/tstep, .false., p_Fflux0)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixQPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, -dscale/tstep, .false., p_Fflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixDPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, -dscale/tstep, .false., p_Fflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixSPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, -dscale/tstep, .false., p_Fflux0)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixQPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, -dscale/tstep, .false., p_Fflux0)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixDPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, -dscale/tstep, .false., p_Fflux0)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixSPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, -dscale/tstep, .false., p_Fflux0)
          end select
        end if

      end if


      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_IMPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble implicit part of raw-antidiffusive fluxes
        !-----------------------------------------------------------------------

        if ((rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_ITERATIVE) .and.&
            iand(ioperationSpec, AFCSTAB_FCTFLUX_REJECTED) .ne. 0) then
          !---------------------------------------------------------------------
          ! Apply the rejected antidiffusive fluxes from the previous limiting
          ! step to the implicit part of the raw-antidiffusive fluxes
          ! --------------------------------------------------------------------

          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
            call output_line('Stabilisation does not provide correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
            call sys_halt()
          end if

          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
            call sys_halt()
          end if

          ! Set pointer
          select case(rafcstab%cdataType)
          case(ST_QUAD)
             call lsyssc_getbase_quad(rafcstab%p_rvectorAlpha, p_Qalpha)
          case(ST_DOUBLE)
             call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
          case(ST_SINGLE)
             call lsyssc_getbase_single(rafcstab%p_rvectorAlpha, p_Falpha)
          end select

          ! Subtract amount of rejected antidiffusion
          select case(rafcstab%cdataType)
          case(ST_QUAD)
             call afcstab_combineFluxes(rafcstab%NEDGE, -1.0_QP, p_Qflux, p_Qflux0, p_Qalpha)
          case(ST_DOUBLE)
             call afcstab_combineFluxes(rafcstab%NEDGE, -1.0_DP, p_Dflux, p_Dflux0, p_Dalpha)
          case(ST_SINGLE)
             call afcstab_combineFluxes(rafcstab%NEDGE, -1.0_SP, p_Fflux, p_Fflux0, p_Falpha)
          end select
        end if

        !-----------------------------------------------------------------------

        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij} = \theta\Delta t d_{ij}(u_i-u_j) $$
           if (buseCallback) then
              call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE,&
                   p_DcoeffsAtEdge, p_Dx, dscale*theta, bclear, p_Dflux)
           else
              select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
              case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
                 call doFluxesByCoeffsQPQP(p_IedgeList, rafcstab%NEDGE,&
                      p_Qcoefficients, p_Qx, dscale*theta, bclear, p_Qflux)
              case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
                 call doFluxesByCoeffsDPQP(p_IedgeList, rafcstab%NEDGE,&
                      p_Qcoefficients, p_Dx, dscale*theta, bclear, p_Qflux)
              case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
                 call doFluxesByCoeffsSPQP(p_IedgeList, rafcstab%NEDGE,&
                      p_Qcoefficients, p_Fx, dscale*theta, bclear, p_Qflux)

              case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
                 call doFluxesByCoeffsQPDP(p_IedgeList, rafcstab%NEDGE,&
                      p_Dcoefficients, p_Qx, dscale*theta, bclear, p_Dflux)
              case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
                 call doFluxesByCoeffsDPDP(p_IedgeList, rafcstab%NEDGE,&
                      p_Dcoefficients, p_Dx, dscale*theta, bclear, p_Dflux)
              case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
                 call doFluxesByCoeffsSPDP(p_IedgeList, rafcstab%NEDGE,&
                      p_Dcoefficients, p_Fx, dscale*theta, bclear, p_Dflux)

              case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
                 call doFluxesByCoeffsQPSP(p_IedgeList, rafcstab%NEDGE,&
                      p_Fcoefficients, p_Qx, dscale*theta, bclear, p_Fflux)
              case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
                 call doFluxesByCoeffsDPSP(p_IedgeList, rafcstab%NEDGE,&
                      p_Fcoefficients, p_Dx, dscale*theta, bclear, p_Fflux)
              case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
                 call doFluxesByCoeffsSPSP(p_IedgeList, rafcstab%NEDGE,&
                      p_Fcoefficients, p_Fx, dscale*theta, bclear, p_Fflux)
              end select
           end if
        end if

        if (bquickAssembly) then
          ! We may check of either the implicit or explicit part are
          ! missing so that some redundant computations may be skipped
          if (theta .ne. 1.0_DP) then
            ! The explicit part of the raw-antidiffusive fluxes exists
            if (theta .ne. 0.0_DP) then
              ! The implicit part of the raw-antidiffusive fluxes
              ! exists; so combine them both into common fluxes
               select case(rafcstab%cdataType)
               case(ST_QUAD)
                  call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_QP, p_Qflux0, p_Qflux)
               case(ST_DOUBLE)
                  call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_DP, p_Dflux0, p_Dflux)
               case(ST_SINGLE)
                  call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_SP, p_Fflux0, p_Fflux)
               end select

            else
              ! The implicit part of the raw-antidiffusive fluxes does
              ! not exists; the fluxes should be cleared so just
              ! overwrite them by the explicit part
              select case(rafcstab%cdataType)
              case(ST_QUAD)
                 call lalg_copyVector(p_Qflux0, p_Qflux)
              case(ST_DOUBLE)
                 call lalg_copyVector(p_Dflux0, p_Dflux)
              case(ST_SINGLE)
                 call lalg_copyVector(p_Fflux0, p_Fflux)
              end select
            end if
            ! if theta = 1 then the explicit part does not exist
          end if
        else
          ! Truely combine both parts of the raw-antidiffusive fluxes
           select case(rafcstab%cdataType)
           case(ST_QUAD)
              call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_QP, p_Qflux0, p_Qflux)
           case(ST_DOUBLE)
              call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_DP, p_Dflux0, p_Dflux)
           case(ST_SINGLE)
              call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_SP, p_Fflux0, p_Fflux)
           end select
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then

          ! Set pointers
           select case(rmatrix%cdataType)
           case(ST_QUAD)
              call lsyssc_getbase_quad(rmatrix, p_Qmatrix)
           case(ST_DOUBLE)
              call lsyssc_getbase_double(rmatrix, p_Dmatrix)
           case(ST_SINGLE)
              call lsyssc_getbase_single(rmatrix, p_Fmatrix)
           case default
              call output_line('Unsupported data type',&
                   OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
              call sys_halt()
           end select

          ! Assemble the implicit part of the mass-antidiffusive fluxes
          ! $$ f_{ij}^m := f_{ij}^m + m_{ij}(u_i^m-u_j^m) $$

          select case(FEAT2_PP_ID3(rmatrix%cdataType,rx%cdataType,rafcstab%cdataType,10))
             ! AFC = QUAD
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixQPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, dscale/tstep, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixDPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, dscale/tstep, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixSPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, dscale/tstep, .false., p_Qflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixQPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, dscale/tstep, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixDPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, dscale/tstep, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixSPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, dscale/tstep, .false., p_Qflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixQPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, dscale/tstep, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixDPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, dscale/tstep, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixSPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, dscale/tstep, .false., p_Qflux)


             ! AFC = DOUBLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixQPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, dscale/tstep, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixDPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, dscale/tstep, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixSPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, dscale/tstep, .false., p_Dflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, dscale/tstep, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, dscale/tstep, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, dscale/tstep, .false., p_Dflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, dscale/tstep, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, dscale/tstep, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, dscale/tstep, .false., p_Dflux)


             ! AFC = SINGLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixQPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, dscale/tstep, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixDPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, dscale/tstep, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixSPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, dscale/tstep, .false., p_Fflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixQPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, dscale/tstep, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixDPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, dscale/tstep, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixSPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, dscale/tstep, .false., p_Fflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixQPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, dscale/tstep, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixDPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, dscale/tstep, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixSPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, dscale/tstep, .false., p_Fflux)
          end select

        end if

      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_LINFCT)

      !-------------------------------------------------------------------------
      ! Linearised FEM-FCT algorithm
      !-------------------------------------------------------------------------

      ! Set pointer
      select case(rafcstab%cdataType)
      case(ST_QUAD)
         call lsyssc_getbase_quad(rafcstab%p_rvectorFlux, p_Qflux)
      case(ST_DOUBLE)
         call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      case(ST_SINGLE)
         call lsyssc_getbase_single(rafcstab%p_rvectorFlux, p_Fflux)
      end select

      ! Use callback routine?
      if (present(fcb_calcFluxFCTSc_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        select case(rgroupFEMSet%cdataType)
        case(ST_QUAD)
           call gfem_getbase_QcoeffsAtEdge(rgroupFEMSet, p_QcoeffsAtEdge)
        case(ST_DOUBLE)
           call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        case(ST_SINGLE)
           call gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)
        case default
           call output_line('Unsupported data type',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
           call sys_halt()
        end select
        buseCallback = .true.
      else

        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if

        select case(rafcstab%cdataType)
        case(ST_QUAD)
           call afcstab_getbase_QcoeffsAtEdge(rafcstab, p_Qcoefficients)
        case(ST_DOUBLE)
           call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        case(ST_SINGLE)
           call afcstab_getbase_FcoeffsAtEdge(rafcstab, p_Fcoefficients)
        end select
        buseCallback = .false.
      end if

      ! Assemble spatial part of raw-antidiffusive fluxes
      if (buseCallback) then
        call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE,&
            p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Dflux)
      else
         
         select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
         case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
            call doFluxesByCoeffsQPQP(p_IedgeList, rafcstab%NEDGE,&
                 p_Qcoefficients, p_Qx, dscale, bclear, p_Qflux)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
            call doFluxesByCoeffsDPQP(p_IedgeList, rafcstab%NEDGE,&
                 p_Qcoefficients, p_Dx, dscale, bclear, p_Qflux)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
            call doFluxesByCoeffsSPQP(p_IedgeList, rafcstab%NEDGE,&
                 p_Qcoefficients, p_Fx, dscale, bclear, p_Qflux)

         case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
            call doFluxesByCoeffsQPDP(p_IedgeList, rafcstab%NEDGE,&
                 p_Dcoefficients, p_Qx, dscale, bclear, p_Dflux)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
            call doFluxesByCoeffsDPDP(p_IedgeList, rafcstab%NEDGE,&
                 p_Dcoefficients, p_Dx, dscale, bclear, p_Dflux)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
            call doFluxesByCoeffsSPDP(p_IedgeList, rafcstab%NEDGE,&
                 p_Dcoefficients, p_Fx, dscale, bclear, p_Dflux)

         case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
            call doFluxesByCoeffsQPSP(p_IedgeList, rafcstab%NEDGE,&
                 p_Fcoefficients, p_Qx, dscale, bclear, p_Fflux)
         case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
            call doFluxesByCoeffsDPSP(p_IedgeList, rafcstab%NEDGE,&
                 p_Fcoefficients, p_Dx, dscale, bclear, p_Fflux)
         case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
            call doFluxesByCoeffsSPSP(p_IedgeList, rafcstab%NEDGE,&
                 p_Fcoefficients, p_Fx, dscale, bclear, p_Fflux)
         end select

      end if

      !-------------------------------------------------------------------------

      if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
        ! Compute fluxes for standard prelimiting based on the
        ! low-order solution which serves as predictor
        select case(rafcstab%cdataType)
        case(ST_QUAD)
           call lsyssc_getbase_quad(rafcstab%p_rvectorFluxPrel, p_QfluxPrel)
        case(ST_DOUBLE)
           call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        case(ST_SINGLE)
           call lsyssc_getbase_single(rafcstab%p_rvectorFluxPrel, p_FfluxPrel)
        end select


        select case(FEAT2_PP_ID2(rx%cdataType,rafcstab%cdataType,10))
        case(FEAT2_PP_ID2(ST_QUAD,ST_QUAD,10))
           call doDifferencesQPQP(p_IedgeList, rafcstab%NEDGE,&
                p_Qx, p_QfluxPrel)
        case(FEAT2_PP_ID2(ST_DOUBLE,ST_QUAD,10))
           call doDifferencesDPQP(p_IedgeList, rafcstab%NEDGE,&
                p_Dx, p_QfluxPrel)
        case(FEAT2_PP_ID2(ST_SINGLE,ST_QUAD,10))
           call doDifferencesSPQP(p_IedgeList, rafcstab%NEDGE,&
                p_Fx, p_QfluxPrel)

        case(FEAT2_PP_ID2(ST_QUAD,ST_DOUBLE,10))
           call doDifferencesQPDP(p_IedgeList, rafcstab%NEDGE,&
                p_Qx, p_DfluxPrel)
        case(FEAT2_PP_ID2(ST_DOUBLE,ST_DOUBLE,10))
           call doDifferencesDPDP(p_IedgeList, rafcstab%NEDGE,&
                p_Dx, p_DfluxPrel)
        case(FEAT2_PP_ID2(ST_SINGLE,ST_DOUBLE,10))
           call doDifferencesSPDP(p_IedgeList, rafcstab%NEDGE,&
                p_Fx, p_DfluxPrel)

        case(FEAT2_PP_ID2(ST_QUAD,ST_SINGLE,10))
           call doDifferencesQPSP(p_IedgeList, rafcstab%NEDGE,&
                p_Qx, p_FfluxPrel)
        case(FEAT2_PP_ID2(ST_DOUBLE,ST_SINGLE,10))
           call doDifferencesDPSP(p_IedgeList, rafcstab%NEDGE,&
                p_Dx, p_FfluxPrel)
        case(FEAT2_PP_ID2(ST_SINGLE,ST_SINGLE,10))
           call doDifferencesSPSP(p_IedgeList, rafcstab%NEDGE,&
                p_Fx, p_FfluxPrel)
        end select


      elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
        ! Make a backup of the spatial part of the raw-antidiffusive
        ! fluxes which are used for minmod prelimiting
        call lsyssc_copyVector(rafcstab%p_rvectorFlux,&
            rafcstab%p_rvectorFluxPrel)
      end if

      !-------------------------------------------------------------------------

      ! Do we have to include mass antidiffusion?
      if (present(rmatrix) .and. present(rxTimeDeriv)) then

        ! Set pointer
        select case(rmatrix%cdataType)
        case(ST_QUAD)
           call lsyssc_getbase_quad(rmatrix, p_Qmatrix)
        case(ST_DOUBLE)
           call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        case(ST_SINGLE)
           call lsyssc_getbase_single(rmatrix, p_Fmatrix)
        case default
           call output_line('Unsupported data type',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
           call sys_halt()
        end select

        select case(rxTimeDeriv%cdataType)
        case(ST_QUAD)
           call lsyssc_getbase_quad(rxTimeDeriv, p_QxTimeDeriv)
        case(ST_DOUBLE)
           call lsyssc_getbase_double(rxTimeDeriv, p_DxTimeDeriv)
        case(ST_SINGLE)
           call lsyssc_getbase_single(rxTimeDeriv, p_FxTimeDeriv)
        case default
           call output_line('Unsupported data type',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
           call sys_halt()
        end select

        ! Apply mass antidiffusion to antidiffusive fluxes based on
        ! the approximation to the time derivative

          select case(FEAT2_PP_ID3(rmatrix%cdataType,rx%cdataType,rafcstab%cdataType,10))
             ! AFC = QUAD
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixQPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_QxTimeDeriv, dscale, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixDPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_QxTimeDeriv, dscale, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixSPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_QxTimeDeriv, dscale, .false., p_Qflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixQPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_DxTimeDeriv, dscale, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixDPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_DxTimeDeriv, dscale, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixSPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_DxTimeDeriv, dscale, .false., p_Qflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixQPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_FxTimeDeriv, dscale, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixDPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_FxTimeDeriv, dscale, .false., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixSPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_FxTimeDeriv, dscale, .false., p_Qflux)


             ! AFC = DOUBLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixQPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_QxTimeDeriv, dscale, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixDPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_QxTimeDeriv, dscale, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixSPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_QxTimeDeriv, dscale, .false., p_Dflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_DxTimeDeriv, dscale, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_DxTimeDeriv, dscale, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_DxTimeDeriv, dscale, .false., p_Dflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_FxTimeDeriv, dscale, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_FxTimeDeriv, dscale, .false., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_FxTimeDeriv, dscale, .false., p_Dflux)


             ! AFC = SINGLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixQPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_QxTimeDeriv, dscale, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixDPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_QxTimeDeriv, dscale, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixSPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_QxTimeDeriv, dscale, .false., p_Fflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixQPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_DxTimeDeriv, dscale, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixDPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_DxTimeDeriv, dscale, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixSPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_DxTimeDeriv, dscale, .false., p_Fflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixQPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_FxTimeDeriv, dscale, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixDPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_FxTimeDeriv, dscale, .false., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixSPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_FxTimeDeriv, dscale, .false., p_Fflux)
          end select

      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)


    case (AFCSTAB_LINFCT_MASS)

      !-------------------------------------------------------------------------
      ! FEM-FCT algorithm for mass antidiffusion
      !-------------------------------------------------------------------------

      if (present(rmatrix)) then

        ! Set pointers
        select case(rmatrix%cdataType)
        case(ST_QUAD)
           call lsyssc_getbase_quad(rmatrix, p_Qmatrix)
        case(ST_DOUBLE)
           call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        case(ST_SINGLE)
           call lsyssc_getbase_single(rmatrix, p_Fmatrix)
        case default
           call output_line('Unsupported data type',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
           call sys_halt()
        end select

        select case(rafcstab%cdataType)
        case(ST_QUAD)
           call lsyssc_getbase_quad(rafcstab%p_rvectorFlux, p_Qflux)
        case(ST_DOUBLE)
           call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        case(ST_SINGLE)
           call lsyssc_getbase_single(rafcstab%p_rvectorFlux, p_Fflux)
        end select
        
        ! Clear vector and assemble antidiffusive fluxes

        select case(FEAT2_PP_ID3(rmatrix%cdataType,rx%cdataType,rafcstab%cdataType,10))
             ! AFC = QUAD
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixQPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, dscale, .true., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixDPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, dscale, .true., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_QUAD,10))
             call doFluxesByMatrixSPQPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, dscale, .true., p_Qflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixQPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, dscale, .true., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixDPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, dscale, .true., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_QUAD,10))
             call doFluxesByMatrixSPDPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, dscale, .true., p_Qflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixQPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, dscale, .true., p_Qflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixDPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, dscale, .true., p_Qflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_QUAD,10))
             call doFluxesByMatrixSPSPQP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, dscale, .true., p_Qflux)


             ! AFC = DOUBLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixQPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, dscale, .true., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixDPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, dscale, .true., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_DOUBLE,10))
             call doFluxesByMatrixSPQPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, dscale, .true., p_Dflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, dscale, .true., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, dscale, .true., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPDPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, dscale, .true., p_Dflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixQPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, dscale, .true., p_Dflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixDPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, dscale, .true., p_Dflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_DOUBLE,10))
             call doFluxesByMatrixSPSPDP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, dscale, .true., p_Dflux)


             ! AFC = SINGLE
          case(FEAT2_PP_ID3(ST_QUAD,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixQPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Qx, dscale, .true., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixDPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Qx, dscale, .true., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_QUAD,ST_SINGLE,10))
             call doFluxesByMatrixSPQPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Qx, dscale, .true., p_Fflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixQPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Dx, dscale, .true., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixDPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Dx, dscale, .true., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_DOUBLE,ST_SINGLE,10))
             call doFluxesByMatrixSPDPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Dx, dscale, .true., p_Fflux)

          case(FEAT2_PP_ID3(ST_QUAD,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixQPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Qmatrix, p_Fx, dscale, .true., p_Fflux)
          case(FEAT2_PP_ID3(ST_DOUBLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixDPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Dmatrix, p_Fx, dscale, .true., p_Fflux)
          case(FEAT2_PP_ID3(ST_SINGLE,ST_SINGLE,ST_SINGLE,10))
             call doFluxesByMatrixSPSPSP(p_IedgeList, rafcstab%NEDGE,&
                  p_Fmatrix, p_Fx, dscale, .true., p_Fflux)
          end select


        ! Set specifiers for raw antidiffusive fluxes
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
        call sys_halt()
      end if


    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the edge-by-edge array DcoefficientsAtEdge

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByCoeffs.h"
#undef TemplateType_Vector
#undef TemplateType_AFC




    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix stored in Dmatrix


#define TemplateType_Matrix QUAD_PREC
#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC
#undef TemplateType_Matrix


#define TemplateType_Matrix DOUBLE_PREC
#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC
#undef TemplateType_Matrix


#define TemplateType_Matrix SINGLE_PREC
#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doFluxesByMatrix.h"
#undef TemplateType_Vector
#undef TemplateType_AFC
#undef TemplateType_Matrix




    !**************************************************************
    ! Assemble raw antidiffusive fluxes with aid of callback function

    subroutine doFluxesByCallbackDP(IedgeList, NEDGE,&
         DcoeffsAtEdge, Dx, dscale, bclear, Dflux)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:), pointer :: DfluxAtEdge

      ! local variables
      integer :: idx,iedge,IEDGEset,IEDGEmax


      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (bclear) then

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,idx,iedge,IEDGEmax)

        ! Allocate temporal memory
        allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(NEDGE, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dflux(IEDGEset:IEDGEmax), rcollection)
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        !$omp end parallel

      else   ! bclear = .false.

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,DfluxAtEdge,idx,iedge,IEDGEmax)

        ! Allocate temporal memory
        allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
        allocate(DfluxAtEdge(p_rperfconfig%NEDGESIM))

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(NEDGE, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxAtEdge(1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Add antidiffusive fluxes
            Dflux(iedge) = Dflux(iedge) + DfluxAtEdge(idx)
          end do
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge, DfluxAtEdge)
        !$omp end parallel

      end if

    end subroutine doFluxesByCallbackDP



    

    !**************************************************************
    ! Assemble solution difference used for classical prelimiting.

#define TemplateType_AFC    QUAD_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    DOUBLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#undef TemplateType_AFC

#define TemplateType_AFC    SINGLE_PREC
#define TemplateType_Vector QUAD_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#define TemplateType_Vector DOUBLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#define TemplateType_Vector SINGLE_PREC
#include "kernel/PDEOperators/afcsc_buildFluxFCTScalar_doDifferences.h"
#undef TemplateType_Vector
#undef TemplateType_AFC



  end subroutine afcsc_buildFluxFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearFCTBlock(rx, theta, tstep, hstep,&
      bclear, rafcstab, rjacobian, rmatrix)

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
    type(t_matrixScalar), intent(in), optional :: rmatrix
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTBlock')
      call sys_halt()

    else

      call afcsc_buildJacLinearFCTScalar(&
          rx%RvectorBlock(1), theta, tstep, hstep, bclear,&
          rafcstab, rjacobian, rmatrix)

    end if
  end subroutine afcsc_buildJacLinearFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearFCTScalar(rx, theta, tstep, hstep,&
      bclear, rafcstab, rjacobian, rmatrix)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation. Note that the velocity is assumed
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
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dflux,p_Dflux0
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_Dx


    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)   .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)


    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_NLINFCT_IMPLICIT)

      ! What kind of matrix are we?
      select case(rjacobian%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        !-------------------------------------------------------------------------
        ! Matrix format 7
        !-------------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobian, p_Kld)

        if (present(rmatrix)) then
          call lsyssc_getbase_double(rmatrix, p_MC)
          call doJacobian_implFCTconsMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kld, p_MC, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kld, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        end if

      case(LSYSSC_MATRIX9)
        !-------------------------------------------------------------------------
        ! Matrix format 9
        !-------------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobian, p_Kld)
        call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)

        if (present(rmatrix)) then
          call lsyssc_getbase_double(rmatrix, p_MC)
          call doJacobian_implFCTconsMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal, p_MC, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        end if


      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTScalar')
        call sys_halt()
      end select

    case default
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for semi-implicit FEM-FCT,
    ! whereby no mass antidiffusion is built into the matrix

    subroutine doJacobian_implFCTnoMass(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

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

    subroutine doJacobian_implFCTconsMass(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, MC, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

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

  end subroutine afcsc_buildJacLinearFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianFCTBlock(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, rmatrix, rcollection, rperfconfig)

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
    include 'intf_calcMatrixSc_sim.inc'

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTBlock')
      call sys_halt()

    else

      call afcsc_buildJacobianFCTScalar(rgroupFEMSet,&
          rx%RvectorBlock(1), fcb_calcMatrixSc_sim,&
          theta, tstep, hstep, bclear, rafcstab, rjacobian,&
          rmatrix, rcollection, rperfconfig)

    end if
  end subroutine afcsc_buildJacobianFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianFCTScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, rmatrix, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation. The velocity is assumed to be
    ! nonlinear/arbitrary. This routine will also work for linear
    ! velocities but then it is inefficient since the solution
    ! perturbation does not affect the velocity.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

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
    include 'intf_calcMatrixSc_sim.inc'

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

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
    real(DP), dimension(:), pointer :: p_Dflux,p_Dflux0,p_Dx,p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_MC,p_Jac
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer :: ndim


    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)   .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
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
!!$          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
!!$      call sys_halt()
!!$    end select

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_NLINFCT_IMPLICIT)

      ! What kind of matrix format are we?
      select case(rjacobian%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        !-----------------------------------------------------------------------
        ! Matrix format 7
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobian, p_Kld)

        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if

        case (NDIM2D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE,  p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if

        case (NDIM3D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
        end select


      case(LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)

        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if

        case (NDIM2D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if

        case (NDIM3D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
        call sys_halt()
      end select

    case default
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 1D,
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_1D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep


        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
    subroutine doJacobian_implFCTconsMass_1D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, MC, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep

        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
    subroutine doJacobian_implFCTnoMass_2D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in)   :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep


        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
    subroutine doJacobian_implFCTconsMass_2D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, MC, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in)   :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep

        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
    subroutine doJacobian_implFCTnoMass_3D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep


        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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
    subroutine doJacobian_implFCTconsMass_3D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, MC, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)

        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)

        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep

        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

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

  end subroutine afcsc_buildJacobianFCTScalar

end module afcstabscalarfct

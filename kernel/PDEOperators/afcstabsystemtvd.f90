!##############################################################################
!# ****************************************************************************
!# <name> afcstabsystemtvd </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the TVD-type
!# algebraic flux correction methodology proposed by Kuzmin, Moeller
!# and Turek in a series of publications. As a starting point for
!# systems of conservation laws, the reader is referred to the book
!# chapter
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
!# afcstabbase.
!#
!# The following routines are available:
!#
!# 1.) afcsys_buildVectorTVD = afcsys_buildVectorTVDScalar /
!#                             afcsys_buildVectorTVDBlock
!#     -> Assembles the vector for AFC stabilisation of TVD type
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) minmod QP / minmodDP / minmodSP
!#     -> Computes minmod(a,b)
!#
!# </purpose>
!##############################################################################

module afcstabsystemtvd

!$use omp_lib
  use afcstabbase
  use afcstabsystem
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

  public :: afcsys_buildVectorTVD

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

  interface afcsys_buildVectorTVD
    module procedure afcsys_buildVectorTVDScalar
    module procedure afcsys_buildVectorTVDBlock
  end interface

contains

  !****************************************************************************

!<subroutine>

  subroutine afcsys_buildVectorTVDBlock(rafcstab, rgroupFEMSet, rx, ndim,&
      fcb_calcFluxSys_sim, fcb_calcCharacteristics_sim, dscale, bclear, ry,&
      rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector for FEM-TVD schemes. If
    ! the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! number of spatial dimensions
    integer, intent(in) :: ndim

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback function to compute local fluxes
    include 'intf_calcFluxSys_sim.inc'

    ! callback function to compute local characteristics
    include 'intf_calcCharacteristics_sim.inc'

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call afcsys_buildVectorTVDScalar(rafcstab, rgroupFEMSet, rx%RvectorBlock(1),&
          ndim, fcb_calcFluxSys_sim, fcb_calcCharacteristics_sim, dscale, bclear,&
          ry%RvectorBlock(1), rcollection, rperfconfig)
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST).eq. 0) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDBlock')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)

    ! Set pointers
    call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)

    ! How many dimensions do we have?
    select case(ndim)
    case (NDIM1D)
      call doLimitTVD_1D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM2D)
      call doLimitTVD_2D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM3D)
      call doLimitTVD_3D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D

    subroutine doLimitTVD_1D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

        ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_1D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D

    subroutine doLimitTVD_2D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_2D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D

    subroutine doLimitTVD_3D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(3, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Z-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(3, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_3D


    !**************************************************************
    ! Assemble the upper and lower bounds Q and the sums of
    ! antidiffusive contributions P for a given set of edges

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doBoundsAndIncrements_sim(idirection, NVAR,&
        DcoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IedgeList, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dpp,Dpm,Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      integer :: idx,ivar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,i,ivar,j)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)

        ! Loop over all characteristic variables
        do ivar = 1, NVAR
          ! Set node orientation
          if (Daux1(ivar) .gt. 0) then
            i = IedgeList(2,idx)
            j = IedgeList(1,idx)
            Dflux(ivar) = -Dflux(ivar)
          else
            i = IedgeList(1,idx)
            j = IedgeList(2,idx)
          end if

          ! Assemble P's and Q's
          if (Dflux(ivar) .gt. 0) then
            Dpp(ivar,i) = Dpp(ivar,i)+Dflux(ivar)
            Dqm(ivar,i) = Dqm(ivar,i)-Dflux(ivar)
            Dqp(ivar,j) = Dqp(ivar,j)+Dflux(ivar)
          else
            Dpm(ivar,i) = Dpm(ivar,i)+Dflux(ivar)
            Dqp(ivar,i) = Dqp(ivar,i)-Dflux(ivar)
            Dqm(ivar,j) = Dqm(ivar,j)+Dflux(ivar)
          end if
        end do

      end do
      !$omp end parallel do

    end subroutine doBoundsAndIncrements_sim


    !**************************************************************
    ! Limit the antidiffusive fluxes and apply them to the vector

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doLimitADFluxes_sim(idirection, NVAR, dscale, Drp, Drm,&
        DcoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IedgeList, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      real(DP), dimension(:,:), intent(in) :: DrighteigenvectorsAtEdge,Drp,Drm
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      real(DP) :: daux
      integer :: idx,ivar,jvar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,daux,i,ivar,j,jvar)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)
        Daux2 = abs(Daux1)*DcharVariablesAtEdge(:,idx)

        ! Get position of nodes
        i = IedgeList(1,idx)
        j = IedgeList(2,idx)

        ! Loop over all characteristic variables
        ! and limit characteristic fluxes
        do ivar = 1, NVAR

          if (Daux1(ivar) .lt. 0) then
            if (Dflux(ivar) .gt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,i)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,i)*Dflux(ivar)
            end if
          else
            if (Dflux(ivar) .lt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,j)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,j)*Dflux(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        do ivar = 1, NVAR
          daux = 0.0_DP
          do jvar = 1, NVAR
            daux = daux+DrighteigenvectorsAtEdge(NVAR*(jvar-1)+ivar,idx)*Daux2(jvar)
          end do
          Dflux(ivar) = dscale*daux
        end do

        ! Apply limited fluxes to global vector
        Dy(i,:) = Dy(i,:)+Dflux
        Dy(j,:) = Dy(j,:)-Dflux

      end do

      !$omp end parallel do

    end subroutine doLimitADFluxes_sim

  end subroutine afcsys_buildVectorTVDBlock

  ! ****************************************************************************

!<subroutine>

  subroutine afcsys_buildVectorTVDScalar(rafcstab, rgroupFEMSet, rx, ndim,&
      fcb_calcFluxSys_sim, fcb_calcCharacteristics_sim, dscale, bclear, ry,&
      rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector for FEM-TVD schemes.
!</description>

!<input>
    ! group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! number of spatial dimensions
    integer, intent(in) :: ndim

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback function to compute local fluxes
    include 'intf_calcFluxSys_sim.inc'

    ! callback function to compute local characteristics
    include 'intf_calcCharacteristics_sim.inc'

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST) .eq. 0) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDScalar')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call gfem_getbase_DcoeffsAtEdge(rgroupFEMset, p_DcoeffsAtEdge)
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

    ! How many dimensions do we have?
    select case(ndim)
    case (NDIM1D)
      call doLimitTVD_1D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM2D)
      call doLimitTVD_2D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM3D)
      call doLimitTVD_3D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D

    subroutine doLimitTVD_1D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameter
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_1D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D

    subroutine doLimitTVD_2D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_2D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D

    subroutine doLimitTVD_3D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(3, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Z-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(3, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_3D

    !**************************************************************
    ! Assemble the upper and lower bounds Q and the sums of
    ! antidiffusive contributions P for a given set of edges

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doBoundsAndIncrements_sim(idirection, NVAR,&
        DcoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IedgeList, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dpp,Dpm,Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      integer :: idx,ivar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,i,ivar,j)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)

        ! Loop over all characteristic variables
        do ivar = 1, NVAR
          ! Set node orientation
          if (Daux1(ivar) .gt. 0) then
            i = IedgeList(2,idx)
            j = IedgeList(1,idx)
            Dflux(ivar) = -Dflux(ivar)
          else
            i = IedgeList(1,idx)
            j = IedgeList(2,idx)
          end if

          ! Assemble P's and Q's
          if (Dflux(ivar) .gt. 0) then
            Dpp(ivar,i) = Dpp(ivar,i)+Dflux(ivar)
            Dqm(ivar,i) = Dqm(ivar,i)-Dflux(ivar)
            Dqp(ivar,j) = Dqp(ivar,j)+Dflux(ivar)
          else
            Dpm(ivar,i) = Dpm(ivar,i)+Dflux(ivar)
            Dqp(ivar,i) = Dqp(ivar,i)-Dflux(ivar)
            Dqm(ivar,j) = Dqm(ivar,j)+Dflux(ivar)
          end if
        end do

      end do
      !$omp end parallel do

    end subroutine doBoundsAndIncrements_sim

    !**************************************************************
    ! Limit the antidiffusive fluxes and apply them to the vector

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doLimitADFluxes_sim(idirection, NVAR, dscale, Drp, Drm,&
        DcoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IedgeList, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      real(DP), dimension(:,:), intent(in) :: DrighteigenvectorsAtEdge,Drp,Drm
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      real(DP) :: daux
      integer :: idx,ivar,jvar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,daux,i,ivar,j,jvar)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)
        Daux2 = abs(Daux1)*DcharVariablesAtEdge(:,idx)

        ! Get position of nodes
        i = IedgeList(1,idx)
        j = IedgeList(2,idx)

        ! Loop over all characteristic variables
        ! and limit characteristic fluxes
        do ivar = 1, NVAR

          if (Daux1(ivar) .lt. 0) then
            if (Dflux(ivar) .gt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,i)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,i)*Dflux(ivar)
            end if
          else
            if (Dflux(ivar) .lt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,j)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,j)*Dflux(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        do ivar = 1, NVAR
          daux = 0.0_DP
          do jvar = 1, NVAR
            daux = daux+DrighteigenvectorsAtEdge(NVAR*(jvar-1)+ivar,idx)*Daux2(jvar)
          end do
          Dflux(ivar) = dscale*daux
        end do

        ! Apply limited fluxes to global vector
        Dy(:,i) = Dy(:,i)+Dflux
        Dy(:,j) = Dy(:,j)-Dflux

      end do
      !$omp end parallel do

    end subroutine doLimitADFluxes_sim

  end subroutine afcsys_buildVectorTVDScalar

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

end module afcstabsystemtvd

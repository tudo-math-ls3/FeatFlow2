!##############################################################################
!# ****************************************************************************
!# <name> eulerlagrange_callback3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible Euler/Navier-Stokes equations in 3D.
!#
!# The following callback functions are available:
!#
!# 1.) eulerlagrange_calcFluxGal3d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#
!# 2.) eulerlagrange_calcFluxGalNoBdr3d_sim
!#     -> Computes inviscid fluxes for standard Galerkin scheme
!#        without assembling the symmetric boundary contribution
!#
!# 3.) eulerlagrange_calcFluxScDiss3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities
!#
!# 4.) eulerlagrange_calcFluxScDissDiSp3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting scalar artificial viscosities based on
!#        dimensional splitting approach
!#
!# 5.) eulerlagrange_calcFluRoeDiss3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities
!#
!# 6.) eulerlagrange_calcFluxRoeDissDiSp3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting tensorial artificial viscosities based on
!#        dimensional splitting approach
!#
!# 7.) eulerlagrange_calcFluxRusDiss3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 8.) eulerlagrange_calcFluxRusDissDiSp3d_sim
!#     -> Computes inviscid fluxes for low-order discretisation
!#        adopting the Rusanov artificial diffusion
!#
!# 9.) eulerlagrange_calcMatDiagMatD3d_sim
!#     -> Computes local matrix for diagonal entry
!#
!# 10.) eulerlagrange_calcMatDiag3d_sim
!#      -> Computes local matrix for diagonal entry
!#
!# 11.) eulerlagrange_calcMatGalMatD3d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 12.) eulerlagrange_calcMatGal3d_sim
!#      -> Computes local matrices for standard Galerkin scheme
!#
!# 13.) eulerlagrange_calcMatScDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 14.) eulerlagrange_calcMatScDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting scalar artificial viscosities
!#
!# 15.) eulerlagrange_calcMatRoeDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 16.) eulerlagrange_calcMatRoeDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting tensorial artificial viscosities
!#
!# 17.) eulerlagrange_calcMatRusDissMatD3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov artificial viscosities
!#
!# 18.) eulerlagrange_calcMatRusDiss3d_sim
!#      -> Computes local matrices for low-order discretisation
!#         adopting the Rusanov flux artificial viscosities
!#
!# 19.) eulerlagrange_calcCharacteristics3d_sim
!#      -> Computes characteristic variables
!#
!# 20.) eulerlagrange_calcFluxFCTScalarDiss3d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting scalar artificial viscosities
!#
!# 21.) eulerlagrange_calcFluxFCTTensorDiss3d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting tensorial artificial viscosities
!#
!# 22.) eulerlagrange_calcFluxFCTRusanov3d
!#      -> Computes inviscid fluxes for FCT algorithm
!#         adopting the Rusanov artificial viscosities
!#
!# 23.) eulerlagrange_trafoFluxDensity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density
!#
!# 24.) eulerlagrange_trafoDiffDensity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density
!#
!# 25.) eulerlagrange_trafoFluxEnergy3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the energy
!#
!# 26.) eulerlagrange_trafoDiffEnergy3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the energy
!#
!# 27.) eulerlagrange_trafoFluxPressure3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the pressure
!#
!# 28.) eulerlagrange_trafoFluxVelocity3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the velocity
!#
!# 29.) eulerlagrange_trafoDiffVelocity3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the velocity
!#
!# 30.) eulerlagrange_trafoFluxMomentum3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the momentum
!#
!# 31.) eulerlagrange_trafoDiffMomentum3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the momentum
!#
!# 32.) eulerlagrange_trafoFluxDenEng3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and energy
!#
!# 33.) eulerlagrange_trafoDiffDenEng3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and energy
!#
!# 34.) eulerlagrange_trafoFluxDenPre3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density and the pessure
!#
!# 35.) eulerlagrange_trafoDiffDenPre3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density and the pessure
!#
!# 36.) eulerlagrange_trafoFluxDenPreVel3d_sim
!#      -> Computes the transformation from conservative fluxes
!#         to fluxes for the density, the pressure and the velocity
!#
!# 37.) eulerlagrange_trafoDiffDenPreVel3d_sim
!#      -> Computes the transformation from conservative solution
!#         differences to differences for the density, the pressure 
!#         and the velocity
!#
!# 38.) eulerlagrange_calcBoundaryvalues3d
!#      -> Computes the boundary values for a given node
!#
!# 39.) eulerlagrange_hadaptCallbackScalar3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in interleave format
!#
!# 40.) eulerlagrange_hadaptCallbackBlock3d
!#      -> Performs application specific tasks in the adaptation
!#         algorithm in 3D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module eulerlagrange_callback3d

  use boundaryfilter
  use collection
  use eulerlagrange_basic
  use flagship_callback
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use storage
  use thermodynamics

  implicit none

  private
  public :: eulerlagrange_calcFluxGal3d_sim
  public :: eulerlagrange_calcFluxGalNoBdr3d_sim
  public :: eulerlagrange_calcFluxScDiss3d_sim
  public :: eulerlagrange_calcFluxScDissDiSp3d_sim
  public :: eulerlagrange_calcFluxRoeDiss3d_sim
  public :: eulerlagrange_calcFluxRoeDissDiSp3d_sim
  public :: eulerlagrange_calcFluxRusDiss3d_sim
  public :: eulerlagrange_calcFluxRusDissDiSp3d_sim
  public :: eulerlagrange_calcMatDiagMatD3d_sim
  public :: eulerlagrange_calcMatDiag3d_sim
  public :: eulerlagrange_calcMatGalMatD3d_sim
  public :: eulerlagrange_calcMatGal3d_sim
  public :: eulerlagrange_calcMatScDissMatD3d_sim
  public :: eulerlagrange_calcMatScDiss3d_sim
  public :: eulerlagrange_calcMatRoeDissMatD3d_sim
  public :: eulerlagrange_calcMatRoeDiss3d_sim
  public :: eulerlagrange_calcMatRusDissMatD3d_sim
  public :: eulerlagrange_calcMatRusDiss3d_sim
  public :: eulerlagrange_calcCharacteristics3d_sim
  public :: eulerlagrange_calcFluxFCTScalarDiss3d
  public :: eulerlagrange_calcFluxFCTTensorDiss3d
  public :: eulerlagrange_calcFluxFCTRusanov3d
  public :: eulerlagrange_trafoFluxDensity3d_sim
  public :: eulerlagrange_trafoFluxEnergy3d_sim
  public :: eulerlagrange_trafoFluxPressure3d_sim
  public :: eulerlagrange_trafoFluxVelocity3d_sim
  public :: eulerlagrange_trafoFluxMomentum3d_sim
  public :: eulerlagrange_trafoFluxDenEng3d_sim
  public :: eulerlagrange_trafoFluxDenPre3d_sim
  public :: eulerlagrange_trafoFluxDenPreVel3d_sim
  public :: eulerlagrange_trafoDiffDensity3d_sim
  public :: eulerlagrange_trafoDiffEnergy3d_sim
  public :: eulerlagrange_trafoDiffPressure3d_sim
  public :: eulerlagrange_trafoDiffVelocity3d_sim
  public :: eulerlagrange_trafoDiffMomentum3d_sim
  public :: eulerlagrange_trafoDiffDenEng3d_sim
  public :: eulerlagrange_trafoDiffDenPre3d_sim
  public :: eulerlagrange_trafoDiffDenPreVel3d_sim
  public :: eulerlagrange_calcBoundaryvalues3d
  public :: eulerlagrange_hadaptCallbackScalar3d
  public :: eulerlagrange_hadaptCallbackBlock3d

contains

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxGal3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the standard
    ! Galerkin discretisation in 3D.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxGalNoBdr3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the TVD
    ! discretisation in 3D. The symmetric boundary contributions
    ! are neglected and incorporated in the antidiffusive flux.
    ! Hence, this is simply the standard Galerkin flux for the
    ! skew-symmetric internal contributions.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxGalNoBdr3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxScDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using scalar dissipation.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxScDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using scalar dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxScDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRoeDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using tensorial dissipation.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRoeDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)
    

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using tensorial dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxRoeDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRusDiss3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using the Rusanov dissipation.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxRusDissDiSp3d_sim(DdataAtEdge, DmatrixCoeffsAtEdge,&
      IdofsAtEdge, dscale, DfluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the inviscid fluxes for the
    ! low-order scheme in 3D using the Rusanov dissipation,
    ! whereby dimensional splitting is employed.
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  !   DIMENSION(ndim,2,nedges)
  ! with ndim the number of spatial dimensions
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  !   DIMENSION(4,nedges)
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
  ! OPTIONAL: collection structure
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Internodal fluxes for all edges under consideration
  !   DIMENSION(nvar,2,nedges)
  ! with nvar the number of variables at each endpoint
  real(DP), dimension(:,:,:), intent(out) :: DfluxesAtEdge
!</output>

!</subroutine>

  end subroutine eulerlagrange_calcFluxRusDissDiSp3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatDiagMatD3d_sim(DdataAtNode, DmatrixCoeffsAtNode,&
      IdofsAtNode, dscale, DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! for the diagonal block of the global operator in 3D
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtNode

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all nodes under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatDiagMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatDiag3d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices for the diagonal
    ! block of the global operator in 3D
!</description>

!<input>
  ! Nodal solution values for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DdataAtNode

  ! Entries of the coefficient matrices for all nodes under consideration
  real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode

  ! Numbers of vertices and matrix entries for all nodes under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtNode

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all nodes under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatDiag3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatGalMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatGalMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatGal3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatGal3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatScDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)
    
!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies scalar artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatScDissMatD3d_sim

!*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatScDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies scalar artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatScDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRoeDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatRoeDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRoeDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies tensorial artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatRoeDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRusDissMatD3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the diagonal of the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatRusDissMatD3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMatRusDiss3d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IdofsAtEdge, dscale,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the Galerkin matrices
    ! and applies the Rusanov artificial viscosities in 3D
!</description>

!<input>
  ! Nodal solution values for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge

  ! Entries of the coefficient matrices for all edges under consideration
  real(DP), dimension(:,:,:), intent(in) ::  DmatrixCoeffsAtEdge

  ! Numbers of vertices and matrix entries for all edges under consideration
  integer, dimension(:,:), intent(in) :: IdofsAtEdge

  ! Scaling parameter
  real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Coefficients of the matrix for all edges under consideration
  real(DP), dimension(:,:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcMatRusDiss3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcCharacteristics3d_sim(Dweight, DdataAtEdge,&
      DcharVariablesAtEdge, DeigenvaluesAtEdge,&
      DrightEigenvectorsAtEdge, DleftEigenvectorsAtEdge, rcollection)

!<description>
    ! This subroutine computes the characteristic variables in 3D
!</description>

!<input>
    ! Weighting coefficient for wave-decomposition
    real(DP), dimension(:), intent(in)  :: Dweight
    
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: Characteristic variables for all edges under consideration
    !   DIMENSION(nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DcharVariablesAtEdge
    
    ! OPTIONAL: Eigenvalues for all edges under consideration
    !   DIMENSION(nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DeigenvaluesAtEdge
    
    ! OPTIONAL: Matrices of left eigenvectors for all edges under consideration
    !   DIMENSION(nvar*nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DleftEigenvectorsAtEdge
    
    ! OPTIONAL: Matrices of right eigenvectors for all edges under consideration
    !   DIMENSION(nvar*nvar,nedge)
    ! with nvar the number of variables at each edge
    real(DP), dimension(:,:), intent(out), optional :: DrightEigenvectorsAtEdge
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcCharacteristics3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTScalarDiss3d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 3D using scalar dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcFluxFCTScalarDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTTensorDiss3d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 3D using tensorial dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcFluxFCTTensorDiss3d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcFluxFCTRusanov3d(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive fluxes for
    ! FCT algorithms in 3D using the Rusanov dissipation.
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretisation
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>

  end subroutine eulerlagrange_calcFluxFCTRusanov3d

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDensity3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)
    end do

  end subroutine eulerlagrange_trafoFluxDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDensity3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)
    end do

  end subroutine eulerlagrange_trafoDiffDensity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxEnergy3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(5,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(5,idx)
    end do

  end subroutine eulerlagrange_trafoFluxEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffEnergy3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed total density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(5,2,idx)-DdataAtEdge(5,1,idx)
    end do

  end subroutine eulerlagrange_trafoDiffEnergy3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxPressure3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the pressure in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      wi = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      wj = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)

      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(1,1,idx) = G1*(0.5_DP*(ui*ui+vi*vi+wi*wi)*DfluxesAtEdge(1,idx)-&
                                              ui*DfluxesAtEdge(2,idx)-&
                                              vi*DfluxesAtEdge(3,idx)-&
                                              wi*DfluxesAtEdge(4,idx)+DfluxesAtEdge(5,idx))
      DtransformedFluxesAtEdge(1,2,idx) =-G1*(0.5_DP*(uj*uj+vj*vj+wj*wj)*DfluxesAtEdge(1,idx)-&
                                              uj*DfluxesAtEdge(2,idx)-&
                                              vj*DfluxesAtEdge(3,idx)-&
                                              wj*DfluxesAtEdge(4,idx)+DfluxesAtEdge(5,idx))
    end do

  end subroutine eulerlagrange_trafoFluxPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffPressure3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the pressure in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: pi,pj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute pressures
      pi = G1*(DdataAtEdge(5,1,idx)-0.5_DP*(&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)+&
          DdataAtEdge(3,1,idx)*DdataAtEdge(3,1,idx)+&
          DdataAtEdge(4,1,idx)*DdataAtEdge(4,1,idx))/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(5,2,idx)-0.5_DP*(&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)+&
          DdataAtEdge(3,2,idx)*DdataAtEdge(3,2,idx)+&
          DdataAtEdge(4,2,idx)*DdataAtEdge(4,2,idx))/DdataAtEdge(1,2,idx))
     
      ! Transformed pressure difference
      DtransformedDataAtEdge(1,idx) = pj-pi
    end do

  end subroutine eulerlagrange_trafoDiffPressure3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxVelocity3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-, y- and z-velocity
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      wi = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      wj = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) = (DfluxesAtEdge(2,idx)-&
                                           ui*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-(DfluxesAtEdge(2,idx)-&
                                           uj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) = (DfluxesAtEdge(3,idx)-&
                                           vi*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-(DfluxesAtEdge(3,idx)-&
                                           vj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) = (DfluxesAtEdge(4,idx)-&
                                           wi*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =-(DfluxesAtEdge(4,idx)-&
                                           wj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)
    end do
    
  end subroutine eulerlagrange_trafoFluxVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffVelocity3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-, y- and z-velocity
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)

      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(3,idx) = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
    end do
    
  end subroutine eulerlagrange_trafoDiffVelocity3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxMomentum3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the x-, y- and z-momentum
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed momentum fluxes in x-direction
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(2,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(2,idx)

      ! Transformed momentum fluxes in y-direction
      DtransformedFluxesAtEdge(2,1,idx) = DfluxesAtEdge(3,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-DfluxesAtEdge(3,idx)

      ! Transformed momentum fluxes in z-direction
      DtransformedFluxesAtEdge(3,1,idx) = DfluxesAtEdge(4,idx)
      DtransformedFluxesAtEdge(3,2,idx) =-DfluxesAtEdge(4,idx)
    end do
    
  end subroutine eulerlagrange_trafoFluxMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffMomentum3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the x-, y- and z-momentum
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

     ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed momentum difference in x-direction
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(2,2,idx)-DdataAtEdge(2,1,idx)

      ! Transformed momentum difference in y-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(3,2,idx)-DdataAtEdge(3,1,idx)

      ! Transformed momentum difference in z-direction
      DtransformedDataAtEdge(3,idx) = DdataAtEdge(4,2,idx)-DdataAtEdge(4,1,idx)
    end do
    
  end subroutine eulerlagrange_trafoDiffMomentum3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenEng3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed total energy fluxes
      DtransformedFluxesAtEdge(2,1,idx) = DfluxesAtEdge(5,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-DfluxesAtEdge(5,idx)
    end do

  end subroutine eulerlagrange_trafoFluxDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenEng3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: idx
    
    do idx = 1, size(DdataAtEdge,3)

      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)

      ! Transformed total energy difference
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(5,2,idx)-DdataAtEdge(5,1,idx)
    end do

  end subroutine eulerlagrange_trafoDiffDenEng3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenPre3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to fluxes for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      wi = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      wj = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)

      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(2,1,idx) = G1*(0.5_DP*(ui*ui+vi*vi+wi*wi)*DfluxesAtEdge(1,idx)-&
                                          ui*DfluxesAtEdge(2,idx)-&
                                          vi*DfluxesAtEdge(3,idx)-&
                                          wi*DfluxesAtEdge(4,idx)+DfluxesAtEdge(5,idx))
      DtransformedFluxesAtEdge(2,2,idx) =-G1*(0.5_DP*(uj*uj+vj*vj+wj*wj)*DfluxesAtEdge(1,idx)-&
                                          uj*DfluxesAtEdge(2,idx)-&
                                          vj*DfluxesAtEdge(3,idx)-&
                                          wj*DfluxesAtEdge(4,idx)+DfluxesAtEdge(5,idx))
    end do

  end subroutine eulerlagrange_trafoFluxDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenPre3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density and energy in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: pi,pj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)

      ! Compute pressures
      pi = G1*(DdataAtEdge(5,1,idx)-0.5_DP*(&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)+&
          DdataAtEdge(3,1,idx)*DdataAtEdge(3,1,idx)+&
          DdataAtEdge(4,1,idx)*DdataAtEdge(4,1,idx))/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(5,2,idx)-0.5_DP*(&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)+&
          DdataAtEdge(3,2,idx)*DdataAtEdge(3,2,idx)+&
          DdataAtEdge(4,2,idx)*DdataAtEdge(4,2,idx))/DdataAtEdge(1,2,idx))
      
      ! Transformed density difference
      DtransformedDataAtEdge(1,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)
      
      ! Transformed pressure difference
      DtransformedDataAtEdge(2,idx) = pj-pi
    end do
    
  end subroutine eulerlagrange_trafoDiffDenPre3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoFluxDenPreVel3d_sim(DdataAtEdge,&
      DfluxesAtEdge, DtransformedFluxesAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation
    ! of the given flux into primitive variables in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
    
    ! Internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(in) :: DfluxesAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Transformed internodal fluxes for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(out) :: DtransformedFluxesAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: ui,uj,vi,vj,wi,wj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)
      
      ! Compute velocities
      ui = DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      vi = DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)
      wi = DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)
      uj = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)
      vj = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)
      wj = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)
      
      ! Transformed density fluxes
      DtransformedFluxesAtEdge(1,1,idx) = DfluxesAtEdge(1,idx)
      DtransformedFluxesAtEdge(1,2,idx) =-DfluxesAtEdge(1,idx)

      ! Transformed velocity fluxes in x-direction
      DtransformedFluxesAtEdge(2,1,idx) = (DfluxesAtEdge(2,idx)-&
                                           ui*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(2,2,idx) =-(DfluxesAtEdge(2,idx)-&
                                           uj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in y-direction
      DtransformedFluxesAtEdge(3,1,idx) = (DfluxesAtEdge(3,idx)-&
                                           vi*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(3,2,idx) =-(DfluxesAtEdge(3,idx)-&
                                           vj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed velocity fluxes in z-direction
      DtransformedFluxesAtEdge(4,1,idx) = (DfluxesAtEdge(4,idx)-&
                                           wi*DfluxesAtEdge(1,idx))/DdataAtEdge(1,1,idx)
      DtransformedFluxesAtEdge(4,2,idx) =-(DfluxesAtEdge(4,idx)-&
                                           wj*DfluxesAtEdge(1,idx))/DdataAtEdge(1,2,idx)

      ! Transformed pressure fluxes
      DtransformedFluxesAtEdge(5,1,idx) = G1*(0.5_DP*(ui*ui+vi*vi+wi*wi)*DfluxesAtEdge(1,idx)-&
                                              ui*DfluxesAtEdge(2,idx)-&
                                              vi*DfluxesAtEdge(3,idx)-&
                                              wi*DfluxesAtEdge(4,idx)+DfluxesAtEdge(5,idx))
      DtransformedFluxesAtEdge(5,2,idx) =-G1*(0.5_DP*(uj*uj+vj*vj+wj*wj)*DfluxesAtEdge(1,idx)-&
                                              uj*DfluxesAtEdge(2,idx)-&
                                              vj*DfluxesAtEdge(3,idx)-&
                                              wj*DfluxesAtEdge(4,idx)+DfluxesAtEdge(5,idx))
    end do

  end subroutine eulerlagrange_trafoFluxDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_trafoDiffDenPreVel3d_sim(DdataAtEdge,&
      DtransformedDataAtEdge, rcollection)

!<description>
    ! This subroutine computes the transformation of
    ! conservative to differences for the density, pressure and velocity in 3D
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    !   DIMENSION(nvar,2,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:,:), intent(in) :: DdataAtEdge
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Difference of transformed solution values for all edges under consideration
    !   DIMENSION(nvar,nedges)
    ! with nvar the number of variables at each endpoint
    real(DP), dimension(:,:), intent(out) :: DtransformedDataAtEdge
!</output>
!</subroutine>

    ! local variables
    real(DP) :: pi,pj
    integer :: idx

    do idx = 1, size(DdataAtEdge,3)

      ! Compute pressures
      pi = G1*(DdataAtEdge(5,1,idx)-0.5_DP*(&
          DdataAtEdge(2,1,idx)*DdataAtEdge(2,1,idx)+&
          DdataAtEdge(3,1,idx)*DdataAtEdge(3,1,idx)+&
          DdataAtEdge(4,1,idx)*DdataAtEdge(4,1,idx))/DdataAtEdge(1,1,idx))
      pj = G1*(DdataAtEdge(5,2,idx)-0.5_DP*(&
          DdataAtEdge(2,2,idx)*DdataAtEdge(2,2,idx)+&
          DdataAtEdge(3,2,idx)*DdataAtEdge(3,2,idx)+&
          DdataAtEdge(4,2,idx)*DdataAtEdge(4,2,idx))/DdataAtEdge(1,2,idx))
      
      ! Transformed density difference
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(1,2,idx)-DdataAtEdge(1,1,idx)

      ! Transformed velocity difference in x-direction
      DtransformedDataAtEdge(2,idx) = DdataAtEdge(2,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(2,1,idx)/DdataAtEdge(1,1,idx)
      
      ! Transformed velocity difference in y-direction
      DtransformedDataAtEdge(3,idx) = DdataAtEdge(3,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(3,1,idx)/DdataAtEdge(1,1,idx)

      ! Transformed velocity difference in z-direction
      DtransformedDataAtEdge(4,idx) = DdataAtEdge(4,2,idx)/DdataAtEdge(1,2,idx)-&
                                      DdataAtEdge(4,1,idx)/DdataAtEdge(1,1,idx)

      ! Transformed pressure difference
      DtransformedDataAtEdge(5,idx) = pj-pi
    end do

  end subroutine eulerlagrange_trafoDiffDenPreVel3d_sim

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcBoundaryvalues3d(DbdrNormal, DpointNormal,&
      DbdrValue, ibdrCondType, Du, Du0, istatus)

!<description>
    ! This subroutine computes the boundary values for a given node in 3D
!</description>

!<input>
    ! normal vector at the boundary
    real(DP), dimension(:), intent(in) :: DbdrNormal

    ! normal vector at the point on the boundary
    real(DP), dimension(:), intent(in) :: DpointNormal

    ! evaluated boundary values
    real(DP), dimension(:), intent(in) :: DbdrValue

    ! initial solution from the previous time step
    real(DP), dimension(:), intent(in) :: Du0

    ! type of boundary condition
    integer, intent(in) :: ibdrCondType
!</input>

!<inputoutput>
    ! computed boundary values
    real(DP), dimension(:), intent(inout) :: Du

    ! OPTIONAL: status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

  end subroutine eulerlagrange_calcBoundaryvalues3d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_hadaptCallbackScalar3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion
      rsolution => rcollection%p_rvectorQuickAccess1

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'eulerlagrange_hadaptCallbackScalar3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
            0.5_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)+&
                    p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      do ivar = 1, NVAR3D
        p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
            0.25_DP*(p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(3)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(4)-1)*NVAR3D+ivar)+&
                     p_Dsolution((rcollection%IquickAccess(5)-1)*NVAR3D+ivar))
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR3D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = &
              p_Dsolution((rcollection%IquickAccess(2)-1)*NVAR3D+ivar)
        end do
      else
        do ivar = 1, NVAR3D
          p_Dsolution((rcollection%IquickAccess(1)-1)*NVAR3D+ivar) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)

    end select

  end subroutine eulerlagrange_hadaptCallbackScalar3d

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_hadaptCallbackBlock3d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 3D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion
      rsolution => rcollection%p_rvectorQuickAccess1

      ! Check if solution is stored in interleave format
      if (rsolution%nblocks .ne. NVAR3D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'eulerlagrange_hadaptCallbackBlock3d')
        call sys_halt()
      end if

      ! Set pointer
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. NVAR3D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            NVAR3D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      neq = rsolution%NEQ/NVAR3D
      do ivar = 1, NVAR3D
        p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) =&
            0.25_DP*(p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(3))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(4))+&
                     p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(5)) )
      end do

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolution%NEQ/NVAR3D
        do ivar = 1, NVAR3D
          p_Dsolution((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback3d(iOperation, rcollection)

    end select

  end subroutine eulerlagrange_hadaptCallbackBlock3d

end module eulerlagrange_callback3d

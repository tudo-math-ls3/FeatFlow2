!##############################################################################
!# ****************************************************************************
!# <name> transport_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 2D.
!#
!# The following general routines are available:
!#
!# 1.) transp_refFuncBdrInt2d_sim
!#     -> Callback routine for the evaluation of the boundary integral
!#        of the target functional for goal-oriented error estimation
!#
!# 2.) transp_errorBdrInt2d_sim
!#     -> Callback routine for the evaluation of the boundary integral
!#        of the error in the target functional for goal-oriented
!#        error estimation
!#
!# 3.) transp_weightFuncBdrInt2d_sim
!#     -> Callback routine for the evaluation of the weights in
!#        the boundary integral of the target functional for
!#        goal-oriented error estimation
!#
!# 4.) transp_calcBilfBdrCond2d
!#     -> Calculates the bilinear form arising from the weak
!#        imposition of boundary conditions in 2D
!#
!# 5.) transp_calcLinfBdrCond2d
!#     -> Calculates the linear form arising from the weak
!#        imposition of boundary conditions in 2D
!#
!#
!# ****************************************************************************
!#
!# The following routines for linear velocity case are available:
!#
!# 1.) transp_calcMatDiagConvP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalConvP2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwConvP2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 2D (primal formulation)
!#
!# 4.) transp_calcMatDiagConvD2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#
!# 5.) transp_calcMatGalConvD2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#
!# 6.) transp_calcMatUpwConvD2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 2D (dual formulation)
!#
!# 7.) transp_coeffVecBdrConvP2d_sim
!#     -> Calculates the coefficients required by the evaluation of the
!#        linear form in 2D (primal formulation)
!#
!# 8.) transp_coeffMatBdrConvP2d_sim
!#     -> Calculates the coefficients required by the evaluation of the
!#        bilinear form in 2D (primal formulation)
!#
!# 9.) transp_coeffVecBdrConvD2d_sim
!#      -> Calculates the coefficients required by the evaluation of the
!#         linear form in 2D (dual formulation)
!#
!# 10.) transp_coeffMatBdrConvD2d_sim
!#     -> Calculates the coefficients required by the evaluation of the
!#        bilinear form in 2D (dual formulation)
!#
!# 11.) transp_calcVecBdrConvP2d_sim
!#      -> Calculates the group finite element coefficients for the
!#         linear form in 2D (primal formulation)
!#
!# 12.) transp_calcMatBdrConvP2d_sim
!#      -> Calculates the group finite element coefficients for the
!#         bilinear form in 2D (primal formulation)
!#
!# 13.) transp_calcVecBdrConvD2d_sim
!#      -> Calculates the group finite element coefficients for the
!#         linear form in 2D (dual formulation)
!#
!# 14.) transp_calcMatBdrConvD2d_sim
!#      -> Calculates the group finite element coefficients for the
!#         bilinear form in 2D (dual formulation)
!#
!# ****************************************************************************
!#
!# The following routines for Burgers` equation in space-time are available:
!#
!# 1.) transp_calcMatDiagSTBurgP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalSTBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwSTBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 4.) transp_coeffVecBdrSTBurgP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrSTBurgP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# ****************************************************************************
!#
!# The following routines for the Buckley-Leverett equation in
!# space-time are available:
!#
!# 1.) transp_calcMatDiagSTBLevP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Buckley-Leverett equation in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalSTBLevP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Buckley-Leverett equation in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwSTBLevP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Buckley-Leverett equation in 2D (primal formulation)
!#
!# 4.) transp_coeffVecBdrSTBLevP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrSTBLevP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# ****************************************************************************
!#
!# The following routines for the Burgers` equation in 2D are available:
!#
!# 1.) transp_calcMatDiagBurgP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 4.) transp_coeffVecBdrBurgP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrBurgP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# </purpose>
!##############################################################################

module transport_callback2d

#include "../../flagship.h"

!$use omp_lib
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use collection
  use cubature
  use derivatives
  use dofpreprocessing
  use domainintegration
  use element
  use feevaluation
  use fparser
  use fsystem
  use genoutput
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use paramlist
  use problem
  use scalarpde
  use spatialdiscretisation
  use storage

  ! Modules from transport model
  use transport_basic

  implicit none

  private

  ! generic routines
  public :: transp_refFuncBdrInt2d_sim
  public :: transp_errorBdrInt2d_sim
  public :: transp_weightFuncBdrInt2d_sim
  public :: transp_calcBilfBdrCond2d
  public :: transp_calcLinfBdrCond2d

  ! linear velocity in 2D - primal formulation
  public :: transp_calcMatDiagConvP2d_sim
  public :: transp_calcMatGalConvP2d_sim
  public :: transp_calcMatUpwConvP2d_sim
  public :: transp_coeffMatBdrConvP2d_sim
  public :: transp_coeffVecBdrConvP2d_sim
  public :: transp_calcMatBdrConvP2d_sim
  public :: transp_calcVecBdrConvP2d_sim

  ! linear velocity in 2D - dual formulation
  public :: transp_calcMatDiagConvD2d_sim
  public :: transp_calcMatGalConvD2d_sim
  public :: transp_calcMatUpwConvD2d_sim
  public :: transp_coeffMatBdrConvD2d_sim
  public :: transp_coeffVecBdrConvD2d_sim
  public :: transp_calcMatBdrConvD2d_sim
  public :: transp_calcVecBdrConvD2d_sim
  
  ! Burgers` equation in space-time - primal formulation
  public :: transp_calcMatDiagSTBurgP2d_sim
  public :: transp_calcMatGalSTBurgP2d_sim
  public :: transp_calcMatUpwSTBurgP2d_sim
  public :: transp_coeffMatBdrSTBurgP2d_sim
  public :: transp_coeffVecBdrSTBurgP2d_sim
  
  ! Buckley-Leverett equation in space-time - primal formulation
  public :: transp_calcMatDiagSTBLevP2d_sim
  public :: transp_calcMatGalSTBLevP2d_sim
  public :: transp_calcMatUpwSTBLevP2d_sim
  public :: transp_coeffMatBdrSTBLevP2d_sim
  public :: transp_coeffVecBdrSTBLevP2d_sim

  ! Burgers` equation in 2D - primal formulation
  public :: transp_calcMatDiagBurgP2d_sim
  public :: transp_calcMatGalBurgP2d_sim
  public :: transp_calcMatUpwBurgP2d_sim
  public :: transp_coeffMatBdrBurgP2d_sim
  public :: transp_coeffVecBdrBurgP2d_sim

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncBdrInt2d_sim(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues,&
      rcollection)

!<description>

    ! This subroutine is called during the calculation of errors with
    ! boundary integrals. It has to compute the values of the function
    !   $$ u {\bf v}\cdot{\bf n} $$
    ! where $u$ is the exact solution and ${\bf v}$ is the exact
    ! velocity vector and ${\bf n}$ denotes the outward unit normal.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   DquickAccess(1): simulation time
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients,Dnormal
    real(DP) :: dtime
    integer :: iel,ipoint,icompFunc,icompVelX,icompVelY
    integer :: npointsPerElement,nelements
    
    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)

    ! This subroutine assumes that the first quick access integer
    ! value holds the number of the reference function.  Moreover,
    ! quick access interger values 3 and 4 hold the numbers of the
    ! functions to be evaluated for the x-velocity and y-velocity.
    icompFunc = rcollection%IquickAccess(1)
    icompVelX = rcollection%IquickAccess(3)
    icompVelY = rcollection%IquickAccess(4)

    ! Get dimensions
    npointsPerElement = size(Dvalues,1)
    nelements         = size(Dvalues,2)

    ! Allocate temporal memory
    allocate(Dcoefficients(npointsPerElement,nelements,3))
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Evaluate the function parser in the cubature points on the
    ! boundary and store the resuls in Dcoefficients(:,:,1:3)
    call fparser_evalFuncBlockByNumber2(p_rfparser, icompFunc,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dcoefficients(:,:,1), (/dtime/))
    call fparser_evalFuncBlockByNumber2(p_rfparser, icompVelX,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dcoefficients(:,:,2), (/dtime/))
    call fparser_evalFuncBlockByNumber2(p_rfparser, icompVelY,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dcoefficients(:,:,3), (/dtime/))

    ! Get the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
                                  Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Compute the expression from the data stored in Dcoefficients
    !
    ! $$ u*(v x n) $$
    !
    ! in each cubature point on each elements
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement      
        Dvalues(ipoint,iel) = Dcoefficients(ipoint,iel,1) *&
                              (Dnormal(ipoint,iel,1)*Dcoefficients(ipoint,iel,2) +&
                               Dnormal(ipoint,iel,2)*Dcoefficients(ipoint,iel,3))
      end do
    end do
    
    ! Free temporal memory
    deallocate(Dcoefficients, Dnormal)

  end subroutine transp_refFuncBdrInt2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_errorBdrInt2d_sim(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues,&
      rcollection)

!<description>

    ! This subroutine is called during the calculation of errors with
    ! boundary integrals. It has to compute the values of the function
    !   $$ u {\bf v}\cdot{\bf n} - u_h {\bf v}_h\cdot{\bf n} $$
    ! where $u$ is the exact solution and $u_h$ is its FE approximation.
    ! Moreover, ${\bf v}$ and ${\bf v}_h$ are the exact and approximate
    ! velocity vectors and ${\bf n}$ denotes the outward unit normal.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   IquickAccess(1):     number of the reference function
    !   IquickAccess(3):     number of the reference x-velocity
    !   IquickAccess(4):     number of the reference y-velocity
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients,Dnormal
    real(DP) :: dtime
    integer :: iel,ipoint,icompFunc,icompVelX,icompVelY
    integer :: npointsPerElement,nelements

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the primal solution vector and the second quick access
    ! vector points to the velocity vector
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)

    ! This subroutine assumes that the first quick access integer
    ! value holds the number of the reference function.  Moreover,
    ! quick access interger values 3 and 4 hold the numbers of the
    ! functions to be evaluated for the x-velocity and y-velocity
    icompFunc = rcollection%IquickAccess(1)
    icompVelX = rcollection%IquickAccess(3)
    icompVelY = rcollection%IquickAccess(4)

    ! Get dimensions
    npointsPerElement = size(Dvalues,1)
    nelements         = size(Dvalues,2)

    ! Allocate temporal memory
    allocate(Dcoefficients(npointsPerElement,nelements,5))
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))
    
    ! Evaluate the FE function in the cubature points on the boundary
    call fevl_evaluate_sim(DER_FUNC, Dvalues,&
        p_rsolution%RvectorBlock(1), Dpoints, Ielements, DpointsRef)
    
    ! Evaluate the velocity field in the cubature points on the boundary
    ! and store the result in Dcoefficients(:,:,1:2)
    call fevl_evaluate_sim(DER_FUNC, Dcoefficients(:,:,1),&
        p_rvelocity%RvectorBlock(1), Dpoints, Ielements, DpointsRef)
    call fevl_evaluate_sim(DER_FUNC, Dcoefficients(:,:,2),&
        p_rvelocity%RvectorBlock(2), Dpoints, Ielements, DpointsRef)

    ! Evaluate the function parser in the cubature points on the
    ! boundary and store the resuls in Dcoefficients(:,:,3:5)
    call fparser_evalFuncBlockByNumber2(p_rfparser, icompFunc,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dcoefficients(:,:,3), (/dtime/))
    call fparser_evalFuncBlockByNumber2(p_rfparser, icompVelX,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dcoefficients(:,:,4), (/dtime/))
    call fparser_evalFuncBlockByNumber2(p_rfparser, icompVelY,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dcoefficients(:,:,5), (/dtime/))

    ! Get the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
                                  Dnormal(:,:,1), Dnormal(:,:,2), 1)

    ! Compute the expression from the data stored in Dcoefficients
    !
    ! $$ u*(v x n) - u_h*(v_h x n) $$
    !
    ! in each cubature point on each elements
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        Dvalues(ipoint,iel) = Dcoefficients(ipoint,iel,3) *&
                              (Dnormal(ipoint,iel,1)*Dcoefficients(ipoint,iel,4) +&
                               Dnormal(ipoint,iel,2)*Dcoefficients(ipoint,iel,5))-&
                              Dvalues(ipoint,iel) *&
                              (Dnormal(ipoint,iel,1)*Dcoefficients(ipoint,iel,1) +&
                               Dnormal(ipoint,iel,2)*Dcoefficients(ipoint,iel,2))
      end do
    end do

    ! Free temporal memory
    deallocate(Dcoefficients, Dnormal)

  end subroutine transp_errorBdrInt2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_weightFuncBdrInt2d_sim(rdiscretisation, DpointsRef,&
      Dpoints, ibct, DpointPar, Ielements, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It
    ! has to compute the values of a weighting function in a couple of
    ! points on a couple of elements. These values are multiplied by
    ! the calculated error.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to
    ! compute simultaneously for all these points.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   IquickAccess(2): number of the weighting function
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the weights in all the
    ! points specified in Dpoints, or the appropriate derivative of
    ! the function, respectively, according to cderivative.
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    integer :: icomp,npointsPerElement,nelements
    real(DP) :: dtime

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)

    ! Moreover, this subroutine assumes that the second quick access
    ! integer value holds the number of the function to be evaluated
    icomp = rcollection%IquickAccess(2)

    ! Get dimensions
    npointsPerElement = size(Dvalues,1)
    nelements         = size(Dvalues,2)
    
    ! Evaluate the function parser in the cubature points on the
    ! boundary and store the resuls in Dvalues
    call fparser_evalFuncBlockByNumber2(p_rfparser, icomp,&
        NDIM2D, npointsPerElement*nelements, Dpoints,&
        npointsPerElement*nelements, Dvalues, (/dtime/))

  end subroutine transp_weightFuncBdrInt2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcBilfBdrCond2d(rproblemLevel,&
      rboundaryCondition, rsolution, ssectionName, dtime, dscale,&
      fcoeff_buildMatrixScBdr2D_sim, bclear, rmatrix, rcollection)

!<description>
    ! This subroutine computes the bilinear form arising from the weak
    ! imposition of boundary conditions in 2D. The following types of
    ! boundary conditions are supported for this application
    !
    ! - (In-)homogeneous Neumann boundary conditions
    ! - Dirichlet boundary conditions
    ! - Robin boundary conditions
    ! - Flux boundary conditions
    ! - (Anti-)periodic boundary conditions
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! callback routine for nonconstant coefficient matrices.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientMatrixScBdr2D.inc'
!</intput>

!<inputoutput>
    ! matrix
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_fparser), pointer :: p_rfparser
    type(t_collection) :: rcollectionTmp
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_bilinearform) :: rform
    character(LEN=SYS_STRLEN) :: sdiffusionName
    real(DP), dimension(NDIM2D,NDIM2D) :: DdiffusionTensor
    real(DP), dimension(1) :: Dunity = (/1.0_DP/)
    integer, dimension(:), pointer :: p_IbdrCondCpIdx,p_IbdrCondType
    integer :: ibdc,isegment,ccubTypeBdr
    integer :: ivelocitytype,velocityfield,idiffusiontype

    ! Evaluate bilinear form for boundary integral and
    ! return if there are no weak boundary conditions
    if (.not.rboundaryCondition%bWeakBdrCond) return

    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameter values from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ccubTypeBdr', ccubTypeBdr)

    
    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Attach solution vector to first quick access vector of the
    ! temporal collection structure
    rcollectionTmp%p_rvectorQuickAccess1 => rsolution

    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        rboundaryCondition%rfparser, .true.)
    

    ! Attach velocity vector (if any) to second quick access vector of
    ! the temporal collection structure
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)
    else
      nullify(rcollectionTmp%p_rvectorQuickAccess2)
    end if


    ! Attach type of diffusion operator (if any) to temporal
    ! collection structure (only required for Dirichlet boundary
    ! conditions). Therefore, here we do some initialisation
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'idiffusiontype', idiffusiontype)

    ! Get function parser from collection
    p_rfparser => collct_getvalue_pars(rcollection,&
        'rfparser', ssectionName=ssectionName)
    
    ! What type of diffusion are we?
    select case(idiffusiontype)
    case default
      DdiffusionTensor = 0.0_DP
      
    case(DIFFUSION_ISOTROPIC)
      DdiffusionTensor = 0.0_DP

      ! Initialise isotropic diffusion tensor D=(/d,0;0,d/)
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(1,1))
      DdiffusionTensor(2,2) = DdiffusionTensor(1,1)

    case(DIFFUSION_ANISOTROPIC)
      ! Initialise anisotropic diffusion tensor D=(/d11,d12;d21,d22/)
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(1,1))
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=2)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(1,2))
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=3)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(2,1))
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=4)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(2,2))     
    end select
    
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrix)

    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)
    
    ! Loop over all boundary components
    do ibdc = 1, rboundaryCondition%iboundarycount
      
      ! Loop over all boundary segments
      do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
        
        ! Check if this segment has weak boundary conditions
        if (iand(p_IbdrCondType(isegment), BDRC_WEAK) .ne. BDRC_WEAK) cycle
        
        ! Prepare further quick access arrays of temporal collection
        ! structure with boundary component, type and maximum expressions
        rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(2) = isegment
        rcollectionTmp%IquickAccess(3) = rboundaryCondition%nmaxExpressions
        
        ! What type of boundary conditions are we?
        select case(iand(p_IbdrCondType(isegment), BDRC_TYPEMASK))

        case (BDRC_DIRICHLET, BDRC_PERIODIC, BDRC_ANTIPERIODIC)

          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%IquickAccess(4) = idiffusiontype
          rcollectionTmp%DquickAccess(3) = DdiffusionTensor(1,1)
          rcollectionTmp%DquickAccess(4) = DdiffusionTensor(1,2)
          rcollectionTmp%DquickAccess(5) = DdiffusionTensor(2,1)
          rcollectionTmp%DquickAccess(6) = DdiffusionTensor(2,2)
          
          ! Initialise the bilinear form
          rform%itermCount = 5
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%Idescriptors(1,2) = DER_DERIV2D_X
          rform%Idescriptors(2,2) = DER_FUNC
          rform%Idescriptors(1,3) = DER_DERIV2D_Y
          rform%Idescriptors(2,3) = DER_FUNC
          rform%Idescriptors(1,4) = DER_FUNC
          rform%Idescriptors(2,4) = DER_DERIV2D_X
          rform%Idescriptors(1,5) = DER_FUNC
          rform%Idescriptors(2,5) = DER_DERIV2D_Y
          
          ! We have no constant coefficients
          rform%ballCoeffConstant = .false.
          rform%BconstantCoeff    = .false.
          
          ! Create boundary region
          call bdrc_createRegion(rboundaryCondition,&
              ibdc, isegment-p_IbdrCondCpIdx(ibdc)+1,&
              rboundaryRegion)
          
          ! Assemble the bilinear form
          if (iand(p_IbdrCondType(ibdc), BDRC_LUMPED) .eq. BDRC_LUMPED) then
            ! ... using lumped assembly
            call bilf_buildMatrixScalarBdr2D(rform, ccubTypeBdr,&
                .false., rmatrix, fcoeff_buildMatrixScBdr2D_sim,&
                rboundaryRegion, rcollectionTmp, BILF_MATC_LUMPED)
          else
            ! ... using standard assembly
            call bilf_buildMatrixScalarBdr2D(rform, ccubTypeBdr,&
                .false., rmatrix, fcoeff_buildMatrixScBdr2D_sim,&
                rboundaryRegion, rcollectionTmp)
          end if

        case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN, BDRC_ROBIN, BDRC_FLUX)
          
          ! Initialise the bilinear form
          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          
          ! We have no constant coefficients
          rform%ballCoeffConstant = .false.
          rform%BconstantCoeff    = .false.
          
          ! Create boundary region
          call bdrc_createRegion(rboundaryCondition,&
              ibdc, isegment-p_IbdrCondCpIdx(ibdc)+1,&
              rboundaryRegion)
          
          ! Assemble the bilinear form
          if (iand(p_IbdrCondType(ibdc), BDRC_LUMPED) .eq. BDRC_LUMPED) then
            ! ... using lumped assembly
            call bilf_buildMatrixScalarBdr2D(rform, ccubTypeBdr,&
                .false., rmatrix, fcoeff_buildMatrixScBdr2D_sim,&
                rboundaryRegion, rcollectionTmp, BILF_MATC_LUMPED)
          else
            ! ... using standard assembly
            call bilf_buildMatrixScalarBdr2D(rform, ccubTypeBdr,&
                .false., rmatrix, fcoeff_buildMatrixScBdr2D_sim,&
                rboundaryRegion, rcollectionTmp)
          end if
          
        case default
          call output_line('Unsupported type of boundary conditions!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCond2d')
          call sys_halt()
          
        end select
        
      end do ! isegment

    end do ! ibdc
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)

  end subroutine transp_calcBilfBdrCond2d

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcLinfBdrCond2d(rproblemLevel, rboundaryCondition,&
      rsolution, ssectionName, dtime, dscale, fcoeff_buildVectorScBdr2D_sim,&
      bclear, rvector, rcollection)

!<description>
    ! This subroutine computes the linear form arising from the weak
    ! imposition of boundary conditions in 2D. The following types of
    ! boundary conditions are supported for this application
    !
    ! - Inhomogeneous Neumann boundary conditions
    ! - Dirichlet boundary conditions
    ! - Robin boundary conditions
    ! - Flux boundary conditions
    ! - (Anti-)periodic boundary conditions
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Whether to clear the vector before calculating the entries.
    ! If .FALSE., the new vector entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient vectors.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientVectorScBdr2D.inc'
!</intput>

!<inputoutput>
    ! vector where to store the linear form
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_fparser), pointer :: p_rfparser
    type(t_collection) :: rcollectionTmp
    type(t_boundaryRegion) :: rboundaryRegion,rboundaryRegionMirror,rregion
    type(t_linearForm) :: rform
    character(LEN=SYS_STRLEN) :: sdiffusionName
    real(DP), dimension(NDIM2D,NDIM2D) :: DdiffusionTensor
    real(DP), dimension(1) :: Dunity = (/1.0_DP/)
    integer, dimension(:), pointer :: p_IbdrCondCpIdx,p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic,p_IbdrCondPeriodic
    integer :: ibdc,isegment,ccubTypeBdr
    integer :: ivelocitytype,velocityfield,idiffusiontype
    
    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    if (.not.rboundaryCondition%bWeakBdrCond) return

    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameter values from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ccubTypeBdr', ccubTypeBdr)
    
    
    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Attach solution vector to first quick access vector of the
    ! temporal collection structure
    rcollectionTmp%p_rvectorQuickAccess1 => rsolution

    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        rboundaryCondition%rfparser, .true.)
    

    ! Attach velocity vector (if any) to second quick access vector of
    ! the temporal collection structure
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess2 =>&
          rproblemLevel%RvectorBlock(velocityfield)
    else
      nullify(rcollectionTmp%p_rvectorQuickAccess2)
    end if

    
    ! Attach type of diffusion operator (if any) to temporal
    ! collection structure (only required for Dirichlet boundary
    ! conditions). Therefore, here we do some initialisation
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'idiffusiontype', idiffusiontype)

    ! Get function parser from collection
    p_rfparser => collct_getvalue_pars(rcollection,&
        'rfparser', ssectionName=ssectionName)
    
    ! What type of diffusion are we?
    select case(idiffusiontype)
    case default
      DdiffusionTensor = 0.0_DP
      
    case(DIFFUSION_ISOTROPIC)
      DdiffusionTensor = 0.0_DP

      ! Initialise isotropic diffusion tensor D=(/d,0;0,d/)
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(1,1))
      DdiffusionTensor(2,2) = DdiffusionTensor(1,1)

    case(DIFFUSION_ANISOTROPIC)
      ! Initialise anisotropic diffusion tensor D=(/d11,d12;d21,d22/)
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=1)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(1,1))
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=2)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(1,2))
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=3)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(2,1))
      call parlst_getvalue_string(p_rparlist,&
          ssectionName, 'sdiffusionname', sdiffusionName, isubString=4)
      call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, DdiffusionTensor(2,2))     
    end select

    
    ! Clear vector?
    if (bclear) call lsysbl_clearVector(rvector)

    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)

    ! Set additional pointers for periodic boundary conditions
    if (rboundaryCondition%bPeriodic) then
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic,&
          p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic,&
          p_IbdrCondPeriodic)
    end if
    
    ! Loop over all boundary components
    do ibdc = 1, rboundaryCondition%iboundarycount
      
      ! Loop over all boundary segments
      do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
        
        ! Check if this segment has weak boundary conditions
        if (iand(p_IbdrCondType(isegment), BDRC_WEAK) .ne. BDRC_WEAK) cycle
        
        ! Prepare further quick access arrays of temporal collection
        ! structure with boundary component, type and maximum expressions
        rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(2) = isegment
        rcollectionTmp%IquickAccess(3) = rboundaryCondition%nmaxExpressions
        
        ! What type of boundary conditions are we?
        select case(iand(p_IbdrCondType(isegment), BDRC_TYPEMASK))
          
        case (BDRC_HOMNEUMANN)
          ! Do nothing for homogeneous Neumann boundary conditions
          ! since the boundary integral vanishes by construction
          
        case (BDRC_DIRICHLET, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          
          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%IquickAccess(4) = idiffusiontype
          rcollectionTmp%DquickAccess(3) = DdiffusionTensor(1,1)
          rcollectionTmp%DquickAccess(4) = DdiffusionTensor(1,2)
          rcollectionTmp%DquickAccess(5) = DdiffusionTensor(2,1)
          rcollectionTmp%DquickAccess(6) = DdiffusionTensor(2,2)
          
          ! Initialise the linear form
          rform%itermCount = 3
          rform%Idescriptors(1) = DER_FUNC
          rform%Idescriptors(2) = DER_DERIV2D_X
          rform%Idescriptors(3) = DER_DERIV2D_Y

          ! Create boundary segment
          call bdrc_createRegion(rboundaryCondition, ibdc,&
              isegment-p_IbdrCondCpIdx(ibdc)+1, rboundaryRegion)
          
          ! Check if special treatment of mirror boundary condition is required
          if ((iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
              (iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
          
            ! Create boundary region for mirror boundary in 01-parametrisation
            call bdrc_createRegion(rboundaryCondition, p_IbdrCompPeriodic(isegment),&
                p_IbdrCondPeriodic(isegment)-p_IbdrCondCpIdx(p_IbdrCompPeriodic(isegment))+1,&
                rboundaryRegionMirror)
            
            ! Attach boundary region to temporal collection structure
            call collct_setvalue_bdreg(rcollectionTmp, 'rboundaryRegionMirror',&
                rboundaryRegionMirror, .true.)
            
            ! In the callback-function, the minimum/maximum parameter
            ! values of the boundary region and its mirrored
            ! counterpartqq are required in length parametrisation to
            ! determine the parameter values of the mirrored cubature
            ! points. Therefore, we make a copy of both boundary
            ! regions, convert them to length parametrisation and
            ! attach the minimum/maximum parameter values to the quick
            ! access arrays of the temporal collection structure.
            rregion = rboundaryRegion
            call boundary_convertRegion(&
                rvector%RvectorBlock(1)%p_rspatialDiscr%p_rboundary,&
                rregion, BDR_PAR_LENGTH)
            
            ! Prepare quick access array of temporal collection structure
            rcollectionTmp%DquickAccess(7) = rregion%dminParam
            rcollectionTmp%DquickAccess(8) = rregion%dmaxParam
            
            rregion = rboundaryRegionMirror
            call boundary_convertRegion(&
                rvector%RvectorBlock(1)%p_rspatialDiscr%p_rboundary,&
                rregion, BDR_PAR_LENGTH)
            
            ! Prepare quick access array of temporal collection structure
            rcollectionTmp%DquickAccess(9)  = rregion%dminParam
            rcollectionTmp%DquickAccess(10) = rregion%dmaxParam
          end if

          ! Assemble the linear form
          call linf_buildVectorScalarBdr2d(rform, ccubTypeBdr,&
              .false., rvector%RvectorBlock(1), fcoeff_buildVectorScBdr2D_sim,&
              rboundaryRegion, rcollectionTmp)

        case (BDRC_INHOMNEUMANN, BDRC_ROBIN, BDRC_FLUX)

          ! Initialise the linear form
          rform%itermCount = 1
          rform%Idescriptors(1) = DER_FUNC
          
          ! Create boundary segment
          call bdrc_createRegion(rboundaryCondition, ibdc,&
              isegment-p_IbdrCondCpIdx(ibdc)+1, rboundaryRegion)
          
          ! Assemble the linear form
          call linf_buildVectorScalarBdr2d(rform, ccubTypeBdr,&
              .false., rvector%RvectorBlock(1), fcoeff_buildVectorScBdr2D_sim,&
              rboundaryRegion, rcollectionTmp)

        case default
          call output_line('Unsupported type of boundary copnditions !',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCond2d')
          call sys_halt()
          
        end select
        
      end do ! isegment

    end do ! ibdc

    ! Release temporal collection structure
    call collct_done(rcollectionTmp)

  end subroutine transp_calcLinfBdrCond2d

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvP2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge
    
    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
      DmatrixAtEdge(1,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DmatrixAtEdge(1,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalConvP2d_sim
 
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    if (dscale .gt. 0.0_DP) then

      do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(-DmatrixAtEdge(2,iedge), 0.0_DP,&
                -DmatrixAtEdge(3,iedge))
      end do

    else

do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{k_{ij},0,k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(DmatrixAtEdge(2,iedge), 0.0_DP,&
                DmatrixAtEdge(3,iedge))
      end do

    end if

  end subroutine transp_calcMatUpwConvP2d_sim
 
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvD2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (p_DvelocityX(InodeList(1,inode))*DcoeffsAtNode(1,inode)&
          +p_DvelocityY(InodeList(1,inode))*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
    
    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DmatrixAtEdge(1,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji] = -v_i*C_{ij}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DmatrixAtEdge(1,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji] = v_i*C_{ji}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
    end do
    
  end subroutine transp_calcMatGalConvD2d_sim
  
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvD2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    if (dscale .gt. 0.0_DP) then

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(-DmatrixAtEdge(2,iedge), 0.0_DP,&
                -DmatrixAtEdge(3,iedge))
      end do

    else

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
        DmatrixAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
        DmatrixAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
        DmatrixAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DcoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DcoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
        DmatrixAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DcoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DcoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{k_{ij},0,k_{ji}\}$
        DmatrixAtEdge(1,iedge) =&
            max(DmatrixAtEdge(2,iedge), 0.0_DP,&
                DmatrixAtEdge(3,iedge))
      end do

    end if

  end subroutine transp_calcMatUpwConvD2d_sim

  !***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvP2d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates. According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   IquickAccess(3):     maximum number of expressions
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for Dirichlet boundary conditions
    !   DquickAccess(3:6):   diffusion tensor D=(/d11,d12;d21,d22/)
    !   IquickAccess(4):     type of diffusion operator
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3:6):   diffusion tensor D=(/d11,d12;d21,d22/)
    !   DquickAccess(7):     minimim parameter of the boundary component
    !   DquickAccess(8):     maximum parameter of the boundary component
    !   DquickAccess(9):     minimim parameter of the mirror boundary component
    !   DquickAccess(10)     maximum parameter of the mirror boundary component
    !   IquickAccess(4):     type of diffusion operator
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:,:,:), allocatable :: Daux,Dnormal
    real(DP), dimension(:,:), allocatable :: Dvalue,DpointParMirror
    real(DP) :: ddiffusion,dgamma,dnv,dpenalty,dscale,dtime
    real(DP) :: dmaxParam,dmaxParamMirror,dminParam,dminParamMirror
    integer :: ibdrtype,iel,ipoint,isegment,nmaxExpr
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    type(t_dofSubset) :: rdofSubset
    real(DP), dimension(:,:,:), pointer :: p_DdofCoords
    real(DP), dimension(:,:), pointer :: p_DdofPosition
    real(DP), dimension(:), pointer :: p_Ddata,p_DvelocityX,p_DvelocityY
    real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
    real(DP), dimension(:), allocatable :: DlocalData    
    integer, dimension(:,:), pointer :: p_IdofsLoc
    logical, dimension(EL_MAXNDER) :: Bder
    integer :: idofe,idofGlob
#endif


#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrConvP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the maximum number of
    ! mathematical expressions
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)
    
    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrConvP2d_sim')
      

    case (BDRC_INHOMNEUMANN, BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann or Robin boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$ \int_{Gamma_N} wg_N ds $$
      ! 
      ! for the Neumann boundary condition
      !
      ! $$ D\nabla u\cdot{\bf n}=g_N $$
      !
      ! or
      !
      ! $$ \int_{Gamma_R} wg_R ds $$
      !
      ! for the Robin boundary condition
      !
      ! $$ \alpha u + (D\nabla u)\cdot {\bf n})=g_R $$

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      ! Initialise subset of degrees of freedom
      call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
      p_IdofsLoc   => rdofSubset%p_IdofsLoc
      p_DdofCoords => rdofSubset%p_DdofCoords
      
      ! Allocate temporal memory for velocity field, the coefficients
      ! at the DOFs and the basis function values
      allocate(Dvalue(rdofSubset%ndofsPerElement,nelements),&
            DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                      elem_getMaxDerivative(rdomainIntSubset%celement),&
                      npointsPerElement,nelements))

      ! Evaluate function values only
      Bder = .false.
      Bder(DER_FUNC) = .true.
      
      ! Evaluate the basis functions at the cubature points
      call elem_generic_sim2 (rdomainIntSubset%celement,&
          rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

      ! Evaluate the function parser for the Neumann/Robin values in the
      ! degrees of freedom on the boundary and store the result in Dvalue
      call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
          NDIM2D, rdofSubset%ndofsPerElement*nelements, p_DdofCoords,&
          rdofSubset%ndofsPerElement*nelements, Dvalue, (/dtime/))
      
      ! Allocate temporal memory for local data
      allocate(DlocalData(1))
      
      ! Loop over the cubature points and interpolate the flux
      ! boundary conditions from the DOFs to the cubature points,
      ! where they are needed by the linear form assembly routine
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Clear local data
          DlocalData = 0.0_DP
          
          do idofe = 1, rdofSubset%ndofsPerElement
            DlocalData(1) = DlocalData(1) + Dvalue(idofe,iel)*&
                DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
          end do
          
          ! Store boundary condition in the cubature points
          Dcoefficients(1,ipoint,iel) = dscale*DlocalData(1)
        end do
      end do

      ! Free temporal memory
      deallocate(DbasTrial,Dvalue,DlocalData)
      
      ! Release subset of degrees of freedom
      call dofprep_doneDofSet(rdofSubset)

#else

      ! Evaluate the function parser for the Neumann/Robin boundary
      ! values in the cubature points on the boundary and store the
      ! result in Dcoefficients(1,:,:).
      call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
          NDIM2D, npointsPerElement*nelements, Dpoints,&
          npointsPerElement*nelements, Dcoefficients(1,:,:), (/dtime/))

      ! Multiply by scaling coefficient
      Dcoefficients = dscale*Dcoefficients

#endif

      
    case (BDRC_DIRICHLET, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Dirichlet/periodic boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$  \int_{\Gamma_D} \epsilon wg_D ds                          $$ (1)
      ! $$ -\int_{\Gamma_D} (\gamma D\nabla w)\cdot{\bf n} g_D ds     $$ (2)
      ! $$ -\int_{\Gamma_D\cap\Gamma_-} ({\bf v}w)\cdot{\bf n} g_D ds $$ (3)
      !
      ! with parameter $\epsilon=C|D|/h$ and $\gamma \in \{-1,1\}$
      ! and primal inflow boundary part
      !
      ! $$\Gamma_- := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} < 0\} $$
      !
      ! For periodic boundary conditions, the "prescribed" Dirichlet
      ! value "g_D" is set to the solution value at the mirror boundary
      
      ! Get penalty and weighting parameters which are assumed
      ! constant for the entire boundary segment
      if ((iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
          (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+1, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dgamma)
        
        ! Get mirrored boundary region from collection structure
        p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
            'rboundaryRegionMirror', ssectionName=trim(rcollection%SquickAccess(1)))
        
        ! Get minimum/maximum parameter values from collection structure
        dminParam       = rcollection%DquickAccess(7)
        dmaxParam       = rcollection%DquickAccess(8)
        dminParamMirror = rcollection%DquickAccess(9)
        dmaxParamMirror = rcollection%DquickAccess(10)
      else
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+3, (/0.0_DP/), dgamma)
      end if

      ! Get norm of diffusion tensor
      ddiffusion = max( abs(rcollection%DquickAccess(3))&
                       +abs(rcollection%DquickAccess(4)),&
                        abs(rcollection%DquickAccess(5))&
                       +abs(rcollection%DquickAccess(6)))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      ! Initialise subset of degrees of freedom
      call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
      p_IdofsLoc     => rdofSubset%p_IdofsLoc
      p_DdofCoords   => rdofSubset%p_DdofCoords
      p_DdofPosition => rdofSubset%p_DdofPosition
      
      ! Allocate temporal memory for normal vector, velocity field,
      ! the coefficients at the DOFs and the basis function values
      allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                Daux(3,rdofSubset%ndofsPerElement,nelements),&
                Dvalue(rdofSubset%ndofsPerElement,nelements),&
              DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                        elem_getMaxDerivative(rdomainIntSubset%celement),&
                        npointsPerElement,nelements))

      ! Do we have to apply special treatment for periodic or
      ! antiperiodic boundary conditions?
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
        
        ! Allocate additional temporal memory
        allocate(DpointParMirror(rdofSubset%ndofsPerElement,nelements))

        ! Rescale parameter values DdofPosition on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dminParamMirror,dmaxParamMirror] $$

        call mprim_linearRescale(p_DdofPosition, dminParam, dmaxParam,&
            dminParamMirror, dmaxParamMirror, DpointParMirror)

        ! Evaluate the solution in the positions of the DOFs on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)

      elseif (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC) then

        ! Allocate additional temporal memory
        allocate(DpointParMirror(rdofSubset%ndofsPerElement,nelements))

        ! Rescale parameter values DdofPosition on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dmaxParamMirror,dminParamMirror] $$

        call mprim_linearRescale(p_DdofPosition, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)

        ! Evaluate the solution in the positions of the DOFs on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)

      else
        ! Evaluate the function parser for the Dirichlet values in the
        ! DOFs on the boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, rdofSubset%ndofsPerElement*nelements, p_DdofCoords,&
            rdofSubset%ndofsPerElement*nelements, Dvalue, (/dtime/))
      end if

      ! Evaluate function values only
      Bder = .false.
      Bder(DER_FUNC) = .true.
      
      ! Evaluate the basis functions at the cubature points
      call elem_generic_sim2 (rdomainIntSubset%celement,&
          rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
    
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
          Dnormal(:,:,1), Dnormal(:,:,2), 1)
     
      ! Assemble the diffusive part of the boundary integral, eq. (2)
      ! and impose penalty parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Compute the coefficient for the first term of the linear
            ! form which accounts for the penalty parameter.
            Daux(1,idofe,iel) = dscale*dpenalty*ddiffusion*Dvalue(idofe,iel)/&
                                rdomainIntSubset%p_DedgeLength(iel)

            ! Compute the coefficients for the second and third terms
            ! of the linear form
            !
            ! $$ -\gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} g_D ds $$
            Daux(2,idofe,iel) = -dscale*dgamma*Dvalue(idofe,iel)*&
                (rcollection%DquickAccess(3)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(idofe,iel,2))

            Daux(3,idofe,iel) = -dscale*dgamma*Dvalue(idofe,iel)*&
                (rcollection%DquickAccess(4)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(idofe,iel,2))
          end do
        end do

      else
        ! Clear all coefficients
        Daux = 0.0_DP
      end if


      ! Assemble the convective part of the boundary integral, eq. (3)
      if (associated(p_rvelocity)) then
        
        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
        
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)
            
            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Check if we are at the primal inflow boundary
            if (dnv .lt. -SYS_EPSREAL_DP)&
                Daux(1,idofe,iel) = Daux(1,idofe,iel)-dscale*dnv*Dvalue(idofe,iel)
          end do
        end do
        
      end if

      ! Allocate temporal memory for local data
      allocate(DlocalData(3))
      
      ! Loop over the cubature points and interpolate the flux
      ! boundary conditions from the DOFs to the cubature points,
      ! where they are needed by the linear form assembly routine
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Clear local data
          DlocalData = 0.0_DP
          
          do idofe = 1, rdofSubset%ndofsPerElement
            DlocalData = DlocalData + Daux(:,idofe,iel)*&
                DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
          end do
          
          ! Store boundary condition in the cubature points
          Dcoefficients(:,ipoint,iel) = DlocalData
        end do
      end do
      
      ! Free temporal memory
      deallocate(Dnormal,DbasTrial,Daux,DlocalData,Dvalue)

      ! Release subset of degrees of freedom
      call dofprep_doneDofSet(rdofSubset)

#else

      ! Allocate temporal memory for normal vector and velocity
      allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                Dvalue(npointsPerElement,nelements))

      ! Do we have to apply special treatment for periodic or
      ! antiperiodic boundary conditions?
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
        
        ! Allocate additional temporal memory
        allocate(DpointParMirror(npointsPerElement,nelements))
        
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dminParamMirror,dmaxParamMirror] $$

        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dminParamMirror, dmaxParamMirror, DpointParMirror)

        ! Evaluate the solution in the cubature points on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)
        
      elseif (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC) then

        ! Allocate additional temporal memory
        allocate(DpointParMirror(npointsPerElement,nelements))
        
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dmaxParamMirror,dminParamMirror] $$
        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)

        ! Evaluate the solution in the cubature points on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)

      else
        ! Evaluate the function parser for the Dirichlet values in the
        ! cubature points on the boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, npointsPerElement*nelements, Dpoints,&
            npointsPerElement*nelements, Dvalue, (/dtime/))
      end if
      
      ! Get the normal vectors in the cubature points on the boundary
      if (npointsPerElement .gt. 1) then
        call boundary_calcNormalVec2D(Dpoints, Dpoints,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
      else
        call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct,&
            DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
            BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
      end if
      
      ! Assemble the diffusive part of the boundary integral, eq. (2)
      ! and impose penalty parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the coefficient for the first term of the linear
            ! form which accounts for the penalty parameter.
            Dcoefficients(1,ipoint,iel) = dscale*dpenalty*ddiffusion*&
                Dvalue(ipoint,iel)/rdomainIntSubset%p_DedgeLength(iel)

            ! Compute the coefficients for the second and third terms
            ! of the linear form
            !
            ! $$ -\gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} g_D ds $$
            Dcoefficients(2,ipoint,iel) = -dscale*dgamma*Dvalue(ipoint,iel)*&
                (rcollection%DquickAccess(3)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(ipoint,iel,2))

            Dcoefficients(3,ipoint,iel) = -dscale*dgamma*Dvalue(ipoint,iel)*&
                (rcollection%DquickAccess(4)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(ipoint,iel,2))
          end do
        end do
        
      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if

      ! Assemble the convective part of the boundary integral, eq. (3)
      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Daux(npointsPerElement,nelements,NDIM2D))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
        
            ! Check if we are at the primal inflow boundary
            if (dnv .lt. -SYS_EPSREAL_DP)&
                Dcoefficients(1,ipoint,iel) =&
                Dcoefficients(1,ipoint,iel) - dscale*dnv*Dvalue(ipoint,iel)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Daux)
      end if
      
      ! Free temporal memory
      deallocate(Dnormal,Dvalue)

#endif

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Assemble the boundary integral term
      !
      ! $$ -\int_{\Gamma_-} wg_F ds $$
      !
      ! at the primal inflow boundary
      !
      ! $$\Gamma_- := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} < 0\} $$
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      if (associated(p_rvelocity)) then

        ! Initialise subset of degrees of freedom
        call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
        p_IdofsLoc   => rdofSubset%p_IdofsLoc
        p_DdofCoords => rdofSubset%p_DdofCoords
        
        ! Allocate temporal memory for normal vector, velocity field,
        ! the coefficients at the DOFs and the basis function values
        allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                  Dvalue(rdofSubset%ndofsPerElement,nelements),&
               DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                         elem_getMaxDerivative(rdomainIntSubset%celement),&
                         npointsPerElement,nelements))

        ! Evaluate function values only
        Bder = .false.
        Bder(DER_FUNC) = .true.
        
        ! Evaluate the basis functions at the cubature points
        call elem_generic_sim2 (rdomainIntSubset%celement,&
            rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
        
        ! Evaluate the function parser for the flux boundary values in
        ! the DOFs on the boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, rdofSubset%ndofsPerElement*nelements, p_DdofCoords,&
            rdofSubset%ndofsPerElement*nelements, Dvalue, (/dtime/))

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)

        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement

            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Check if we are at the primal inflow boundary
            if (dnv .lt. -SYS_EPSREAL_DP) then
              ! Multiply by normal velocity
              Dvalue(idofe,iel) = -dscale*dnv*Dvalue(idofe,iel)
            else
              ! Set zero coefficient
              Dvalue(idofe,iel) = 0.0_DP
            end if
          end do
        end do
          
        ! Allocate temporal memory for the local data
        allocate(DlocalData(1))

        ! Loop over the cubature points and interpolate the flux
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Clear local data
            DlocalData = 0.0_DP
            
            do idofe = 1, rdofSubset%ndofsPerElement
              DlocalData(1) = DlocalData(1) + Dvalue(idofe,iel)*&
                  DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
            end do
            
            ! Store flux boundary condition in the cubature points
            Dcoefficients(1,ipoint,iel) = DlocalData(1)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Dnormal,DbasTrial,Dvalue)
        
        ! Release subset of degrees of freedom
        call dofprep_doneDofSet(rdofSubset)

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if

#else

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                    Daux(npointsPerElement,nelements,NDIM2D),&
                  Dvalue(npointsPerElement,nelements))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Evaluate the function parser in the cubature points on the
        ! boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, npointsPerElement*nelements, Dpoints,&
            npointsPerElement*nelements, Dvalue, (/dtime/))

        ! Get the normal vectors in the cubature points on the boundary
        if (npointsPerElement .gt. 1) then
          call boundary_calcNormalVec2D(Dpoints, Dpoints,&
              Dnormal(:,:,1), Dnormal(:,:,2), 1)
        else
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct,&
              DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
              BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
        end if

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
        
            ! Check if we are at the primal inflow boundary
            if (dnv .lt. -SYS_EPSREAL_DP) then
              ! Multiply by scaling coefficient and normal velocity
              Dcoefficients(1,ipoint,iel) = -dscale*dnv*Dvalue(ipoint,iel)
            else
              ! Set zero coefficient
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do

        ! Free temporal memory
        deallocate(Daux,Dnormal,Dvalue)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#endif
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvP2d_sim')
      call sys_halt()
      
    end select
    
  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate the solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    !***************************************************************************
    
    subroutine doEvaluateAtBdr2d(iderType, n, Dvalues, rvectorScalar,&
        DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n
      real(DP), dimension(n), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdr2d
    
  end subroutine transp_coeffVecBdrConvP2d_sim

  !***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvD2d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   IquickAccess(3):     maximum number of expressions
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for Dirichlet boundary conditions
    !   DquickAccess(3:6):   diffusion tensor D=(/d11,d12;d21,d22/)
    !   IquickAccess(4):     type of diffusion operator
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3:6):   diffusion tensor D=(/d11,d12;d21,d22/)
    !   DquickAccess(7):     minimim parameter of the boundary component
    !   DquickAccess(8):     maximum parameter of the boundary component
    !   DquickAccess(9):     minimim parameter of the mirror boundary component
    !   DquickAccess(10)     maximum parameter of the mirror boundary component
    !   IquickAccess(4):     type of diffusion operator
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:,:,:), allocatable :: Daux,Dnormal
    real(DP), dimension(:,:), allocatable :: Dvalue,DpointParMirror
    real(DP) :: ddiffusion,dgamma,dnv,dpenalty,dscale,dtime
    real(DP) :: dmaxParam,dmaxParamMirror,dminParam,dminParamMirror
    integer :: ibdrtype,iel,ipoint,isegment,nmaxExpr
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    type(t_dofSubset) :: rdofSubset
    real(DP), dimension(:,:,:), pointer :: p_DdofCoords
    real(DP), dimension(:,:), pointer :: p_DdofPosition
    real(DP), dimension(:), pointer :: p_Ddata,p_DvelocityX,p_DvelocityY
    real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
    real(DP), dimension(:), allocatable :: DlocalData
    integer, dimension(:,:), pointer :: p_IdofsLoc
    logical, dimension(EL_MAXNDER) :: Bder
    integer :: idofe,idofGlob
#endif


#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrConvD2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the maximum number of
    ! mathematical expressions
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)
    
    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrConvD2d_sim')


    case (BDRC_INHOMNEUMANN, BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann or Robin boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$ -\int_{Gamma_N} wg_N ds $$
      ! 
      ! for the Neumann boundary condition
      !
      ! $$ D\nabla u\cdot{\bf n}=g_N $$
      !
      ! or
      !
      ! $$ -\int_{Gamma_R} wg_R ds $$
      !
      ! for the Robin boundary condition
      !
      ! $$ \alpha u + (D\nabla u)\cdot {\bf n})=g_R $$

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      ! Initialise subset of degrees of freedom
      call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
      p_IdofsLoc   => rdofSubset%p_IdofsLoc
      p_DdofCoords => rdofSubset%p_DdofCoords
      
      ! Allocate temporal memory for velocity field, the coefficients
      ! at the DOFs and the basis function values
      allocate(Dvalue(rdofSubset%ndofsPerElement,nelements),&
            DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                      elem_getMaxDerivative(rdomainIntSubset%celement),&
                      npointsPerElement,nelements))

      ! Evaluate function values only
      Bder = .false.
      Bder(DER_FUNC) = .true.
      
      ! Evaluate the basis functions at the cubature points
      call elem_generic_sim2 (rdomainIntSubset%celement,&
          rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

      ! Evaluate the function parser for the Neumann/Robin values in the
      ! degrees of freedom on the boundary and store the result in Dvalue
      call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
          NDIM2D, rdofSubset%ndofsPerElement*nelements, p_DdofCoords,&
          rdofSubset%ndofsPerElement*nelements, Dvalue, (/dtime/))
      
      ! Allocate temporal memory for local data
      allocate(DlocalData(1))
      
      ! Loop over the cubature points and interpolate the flux
      ! boundary conditions from the DOFs to the cubature points,
      ! where they are needed by the linear form assembly routine
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Clear local data
          DlocalData = 0.0_DP
          
          do idofe = 1, rdofSubset%ndofsPerElement
            DlocalData(1) = DlocalData(1) + Dvalue(idofe,iel)*&
                DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
          end do
          
          ! Store boundary condition in the cubature points
          Dcoefficients(1,ipoint,iel) = -dscale*DlocalData(1)
        end do
      end do

      ! Free temporal memory
      deallocate(DbasTrial,Dvalue,DlocalData)
      
      ! Release subset of degrees of freedom
      call dofprep_doneDofSet(rdofSubset)

#else

      ! Evaluate the function parser for the Neumann/Robin boundary
      ! values in the cubature points on the boundary and store the
      ! result in Dcoefficients(1,:,:).
      call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
          NDIM2D, npointsPerElement*nelements, Dpoints,&
          npointsPerElement*nelements, Dcoefficients(1,:,:), (/dtime/))

      ! Multiply by scaling coefficient
      Dcoefficients = -dscale*Dcoefficients

#endif

      
    case (BDRC_DIRICHLET, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Dirichlet/periodic boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$ -\int_{\Gamma_D} \epsilon wg_D ds                          $$ (1)
      ! $$ +\int_{\Gamma_D} (\gamma D\nabla w)\cdot{\bf n} g_D ds     $$ (2)
      ! $$ +\int_{\Gamma_D\cap\Gamma_+} ({\bf v}w)\cdot{\bf n} g_D ds $$ (3)
      !
      ! with parameter $\epsilon=C|D|/h$ and $\gamma \in \{-1,1\}$
      ! and dual inflow boundary part
      !
      ! $$\Gamma_+ := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} > 0\} $$
      !
      ! For periodic boundary conditions, the "prescribed" Dirichlet
      ! value "g_D" is set to the solution value at the mirror boundary
      
      ! Get penalty and weighting parameters which are assumed
      ! constant for the entire boundary segment
      if ((iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
          (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+1, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dgamma)
        
        ! Get mirrored boundary region from collection structure
        p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
            'rboundaryRegionMirror', ssectionName=trim(rcollection%SquickAccess(1)))
        
        ! Get minimum/maximum parameter values from collection structure
        dminParam       = rcollection%DquickAccess(7)
        dmaxParam       = rcollection%DquickAccess(8)
        dminParamMirror = rcollection%DquickAccess(9)
        dmaxParamMirror = rcollection%DquickAccess(10)
      else
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+3, (/0.0_DP/), dgamma)
      end if

      ! Get norm of diffusion tensor
      ddiffusion = max( abs(rcollection%DquickAccess(3))&
                       +abs(rcollection%DquickAccess(4)),&
                        abs(rcollection%DquickAccess(5))&
                       +abs(rcollection%DquickAccess(6)))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      ! Initialise subset of degrees of freedom
      call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
      p_IdofsLoc     => rdofSubset%p_IdofsLoc
      p_DdofCoords   => rdofSubset%p_DdofCoords
      p_DdofPosition => rdofSubset%p_DdofPosition
      
      ! Allocate temporal memory for normal vector, velocity field,
      ! the coefficients at the DOFs and the basis function values
      allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                Daux(3,rdofSubset%ndofsPerElement,nelements),&
                Dvalue(rdofSubset%ndofsPerElement,nelements),&
              DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                        elem_getMaxDerivative(rdomainIntSubset%celement),&
                        npointsPerElement,nelements))

      ! Do we have to apply special treatment for periodic or
      ! antiperiodic boundary conditions?
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
        
        ! Allocate additional temporal memory
        allocate(DpointParMirror(rdofSubset%ndofsPerElement,nelements))

        ! Rescale parameter values DdofPosition on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dminParamMirror,dmaxParamMirror] $$

        call mprim_linearRescale(p_DdofPosition, dminParam, dmaxParam,&
            dminParamMirror, dmaxParamMirror, DpointParMirror)

        ! Evaluate the solution in the positions of the DOFs on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)

      elseif (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC) then

        ! Allocate additional temporal memory
        allocate(DpointParMirror(rdofSubset%ndofsPerElement,nelements))

        ! Rescale parameter values DdofPosition on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dmaxParamMirror,dminParamMirror] $$

        call mprim_linearRescale(p_DdofPosition, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)

        ! Evaluate the solution in the positions of the DOFs on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)

      else
        ! Evaluate the function parser for the Dirichlet values in the
        ! DOFs on the boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, rdofSubset%ndofsPerElement*nelements, p_DdofCoords,&
            rdofSubset%ndofsPerElement*nelements, Dvalue, (/dtime/))
      end if

      ! Evaluate function values only
      Bder = .false.
      Bder(DER_FUNC) = .true.
      
      ! Evaluate the basis functions at the cubature points
      call elem_generic_sim2 (rdomainIntSubset%celement,&
          rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
          Dnormal(:,:,1), Dnormal(:,:,2), 1)
     
      ! Assemble the diffusive part of the boundary integral, eq. (2)
      ! and impose penalty parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Compute the coefficient for the first term of the linear
            ! form which accounts for the penalty parameter.
            Daux(1,idofe,iel) = -dscale*dpenalty*ddiffusion*&
                Dvalue(idofe,iel)/rdomainIntSubset%p_DedgeLength(iel)

            ! Compute the coefficients for the second and third terms
            ! of the linear form
            !
            ! $$ \gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} g_D ds $$
            Daux(2,idofe,iel) = dscale*dgamma*Dvalue(idofe,iel)*&
                (rcollection%DquickAccess(3)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(idofe,iel,2))

            Daux(3,idofe,iel) = dscale*dgamma*Dvalue(idofe,iel)*&
                (rcollection%DquickAccess(4)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(idofe,iel,2))
          end do
        end do

      else
        ! Clear all coefficients
        Daux = 0.0_DP
      end if


      ! Assemble the convective part of the boundary integral, eq. (3)
      if (associated(p_rvelocity)) then
        
        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
        
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)
            
            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Check if we are at the dual inflow boundary
            if (dnv .gt. SYS_EPSREAL_DP)&
                Daux(1,idofe,iel) = Daux(1,idofe,iel)+dscale*dnv*Dvalue(idofe,iel)
          end do
        end do
        
      end if

      ! Allocate temporal memory for local data
      allocate(DlocalData(3))
      
      ! Loop over the cubature points and interpolate the flux
      ! boundary conditions from the DOFs to the cubature points,
      ! where they are needed by the linear form assembly routine
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Clear local data
          DlocalData = 0.0_DP
          
          do idofe = 1, rdofSubset%ndofsPerElement
            DlocalData = DlocalData + Daux(:,idofe,iel)*&
                DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
          end do
          
          ! Store boundary condition in the cubature points
          Dcoefficients(:,ipoint,iel) = DlocalData
        end do
      end do
      
      ! Free temporal memory
      deallocate(Dnormal,DbasTrial,Daux,DlocalData,Dvalue)

      ! Release subset of degrees of freedom
      call dofprep_doneDofSet(rdofSubset)

#else

      ! Allocate temporal memory for normal vector and velocity
      allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                Dvalue(npointsPerElement,nelements))

      ! Do we have to apply special treatment for periodic or
      ! antiperiodic boundary conditions?
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
        
        ! Allocate additional temporal memory
        allocate(DpointParMirror(npointsPerElement,nelements))
        
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dminParamMirror,dmaxParamMirror] $$

        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dminParamMirror, dmaxParamMirror, DpointParMirror)

        ! Evaluate the solution in the cubature points on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)
        
      elseif (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC) then

        ! Allocate additional temporal memory
        allocate(DpointParMirror(npointsPerElement,nelements))
        
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region using linear mapping
        !
        ! $$ m : [dminParam,dmaxParam] -> [dmaxParamMirror,dminParamMirror] $$
        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)

        ! Evaluate the solution in the cubature points on the
        ! mirrored (!) boundary and store the result in Dvalue
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Dvalue, p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)
        
        ! Deallocate additional temporal memory
        deallocate(DpointParMirror)

      else
        ! Evaluate the function parser for the Dirichlet values in the
        ! cubature points on the boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, npointsPerElement*nelements, Dpoints,&
            npointsPerElement*nelements, Dvalue, (/dtime/))
      end if
      
      ! Get the normal vectors in the cubature points on the boundary
      if (npointsPerElement .gt. 1) then
        call boundary_calcNormalVec2D(Dpoints, Dpoints,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
      else
        call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct,&
            DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
            BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
      end if

      ! Assemble the diffusive part of the boundary integral, eq. (2)
      ! and impose penalty parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
                        
            ! Compute the coefficient for the first term of the linear
            ! form which accounts for the penalty parameter.
            Dcoefficients(1,ipoint,iel) = -dscale* dpenalty*ddiffusion*&
                Dvalue(ipoint,iel)/rdomainIntSubset%p_DedgeLength(iel)

            ! Compute the coefficients for the second and third terms
            ! of the linear form
            !
            ! $$ \gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} g_D ds $$
            Dcoefficients(2,ipoint,iel) = dscale*dgamma*Dvalue(ipoint,iel)*&
                (rcollection%DquickAccess(3)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(ipoint,iel,2))

            Dcoefficients(3,ipoint,iel) = dscale*dgamma*Dvalue(ipoint,iel)*&
                (rcollection%DquickAccess(4)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(ipoint,iel,2))
          end do
        end do

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if

      ! Assemble the convective part of the boundary integral, eq. (3)
      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Daux(npointsPerElement,nelements,NDIM2D))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
        
            ! Check if we are at the dual inflow boundary
            if (dnv .gt. SYS_EPSREAL_DP)&
                Dcoefficients(1,ipoint,iel) =&
                Dcoefficients(1,ipoint,iel) + dscale*dnv*Dvalue(ipoint,iel)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Daux)
      end if
      
      ! Free temporal memory
      deallocate(Dnormal,Dvalue)

#endif

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Assemble the boundary integral term
      !
      ! $$ +\int_{\Gamma_+} wg_F ds $$
      !
      ! at the dual inflow boundary
      !
      ! $$\Gamma_+ := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} > 0\} $$
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      if (associated(p_rvelocity)) then

        ! Initialise subset of degrees of freedom
        call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
        p_IdofsLoc   => rdofSubset%p_IdofsLoc
        p_DdofCoords => rdofSubset%p_DdofCoords
        
        ! Allocate temporal memory for normal vector, velocity field,
        ! the coefficients at the DOFs and the basis function values
        allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                  Dvalue(rdofSubset%ndofsPerElement,nelements),&
               DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                         elem_getMaxDerivative(rdomainIntSubset%celement),&
                         npointsPerElement,nelements))

        ! Evaluate function values only
        Bder = .false.
        Bder(DER_FUNC) = .true.
        
        ! Evaluate the basis functions at the cubature points
        call elem_generic_sim2 (rdomainIntSubset%celement,&
            rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
  
        ! Evaluate the function parser for the flux boundary values in
        ! the DOFs on the boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, rdofSubset%ndofsPerElement*nelements, p_DdofCoords,&
            rdofSubset%ndofsPerElement*nelements, Dvalue, (/dtime/))
        
        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
        
        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement

            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Check if we are at the dual inflow boundary
            if (dnv .gt. SYS_EPSREAL_DP) then
              ! Multiply by normal velocity
              Dvalue(idofe,iel) = dscale*dnv*Dvalue(idofe,iel)
            else
              ! Set zero coefficient
              Dvalue(idofe,iel) = 0.0_DP
            end if
          end do
        end do
          
        ! Allocate temporal memory for the local data
        allocate(DlocalData(1))

        ! Loop over the cubature points and interpolate the flux
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Clear local data
            DlocalData = 0.0_DP
            
            do idofe = 1, rdofSubset%ndofsPerElement
              DlocalData(1) = DlocalData(1) + Dvalue(idofe,iel)*&
                  DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
            end do
            
            ! Store flux boundary condition in the cubature points
            Dcoefficients(1,ipoint,iel) = DlocalData(1)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Dnormal,DbasTrial,Dvalue)
        
        ! Release subset of degrees of freedom
        call dofprep_doneDofSet(rdofSubset)

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if

#else

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                    Daux(npointsPerElement,nelements,NDIM2D),&
                  Dvalue(npointsPerElement,nelements))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Evaluate the function parser in the cubature points on the
        ! boundary and store the result in Dvalue
        call fparser_evalFuncBlockByNumber2(p_rfparser, nmaxExpr*(isegment-1)+1,&
            NDIM2D, npointsPerElement*nelements, Dpoints,&
            npointsPerElement*nelements, Dvalue, (/dtime/))
        
        ! Get the normal vectors in the cubature points on the boundary
        if (npointsPerElement .gt. 1) then
          call boundary_calcNormalVec2D(Dpoints, Dpoints,&
              Dnormal(:,:,1), Dnormal(:,:,2), 1)
        else
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct,&
              DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
              BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
        end if

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
        
            ! Check if we are at the dual inflow boundary
            if (dnv .gt. SYS_EPSREAL_DP) then
              ! Multiply by scaling coefficient and normal velocity
              Dcoefficients(1,ipoint,iel) = dscale*dnv*DValue(ipoint,iel)
            else
              ! Set zero coefficient
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do

        ! Free temporal memory
        deallocate(Daux,Dnormal,Dvalue)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#endif
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvD2d_sim')
      call sys_halt()
      
    end select
    
  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate the solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    !***************************************************************************
    
    subroutine doEvaluateAtBdr2d(iderType, n, Dvalues, rvectorScalar,&
        DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n
      real(DP), dimension(n), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdr2d
    
  end subroutine transp_coeffVecBdrConvD2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvP2d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to
    ! compute the coefficients in front of the terms of the bilinear
    ! form.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! bilinear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   IquickAccess(3):     maximum number of expressions
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for Dirichlet and periodic boundary conditions
    !   DquickAccess(3:6):   diffusion tensor D=(/d11,d12;d21,d22/)
    !   IquickAccess(4):     type of diffusion operator
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), allocatable :: Daux,Dnormal
    real(DP) :: dalpha,ddiffusion,dgamma,dnv,dpenalty,dscale,dtime
    integer :: ibdrtype,iel,ipoint,isegment,nmaxExpr
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    type(t_dofSubset) :: rdofSubset
    real(DP), dimension(:), pointer :: p_Ddata,p_DvelocityX,p_DvelocityY
    real(DP), dimension(:,:,:), pointer :: p_DdofCoords
    real(DP), dimension(:), allocatable :: DlocalData
    real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
    integer, dimension(:,:), pointer :: p_IdofsLoc
    logical, dimension(EL_MAXNDER) :: Bder
    integer :: idofe,idofGlob
#endif
    

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrConvP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the maximum number of
    ! mathematical expressions
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)
    
    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN, BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary or Robin boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$ -\int_{\Gamma_N} w({\bf v}u)\cdot{\bf n} ds $$
      !
      ! for all both types of Neumann boundary conditions.
      ! Additionally, assemble the boundary integral term
      !
      ! $$ -\int_{\Gamma_R}w[\alpha u+(\bfv u)\cdot\bfn] ds $$
      !
      ! for the Robin boundary condition
      !
      ! $$ \alpha u + (D\nabla u)\cdot {\bf n})=g_R $$
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      if (associated(p_rvelocity)) then
        
        ! Initialise subset of degrees of freedom
        call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
        p_IdofsLoc   => rdofSubset%p_IdofsLoc
        p_DdofCoords => rdofSubset%p_DdofCoords

        ! Allocate temporal memory for normal vector, velocity field,
        ! the coefficients at the DOFs and the basis function values
        allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                  Daux(1,rdofSubset%ndofsPerElement,nelements),&
               DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                         elem_getMaxDerivative(rdomainIntSubset%celement),&
                         npointsPerElement,nelements))

        ! Evaluate function values only
        Bder = .false.
        Bder(DER_FUNC) = .true.
        
        ! Evaluate the basis functions at the cubature points
        call elem_generic_sim2 (rdomainIntSubset%celement,&
            rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
       
        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
        
        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        ! Multiply the velocity vector with the normal in each degree
        ! of freedom to get the normal velocity in the degrees of
        ! freedom which is then interpolated to the cubature points.
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Scale normal velocity by the scaling parameter
            Daux(1,idofe,iel) = -dscale*dnv
          end do
        end do
        
        ! Allocate temporal memory for local data
        allocate(DlocalData(1))

        ! Loop over the cubature points and interpolate the Neumann
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the assembly routine
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Clear local data
            DlocalData = 0.0_DP
            
            do idofe = 1, rdofSubset%ndofsPerElement
              DlocalData(1) = DlocalData(1) + Daux(1,idofe,iel)*&
                  DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
            end do
            
            ! Store boundary condition in the cubature points
            Dcoefficients(1,ipoint,iel) = DlocalData(1)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Dnormal,DbasTrial,Daux,DlocalData)

        ! Release subset of degrees of freedom
        call dofprep_doneDofSet(rdofSubset)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#else

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                    Daux(npointsPerElement,nelements,NDIM2D))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        if (npointsPerElement .gt. 1) then
          call boundary_calcNormalVec2D(Dpoints, Dpoints,&
              Dnormal(:,:,1), Dnormal(:,:,2), 1)
        else
          call boundary_getNormalVec2D(rdiscretisationTest%p_rboundary, ibct,&
              DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
              BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
        end if

        ! Multiply the velocity vector with the normal in each
        ! cubature point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
            
            ! Scale normal velocity by scaling parameter
            Dcoefficients(1,ipoint,iel) = -dscale*dnv
          end do
        end do

        ! Free temporal memory
        deallocate(Daux,Dnormal)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#endif
      
      ! Do we have to prescribe Robin boundary conditions?
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ROBIN) then

        ! Evaluate the function parser for the Robin value alpha
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dalpha)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Also apply Robin boundary condition
            Dcoefficients(1,ipoint,iel) = Dcoefficients(1,ipoint,iel)&
                                        - dscale*dalpha
          end do
        end do
      end if


    case (BDRC_DIRICHLET, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Dirichlet/periodic boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$ -\int_{\Gamma_D} \epsilon wu ds                        $$ (1)
      ! $$ -\int_{\Gamma_D} w({\bf v}u-D\nabla u)\cdot{\bf n} ds  $$ (2)
      ! $$ +\int_{\Gamma_D} (\gamma D\nabla w)\cdot{\bf n}u ds    $$ (3)
      ! $$ +\int_{\Gamma_D\cap\Gamma_-} (\bf v}w)\cdot{\bf n}u ds $$ (4)
      !
      ! with parameter $\epsilon=C|D|/h$ and $\gamma \in \{-1,1\}$
      ! and primal inflow boundary part
      !
      ! $$\Gamma_- := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} < 0\} $$

      ! Get penalty and weighting parameters which are assumed
      ! constant for the entire boundary segment
      if ((iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
          (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+1, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dgamma)
      else
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+3, (/0.0_DP/), dgamma)
      end if

      ! Get norm of diffusion tensor
      ddiffusion = max( abs(rcollection%DquickAccess(3))&
                       +abs(rcollection%DquickAccess(4)),&
                        abs(rcollection%DquickAccess(5))&
                       +abs(rcollection%DquickAccess(6)))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      ! Initialise subset of degrees of freedom
      call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
      p_IdofsLoc   => rdofSubset%p_IdofsLoc
      p_DdofCoords => rdofSubset%p_DdofCoords
      
      ! Allocate temporal memory for normal vector, velocity field,
      ! the coefficients at the DOFs and the basis function values
      allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                Daux(5,rdofSubset%ndofsPerElement,nelements),&
             DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                       elem_getMaxDerivative(rdomainIntSubset%celement),&
                       npointsPerElement,nelements))

      ! Evaluate function values only
      Bder = .false.
      Bder(DER_FUNC) = .true.
      
      ! Evaluate the basis functions at the cubature points
      call elem_generic_sim2 (rdomainIntSubset%celement,&
          rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
          Dnormal(:,:,1), Dnormal(:,:,2), 1)

      ! Assemble the diffusive part of the boundary integral, eq. (2) and (3) 
      ! and impose penalry parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Compute the coefficient for the first term of the
            ! bilinear form which accounts for the penalty parameter.
            Daux(1,idofe,iel) = -dscale*dpenalty*ddiffusion/&
                                 rdomainIntSubset%p_DedgeLength(iel)
            
            ! Compute coefficients for the second and third terms of
            ! the bilinear form
            !
            ! $$ \int_{\Gamma_D} w(D\nabla u)\cdot{\bf n} ds $$
            Daux(2,idofe,iel) = dscale*&
                (rcollection%DquickAccess(3)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(idofe,iel,2))
            
            Daux(3,idofe,iel) = dscale*&
                (rcollection%DquickAccess(4)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(idofe,iel,2))

            ! Compute coefficients for the fourth and fifth terms of
            ! the bilinear form
            !
            ! $$ \gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} u ds $$
            Daux(4,idofe,iel) = dgamma*Daux(2,idofe,iel)
            Daux(5,idofe,iel) = dgamma*Daux(3,idofe,iel)
          end do
        end do

      else
        ! Clear all coefficients
        Daux = 0.0_DP
      end if
      
      ! Assemble the convective part of the boundary integral, eq. (2) and (4)
      if (associated(p_rvelocity)) then

        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        ! Multiply the velocity vector with the normal in each degree
        ! of freedom to get the normal velocity in the degrees of
        ! freedom which is then interpolated to the cubature points.
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Scale normal velocity by scaling parameter
            if (dnv .gt. SYS_EPSREAL_DP) then
              Daux(1,idofe,iel) = Daux(1,idofe,iel)-dscale*dnv
            end if
          end do
        end do

      end if
          
      ! Allocate temporal memory for the local data
      allocate(DlocalData(5))

      ! Loop over the cubature points and interpolate the boundary
      ! conditions from the DOFs to the cubature points, where
      ! they are needed by the linear form assembly routine
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Clear local data
          DlocalData = 0.0_DP
          
          do idofe = 1, rdofSubset%ndofsPerElement
            DlocalData = DlocalData + Daux(:,idofe,iel)*&
                DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
          end do
          
          ! Store boundary condition in the cubature points
          Dcoefficients(:,ipoint,iel) = DlocalData
        end do
      end do
      
      ! Free temporal memory
      deallocate(Dnormal,DbasTrial,Daux,DlocalData)
      
      ! Release subset of degrees of freedom
      call dofprep_doneDofSet(rdofSubset)
      
#else
      
      ! Allocate temporal memory for normal vector and velocity
      allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

      ! Get the normal vectors in the cubature points on the boundary
      if (npointsPerElement .gt. 1) then
        call boundary_calcNormalVec2D(Dpoints, Dpoints,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
      else
        call boundary_getNormalVec2D(rdiscretisationTest%p_rboundary, ibct,&
            DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
            BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
      end if
      
      ! Assemble the diffusive part of the boundary integral, eq. (2) and (3) 
      ! and impose penalry parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the coefficient for the first term of the
            ! bilinear form which accounts for the penalty parameter.
            Dcoefficients(1,ipoint,iel) = -dscale*dpenalty*ddiffusion/&
                                           rdomainIntSubset%p_DedgeLength(iel)

            ! Compute coefficients for the second and third terms of
            ! the bilinear form
            !
            ! $$ \int_{\Gamma_D} w(D\nabla u)\cdot{\bf n} ds $$
            Dcoefficients(2,ipoint,iel) = dscale*&
                (rcollection%DquickAccess(3)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(ipoint,iel,2))
            
            Dcoefficients(3,ipoint,iel) = dscale*&
                (rcollection%DquickAccess(4)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(ipoint,iel,2))

            ! Compute coefficients for the fourth and fifth terms of
            ! the bilinear form
            !
            ! $$ \gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} u ds $$
            Dcoefficients(4,ipoint,iel) = dgamma*Dcoefficients(2,ipoint,iel)
            Dcoefficients(5,ipoint,iel) = dgamma*Dcoefficients(3,ipoint,iel)
          end do
        end do

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if
      
      ! Assemble the convective part of the boundary integral, eq. (2) and (4)
      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for velocity
        allocate(Daux(npointsPerElement,nelements,NDIM2D))
        
        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Multiply the velocity vector with the normal in each
        ! cubature point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
            
            ! Scale normal velocity by scaling parameter and update
            ! boundary condition in the cubature points
            if (dnv .gt. SYS_EPSREAL_DP) then
              Dcoefficients(1,ipoint,iel) = Dcoefficients(1,ipoint,iel)-dscale*dnv
            end if
          end do
        end do

        ! Free temporal memory
        deallocate(Daux)
      end if

      ! Free temporal memory
      deallocate(Dnormal)

#endif
      
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      !
      ! Assemble the boundary integral term
      !
      ! $$ -\int_{\Gamma_+} w({\bf v}u)\cdot{\bf n} ds $$
      !
      ! at the primal outflow boundary
      !
      ! $$\Gamma_+ := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} > 0\} $$

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      if (associated(p_rvelocity)) then

        ! Initialise subset of degrees of freedom
        call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
        p_IdofsLoc   => rdofSubset%p_IdofsLoc
        p_DdofCoords => rdofSubset%p_DdofCoords

        ! Allocate temporal memory for normal vector, velocity field,
        ! the coefficients at the DOFs and the basis function values
        allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                  Daux(1,rdofSubset%ndofsPerElement,nelements),&
               DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                         elem_getMaxDerivative(rdomainIntSubset%celement),&
                         npointsPerElement,nelements))

        ! Evaluate function values only
        Bder = .false.
        Bder(DER_FUNC) = .true.
        
        ! Evaluate the basis functions at the cubature points
        call elem_generic_sim2 (rdomainIntSubset%celement,&
            rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
     
        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)

        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        ! Multiply the velocity vector with the normal in each degree
        ! of freedom to get the normal velocity in the degrees of
        ! freedom which is then interpolated to the cubature points.
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement

            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Only at the primal outflow boundary:
            ! Scale normal velocity by scaling parameter
            if (dnv .gt. SYS_EPSREAL_DP) then
              Daux(1,idofe,iel) = -dscale*dnv
            else
              Daux(1,idofe,iel) = 0.0_DP
            end if
          end do
        end do
        
        ! Allocate temporal memory for local data
        allocate(DlocalData(1))

        ! Loop over the cubature points and interpolate the Neumann
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Clear local data
            DlocalData = 0.0_DP
            
            do idofe = 1, rdofSubset%ndofsPerElement
              DlocalData(1) = DlocalData(1) + Daux(1,idofe,iel)*&
                  DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
            end do
            
            ! Store Neumann boundary condition in the cubature points
            Dcoefficients(1,ipoint,iel) = DlocalData(1)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Dnormal,DbasTrial,Daux,DlocalData)

        ! Release subset of degrees of freedom
        call dofprep_doneDofSet(rdofSubset)

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if

#else

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                    Daux(npointsPerElement,nelements,NDIM2D))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        if (npointsPerElement .gt. 1) then
          call boundary_calcNormalVec2D(Dpoints, Dpoints,&
              Dnormal(:,:,1), Dnormal(:,:,2), 1)
        else
          call boundary_getNormalVec2D(rdiscretisationTest%p_rboundary, ibct,&
              DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
              BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
        end if

        ! Multiply the velocity vector with the normal in each
        ! cubature point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
            
            ! Only at the primal outflow boundary:
            ! Scale normal velocity by scaling parameter
            if (dnv .gt. SYS_EPSREAL_DP) then
              Dcoefficients(1,ipoint,iel) = -dscale*dnv
            else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do

        ! Free temporal memory
        deallocate(Daux,Dnormal)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#endif
    

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvP2d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffMatBdrConvP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvD2d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to
    ! compute the coefficients in front of the terms of the bilinear
    ! form.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! bilinear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   IquickAccess(3):     maximum number of expressions
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for Dirichlet and periodic boundary conditions
    !   DquickAccess(3:6):   diffusion tensor D=(/d11,d12;d21,d22/)
    !   IquickAccess(4):     type of diffusion operator
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), allocatable :: Daux,Dnormal
    real(DP) :: dalpha,ddiffusion,dgamma,dnv,dpenalty,dscale,dtime
    integer :: ibdrtype,iel,ipoint,isegment,nmaxExpr
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    type(t_dofSubset) :: rdofSubset
    real(DP), dimension(:), pointer :: p_Ddata,p_DvelocityX,p_DvelocityY
    real(DP), dimension(:,:,:), pointer :: p_DdofCoords
    real(DP), dimension(:), allocatable :: DlocalData
    real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
    integer, dimension(:,:), pointer :: p_IdofsLoc
    logical, dimension(EL_MAXNDER) :: Bder
    integer :: idofe,idofGlob
#endif
    

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrConvD2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the maximum number of
    ! mathematical expressions
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    nmaxExpr = rcollection%IquickAccess(3)
    
    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN, BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary or Robin boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$ \int_{\Gamma_N} w({\bf v}u)\cdot{\bf n} ds $$
      !
      ! for all both types of Neumann boundary conditions.
      ! Additionally, assemble the boundary integral term
      !
      ! $$ \int_{\Gamma_R}w[\alpha u+(\bfv u)\cdot\bfn] ds $$
      !
      ! for the Robin boundary condition
      !
      ! $$ \alpha u + (D\nabla u)\cdot {\bf n})=g_R $$
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      if (associated(p_rvelocity)) then

        ! Initialise subset of degrees of freedom
        call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
        p_IdofsLoc   => rdofSubset%p_IdofsLoc
        p_DdofCoords => rdofSubset%p_DdofCoords

        ! Allocate temporal memory for normal vector, velocity field,
        ! the coefficients at the DOFs and the basis function values
        allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                  Daux(1,rdofSubset%ndofsPerElement,nelements),&
               DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                         elem_getMaxDerivative(rdomainIntSubset%celement),&
                         npointsPerElement,nelements))

        ! Evaluate function values only
        Bder = .false.
        Bder(DER_FUNC) = .true.
        
        ! Evaluate the basis functions at the cubature points
        call elem_generic_sim2 (rdomainIntSubset%celement,&
            rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
       
        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
        
        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        ! Multiply the velocity vector with the normal in each degree
        ! of freedom to get the normal velocity in the degrees of
        ! freedom which is then interpolated to the cubature points.
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Scale normal velocity by the scaling parameter
            Daux(1,idofe,iel) = dscale*dnv
          end do
        end do
        
        ! Allocate temporal memory for local data
        allocate(DlocalData(1))

        ! Loop over the cubature points and interpolate the Neumann
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the assembly routine
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Clear local data
            DlocalData = 0.0_DP
            
            do idofe = 1, rdofSubset%ndofsPerElement
              DlocalData(1) = DlocalData(1) + Daux(1,idofe,iel)*&
                  DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
            end do
            
            ! Store Neumann boundary condition in the cubature points
            Dcoefficients(1,ipoint,iel) = DlocalData(1)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Dnormal,DbasTrial,Daux,DlocalData)

        ! Release subset of degrees of freedom
        call dofprep_doneDofSet(rdofSubset)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#else

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                    Daux(npointsPerElement,nelements,NDIM2D))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        if (npointsPerElement .gt. 1) then
          call boundary_calcNormalVec2D(Dpoints, Dpoints,&
              Dnormal(:,:,1), Dnormal(:,:,2), 1)
        else
          call boundary_getNormalVec2D(rdiscretisationTest%p_rboundary, ibct,&
              DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
              BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
        end if

        ! Multiply the velocity vector with the normal in each
        ! cubature point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
            
            ! Scale normal velocity by scaling parameter
            Dcoefficients(1,ipoint,iel) = dscale*dnv
          end do
        end do

        ! Free temporal memory
        deallocate(Daux,Dnormal)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#endif
      
      ! Do we have to prescribe Robin boundary conditions?
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ROBIN) then
        
        ! Evaluate the function parser for the Robin value alpha
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dalpha)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Also apply Robin boundary condition
            Dcoefficients(1,ipoint,iel) = Dcoefficients(1,ipoint,iel)&
                                        + dscale*dalpha
          end do
        end do
      end if

    
    case (BDRC_DIRICHLET, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Dirichlet/periodic boundary conditions:
      !
      ! Assemble the boundary integral term
      !
      ! $$  \int_{\Gamma_D} \epsilon wu ds                        $$ (1)
      ! $$ +\int_{\Gamma_D} w({\bf v}u-D\nabla u)\cdot{\bf n} ds  $$ (2)
      ! $$ -\int_{\Gamma_D} (\gamma D\nabla w)\cdot{\bf n}u ds    $$ (3)
      ! $$ -\int_{\Gamma_D\cap\Gamma_+} (\bf v}w)\cdot{\bf n}u ds $$ (4)
      !
      ! with parameter $\epsilon=C|D|/h$ and $\gamma \in \{-1,1\}$
      ! and dual inflow boundary part
      !
      ! $$\Gamma_+ := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} > 0\} $$

      ! Get penalty and weighting parameters which are assumed
      ! constant for the entire boundary segment
      if ((iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
          (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+1, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dgamma)
      else
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+2, (/0.0_DP/), dpenalty)
        call fparser_evalFunction(p_rfparser,&
            nmaxExpr*(isegment-1)+3, (/0.0_DP/), dgamma)
      end if

      ! Get norm of diffusion tensor
      ddiffusion = max( abs(rcollection%DquickAccess(3))&
                       +abs(rcollection%DquickAccess(4)),&
                        abs(rcollection%DquickAccess(5))&
                       +abs(rcollection%DquickAccess(6)))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      ! Initialise subset of degrees of freedom
      call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
      p_IdofsLoc   => rdofSubset%p_IdofsLoc
      p_DdofCoords => rdofSubset%p_DdofCoords
      
      ! Allocate temporal memory for normal vector, velocity field,
      ! the coefficients at the DOFs and the basis function values
      allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                Daux(5,rdofSubset%ndofsPerElement,nelements),&
             DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                       elem_getMaxDerivative(rdomainIntSubset%celement),&
                       npointsPerElement,nelements))

      ! Evaluate function values only
      Bder = .false.
      Bder(DER_FUNC) = .true.
      
      ! Evaluate the basis functions at the cubature points
      call elem_generic_sim2 (rdomainIntSubset%celement,&
          rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
          Dnormal(:,:,1), Dnormal(:,:,2), 1)

      ! Assemble the diffusive part of the boundary integral, eq. (2) and (3) 
      ! and impose penalry parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Compute the coefficient for the first term of the
            ! bilinear form which accounts for the penalty parameter.
            Daux(1,idofe,iel) = dscale*dpenalty*ddiffusion/&
                                rdomainIntSubset%p_DedgeLength(iel)
            
            ! Compute coefficients for the second and third terms of
            ! the bilinear form
            !
            ! $$ -\int_{\Gamma_D} w(D\nabla u)\cdot{\bf n} ds $$
            Daux(2,idofe,iel) = -dscale*&
                (rcollection%DquickAccess(3)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(idofe,iel,2))
            
            Daux(3,idofe,iel) = -dscale*&
                (rcollection%DquickAccess(4)*Dnormal(idofe,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(idofe,iel,2))

            ! Compute coefficients for the fourth and fifth terms of
            ! the bilinear form
            !
            ! $$ -\gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} u ds $$
            Daux(4,idofe,iel) = dgamma*Daux(2,idofe,iel)
            Daux(5,idofe,iel) = dgamma*Daux(3,idofe,iel)
          end do
        end do

      else
        ! Clear all coefficients
        Daux = 0.0_DP
      end if
      
      ! Assemble the convective part of the boundary integral, eq. (2) and (4)
      if (associated(p_rvelocity)) then

        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        ! Multiply the velocity vector with the normal in each degree
        ! of freedom to get the normal velocity in the degrees of
        ! freedom which is then interpolated to the cubature points.
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            
            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Scale normal velocity by scaling parameter
            if (dnv .lt. -SYS_EPSREAL_DP) then
              Daux(1,idofe,iel) = Daux(1,idofe,iel)+dscale*dnv
            end if
          end do
        end do

      end if
          
      ! Allocate temporal memory for the local data
      allocate(DlocalData(5))

      ! Loop over the cubature points and interpolate the boundary
      ! conditions from the DOFs to the cubature points, where
      ! they are needed by the linear form assembly routine
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Clear local data
          DlocalData = 0.0_DP
          
          do idofe = 1, rdofSubset%ndofsPerElement
            DlocalData = DlocalData + Daux(:,idofe,iel)*&
                DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
          end do
          
          ! Store boundary condition in the cubature points
          Dcoefficients(:,ipoint,iel) = DlocalData
        end do
      end do
      
      ! Free temporal memory
      deallocate(Dnormal,DbasTrial,Daux,DlocalData)
      
      ! Release subset of degrees of freedom
      call dofprep_doneDofSet(rdofSubset)
      
#else
      
      ! Allocate temporal memory for normal vector and velocity
      allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

      ! Get the normal vectors in the cubature points on the boundary
      if (npointsPerElement .gt. 1) then
        call boundary_calcNormalVec2D(Dpoints, Dpoints,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)
      else
        call boundary_getNormalVec2D(rdiscretisationTest%p_rboundary, ibct,&
            DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
            BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
      end if
      
      ! Assemble the diffusive part of the boundary integral, eq. (2) and (3) 
      ! and impose penalry parameter $C|D|/h, eq. (1)
      if (ddiffusion .gt. 0.0_DP) then

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the coefficient for the first term of the
            ! bilinear form which accounts for the penalty parameter.
            Dcoefficients(1,ipoint,iel) = dscale*dpenalty*ddiffusion/&
                                          rdomainIntSubset%p_DedgeLength(iel)

            ! Compute coefficients for the second and third terms of
            ! the bilinear form
            !
            ! $$ -\int_{\Gamma_D} w(D\nabla u)\cdot{\bf n} ds $$
            Dcoefficients(2,ipoint,iel) = -dscale*&
                (rcollection%DquickAccess(3)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(5)*Dnormal(ipoint,iel,2))
            
            Dcoefficients(3,ipoint,iel) = -dscale*&
                (rcollection%DquickAccess(4)*Dnormal(ipoint,iel,1)+&
                 rcollection%DquickAccess(6)*Dnormal(ipoint,iel,2))

            ! Compute coefficients for the fourth and fifth terms of
            ! the bilinear form
            !
            ! $$ -\gamma \int_{\Gamma_D} (D\nabla w)\cdot{\bf n} u ds $$
            Dcoefficients(4,ipoint,iel) = dgamma*Dcoefficients(2,ipoint,iel)
            Dcoefficients(5,ipoint,iel) = dgamma*Dcoefficients(3,ipoint,iel)
          end do
        end do

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if
      
      ! Assemble the convective part of the boundary integral, eq. (2) and (4)
      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for velocity
        allocate(Daux(npointsPerElement,nelements,NDIM2D))
        
        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Multiply the velocity vector with the normal in each
        ! cubature point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
            
            ! Scale normal velocity by scaling parameter and update
            ! boundary condition in the cubature points
            if (dnv .lt. -SYS_EPSREAL_DP) then
              Dcoefficients(1,ipoint,iel) = Dcoefficients(1,ipoint,iel)+dscale*dnv
            end if
          end do
        end do

        ! Free temporal memory
        deallocate(Daux)
      end if

      ! Free temporal memory
      deallocate(Dnormal)

#endif
      
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      !
      ! Assemble the boundary integral term
      !
      ! $$ \int_{\Gamma_-} w({\bf v}u)\cdot{\bf n} ds $$
      !
      ! at the dual outflow boundary
      !
      ! $$\Gamma_- := \{{\bf x}\in\Gamma : {\bf v}\cdot{\bf n} < 0\} $$

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY

      if (associated(p_rvelocity)) then

        ! Initialise subset of degrees of freedom
        call dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset)
        p_IdofsLoc   => rdofSubset%p_IdofsLoc
        p_DdofCoords => rdofSubset%p_DdofCoords

        ! Allocate temporal memory for normal vector, velocity field,
        ! the coefficients at the DOFs and the basis function values
        allocate(Dnormal(rdofSubset%ndofsPerElement,nelements,NDIM2D),&
                  Daux(1,rdofSubset%ndofsPerElement,nelements),&
               DbasTrial(elem_igetNDofLoc(rdomainIntSubset%celement),&
                         elem_getMaxDerivative(rdomainIntSubset%celement),&
                         npointsPerElement,nelements))

        ! Evaluate function values only
        Bder = .false.
        Bder(DER_FUNC) = .true.
        
        ! Evaluate the basis functions at the cubature points
        call elem_generic_sim2 (rdomainIntSubset%celement,&
            rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
     
        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, p_DdofCoords,&
            Dnormal(:,:,1), Dnormal(:,:,2), 1)

        ! Set pointers
        call lsysbl_getbase_double(p_rsolution, p_Ddata)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
        call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)

        ! Multiply the velocity vector with the normal in each degree
        ! of freedom to get the normal velocity in the degrees of
        ! freedom which is then interpolated to the cubature points.
        do iel = 1, nelements
          do idofe = 1, rdofSubset%ndofsPerElement

            ! Get global DOF number
            idofGlob = IdofsTest(p_IdofsLoc(idofe,iel),iel)

            ! Compute the normal velocity
            dnv = Dnormal(idofe,iel,1)*p_DvelocityX(idofGlob)+&
                  Dnormal(idofe,iel,2)*p_DvelocityY(idofGlob)

            ! Only at the dual outflow boundary:
            ! Scale normal velocity by scaling parameter
            if (dnv .lt. -SYS_EPSREAL_DP) then
              Daux(1,idofe,iel) = dscale*dnv
            else
              Daux(1,idofe,iel) = 0.0_DP
            end if
          end do
        end do
        
        ! Allocate temporal memory for local data
        allocate(DlocalData(1))

        ! Loop over the cubature points and interpolate the Neumann
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Clear local data
            DlocalData = 0.0_DP
            
            do idofe = 1, rdofSubset%ndofsPerElement
              DlocalData(1) = DlocalData(1) + Daux(1,idofe,iel)*&
                  DbasTrial(p_IdofsLoc(idofe,iel),DER_FUNC,ipoint,iel)
            end do
            
            ! Store Neumann boundary condition in the cubature points
            Dcoefficients(1,ipoint,iel) = DlocalData(1)
          end do
        end do
        
        ! Free temporal memory
        deallocate(Dnormal,DbasTrial,Daux,DlocalData)

        ! Release subset of degrees of freedom
        call dofprep_doneDofSet(rdofSubset)

      else
        ! Clear all coefficients
        Dcoefficients = 0.0_DP
      end if

#else

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory for normal vector and velocity
        allocate(Dnormal(npointsPerElement,nelements,NDIM2D),&
                    Daux(npointsPerElement,nelements,NDIM2D))

        ! Evaluate the velocity field in the cubature points on the
        ! boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        if (npointsPerElement .gt. 1) then
          call boundary_calcNormalVec2D(Dpoints, Dpoints,&
              Dnormal(:,:,1), Dnormal(:,:,2), 1)
        else
          call boundary_getNormalVec2D(rdiscretisationTest%p_rboundary, ibct,&
              DpointPar, Dnormal(:,:,1), Dnormal(:,:,2),&
              BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
        end if

        ! Multiply the velocity vector with the normal in each
        ! cubature point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Compute the normal velocity
            dnv = Dnormal(ipoint,iel,1)*Daux(ipoint,iel,1) +&
                  Dnormal(ipoint,iel,2)*Daux(ipoint,iel,2)
            
            ! Only at the dual outflow boundary:
            ! Scale normal velocity by scaling parameter
            if (dnv .lt. -SYS_EPSREAL_DP) then
              Dcoefficients(1,ipoint,iel) = dscale*dnv
            else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do

        ! Free temporal memory
        deallocate(Daux,Dnormal)

      else
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
      end if

#endif
    

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvD2d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrConvD2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcMatBdrConvP2d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local matrix entries
    ! DmatrixAtNode for the node $i$.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    !   DIMENSION(ncoeffs,nnodes)
    ! with ncoeffs the number of matrix coefficients at the node
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP) :: dtime,dnv
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcMatBdrConvP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
    
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)

      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in node
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)

          ! Scale normal velocity by scaling parameter
          DmatrixAtNode(1,inode) = -dscale*dnv
        end do
        
      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

      ! Loop over all noodes at the boundary
      do inode =  1, nnodes

        ! Impose Dirichlet boundary conditions via penalty method
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode)) .gt. SYS_EPSREAL_DP) then
          DmatrixAtNode(1,inode) = -dscale*BDRC_DIRICHLET_PENALTY
        else
          DmatrixAtNode(1,inode) = 0.0_DP
        end if
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form

      DmatrixAtNode = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcMatBdrConvP2d_sim')


    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions
      
      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes

          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)

          ! Check if node i is at the primal outflow boundary
          if (dnv .gt. SYS_EPSREAL_DP) then
            DmatrixAtNode(1,inode) = -dscale*dnv
          else
            DmatrixAtNode(1,inode) = 0.0_DP
          end if
        end do

      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcMatBdrConvP2d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_calcMatBdrConvP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVecBdrConvP2d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DvectorAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local vector entries
    ! DvectorAtNode for the node $i$.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the vector for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(out) :: DvectorAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime,dnv,dval
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcVecBdrConvP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
        
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)
    
    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))
      
    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      
      DvectorAtNode = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcVecBdrConvP2d_sim')
      
      
    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Loop over all noodes at the boundary
      do inode =  1,  nnodes
        
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Multiply by scaling coefficient
          DvectorAtNode(inode) = dscale*dval
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ({\bf v}u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Loop over all noodes at the boundary
      do inode = 1, nnodes

        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Impose Dirichlet value via penalty method
          DvectorAtNode(inode) = dscale*dval*BDRC_DIRICHLET_PENALTY
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      if (associated(p_rvelocity)) then
        
        ! Loop over all noodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)
          
          if (abs(dnv) .gt. SYS_EPSREAL_DP) then
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at Robin boundary
            DvectorAtNode(inode) = -dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      if (associated(p_rvelocity)) then
        
        ! Loop over all nodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)

          ! Check if node i is at the primal outflow boundary
          if (dnv .lt. -SYS_EPSREAL_DP) then
            
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)

            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at primal outflow boundary
            DvectorAtNode(inode) = -dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      

    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      print *, "Periodic boundary conditions are not implemented yet!"
      stop

    end select

  end subroutine transp_calcVecBdrConvP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcMatBdrConvD2d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DmatrixAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local matrix entries
    ! DmatrixAtNode for the node $i$.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    !   DIMENSION(ncoeffs,nnodes)
    ! with ncoeffs the number of matrix coefficients at the node
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP) :: dtime,dnv
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcMatBdrConvP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
    
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)

      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in node
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)

          ! Scale normal velocity by scaling parameter
          DmatrixAtNode(1,inode) = dscale*dnv
        end do
        
      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

      ! Loop over all noodes at the boundary
      do inode =  1, nnodes

        ! Impose Dirichlet boundary conditions via penalty method
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode)) .gt. SYS_EPSREAL_DP) then
          DmatrixAtNode(1,inode) = dscale*BDRC_DIRICHLET_PENALTY
        else
          DmatrixAtNode(1,inode) = 0.0_DP
        end if
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form

      DmatrixAtNode = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcMatBdrConvP2d_sim')


    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions
      
      if (associated(p_rvelocity)) then

        ! Loop over all noodes at the boundary
        do inode =  1, nnodes

          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)

          ! Check if node i is at the dual outflow boundary
          if (dnv .lt. -SYS_EPSREAL_DP) then
            DmatrixAtNode(1,inode) = dscale*dnv
          else
            DmatrixAtNode(1,inode) = 0.0_DP
          end if
        end do

      else
        ! Clear coefficients for zero velocity
        DmatrixAtNode = 0.0_DP
      end if
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcMatBdrConvD2d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_calcMatBdrConvD2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVecBdrConvD2d_sim(DdataAtNode, DcoeffsAtNode,&
      InodeList, dscale, nnodes, DvectorAtNode, rcollection)

!<description>
    ! Given the solution data DdataAtNode and auxiliary coefficients
    ! DcoeffsAtNode this subroutine computes the local vector entries
    ! DvectorAtNode for the node $i$.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    !   DIMENSION(ndim,nnodes)
    ! with ndim the number of spatial dimensions
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of nodes and matrix entries for all nodes under consideration
    !   DIMENSION(2,nnodes)
    integer, dimension(:,:), intent(in) :: InodeList
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale
    
    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: collection structure
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: coordinates of the degrees of freedom
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the vector for all nodes under consideration
    !   DIMENSION(nnodes)
    real(DP), dimension(:), intent(out) :: DvectorAtNode
!</output>

!</subroutine>

    ! local variable
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rvelocity,p_rdofCoords
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime,dnv,dval
    integer :: inode,ibdrtype,isegment,i

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_calcVecBdrConvD2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! points to the coordinates of the degrees of freedom and to the
    ! velocity field (if any)
    p_rdofCoords => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Set pointers
    call lsyssc_getbase_double(p_rdofCoords%RvectorBlock(1), p_DdofCoords)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
        
    ! The first quick access double values hold the simulation time
    dtime  = rcollection%DquickAccess(1)
    
    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))
      
    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.

      DvectorAtNode = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcVecBdrConvP2d_sim')
      
      
    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Loop over all nodes at the boundary
      do inode =  1,  nnodes
        
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Multiply by scaling coefficient
          DvectorAtNode(inode) = dscale*dval
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ({\bf v}u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Loop over all nodes at the boundary
      do inode = 1, nnodes
        
        if (abs(DcoeffsAtNode(1,inode))+&
            abs(DcoeffsAtNode(2,inode)) .gt. SYS_EPSREAL_DP) then
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Set values for function parser
          Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
          Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Impose Dirichlet value via penalty method
          DvectorAtNode(inode) = -dscale*dval*BDRC_DIRICHLET_PENALTY
        else
          DvectorAtNode(inode) = 0.0_DP
        end if
      end do

      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ ({\bf v}u-d\nabla u)\cdot{\bf n} = ({\bf v}g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      if (associated(p_rvelocity)) then
        
        ! Loop over all nodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)
          
          if (abs(dnv) .gt. SYS_EPSREAL_DP) then
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at Robin boundary
            DvectorAtNode(inode) = dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ ({\bf v}u-d\nabla u)\cdot{\bf n} = ({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      if (associated(p_rvelocity)) then
        
        ! Loop over all nodes at the boundary
        do inode =  1,  nnodes
          
          ! Get global node number
          i = InodeList(1,inode)
          
          ! Compute normal velocity in nodes i
          dnv = DcoeffsAtNode(1,inode)*p_DvelocityX(i)+&
                DcoeffsAtNode(2,inode)*p_DvelocityY(i)

          ! Check if node i is at the dual outflow boundary
          if (dnv .gt. SYS_EPSREAL_DP) then
            
            ! Set values for function parser
            Dvalue(1) = p_DdofCoords((i-1)*NDIM2D+1)
            Dvalue(2) = p_DdofCoords((i-1)*NDIM2D+2)

            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Set value at dual outflow boundary
            DvectorAtNode(inode) = dscale*dnv*dval
          else
            DvectorAtNode(inode) = 0.0_DP
          end if
        end do
        
      else
        ! Clear coefficients for zero velocity
        DvectorAtNode = 0.0_DP
      end if
      

    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      print *, "Periodic boundary conditions are not implemented yet!"
      stop

    end select

  end subroutine transp_calcVecBdrConvD2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagSTBurgP2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for space-time formulation of the
    ! one-dimensional primal Burgers equation $du/dt+df(u)/dx=0$,
    ! whereby the flux function is given by $f(u)=0.5*u^2$.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>
    
    ! local variables
    integer :: inode

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficients $k_{ii} = [0.5*u_i,1]*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          (0.5*DdataAtNode(inode)*DcoeffsAtNode(1,inode)&
          +                       DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficients $k_{ii} = -[0.5*u_i,1]*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          (0.5*DdataAtNode(inode)*DcoeffsAtNode(1,inode)&
          +                       DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagSTBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalSTBurgP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! one-dimensional primal Burgers equation $du/dt+df(u)/dx=0$,
    ! whereby the flux function is given by $f(u)=0.5*u^2$.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
    
    ! local variables
    integer :: iedge

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = [0.5*u_j,1]*C_{ji}$
      DmatrixAtEdge(1,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DcoeffsAtEdge(1,2,iedge)&
          +                            DcoeffsAtEdge(2,2,iedge))

      ! Compute convective coefficient $k_{ji} = [0.5*u_i,1]*C_{ij}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DcoeffsAtEdge(1,1,iedge)&
          +                            DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[0.5*u_j,1]*C_{ij}$
      DmatrixAtEdge(1,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DcoeffsAtEdge(1,1,iedge)&
          +                            DcoeffsAtEdge(2,1,iedge))

      ! Compute convective coefficient $k_{ji} = -[0.5*u_i,1]*C_{ji}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DcoeffsAtEdge(1,2,iedge)&
          +                            DcoeffsAtEdge(2,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalSTBurgP2d_sim
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwSTBurgP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! one-dimensional primal Burgers equation $du/dt+df(u)/dx=0$,
    ! whereby the flux function is given by $f(u)=0.5*u^2$.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
    
    ! local variables
    integer :: iedge

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = [0.5*u_j,1]*C_{ji}$
      DmatrixAtEdge(2,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DcoeffsAtEdge(1,2,iedge)&
          +                            DcoeffsAtEdge(2,2,iedge))

      ! Compute convective coefficient $k_{ji} = [0.5*u_i,1]*C_{ij}$
      DmatrixAtEdge(3,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DcoeffsAtEdge(1,1,iedge)&
          +                            DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[0.5*u_j,1]*C_{ij}$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DcoeffsAtEdge(1,1,iedge)&
          +                            DcoeffsAtEdge(2,1,iedge))

      ! Compute convective coefficient $k_{ji} = -[0.5*u_i,1]*C_{ji}$
      DmatrixAtEdge(3,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DcoeffsAtEdge(1,2,iedge)&
          +                            DcoeffsAtEdge(2,2,iedge))
#endif

      ! Compute artificial diffusion coefficient
      !   $d_{ij} = abs(v_{ij}*0.5*(cx_{ij}-cx_{ji}) + max\{cy_{ij}, cy_{ji}\}$,
      ! where
      !   $v_{ij} = 0.5*(f(u_j)-f(u_i))/(u_j-u_i) = 0.5*(u_i+u_j)$
      DmatrixAtEdge(1,iedge) = dscale*&
          abs(0.25_DP*(DdataAtEdge(1,iedge)+DdataAtEdge(2,iedge))*&
              (DcoeffsAtEdge(1,1,iedge)-DcoeffsAtEdge(1,2,iedge)))+&
          max(-DcoeffsAtEdge(2,1,iedge),-DcoeffsAtEdge(2,2,iedge))
    end do

  end subroutine transp_calcMatUpwSTBurgP2d_sim

  !***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrSTBurgP2d_sim(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! 1D Burgers equation in space and time.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:), pointer :: DnormalX,DnormalY
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrSTBurgP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))
    
    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP

      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrSTBurgP2d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale*Dcoefficients(1,ipoint,iel)
        end do
      end do


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow [0.5*u,1]*u\cdot{\bf n}=[0.5*g,1]*g\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale*dval*BDRC_DIRICHLET_PENALTY
        end do
      end do

      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ -([0.5*u,1]*u-d\nabla u)\cdot{\bf n} = -([0.5*g,1]*g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,1))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,1))
          
          ! Compute the normal velocity and impose Dirichlet boundary condition
          dnv = DnormalX(ipoint,iel)*0.5*Daux(ipoint,iel,1) + DnormalY(ipoint,iel)
          Dcoefficients(1,ipoint,iel) = -dscale*dnv*Daux(ipoint,iel,1)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, DnormalX, DnormalY)


    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u,1]*u-d\nabla u)\cdot{\bf n} = -([0.5*g,1]*g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,2))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))

      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,:1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2). Multiply the velocity vector [0.5*u,1]
      ! with the normal in each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,2))

          ! Compute the normal velocity
          dnv = DnormalX(ipoint,iel)*0.5*Daux(ipoint,iel,1) + DnormalY(ipoint,iel)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            ! Compute the prescribed normal velocity
            dnv = DnormalX(ipoint,iel)*0.5*Daux(ipoint,iel,2) + DnormalY(ipoint,iel)
            Dcoefficients(1,ipoint,iel) = -dscale*dnv*Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, DnormalX, DnormalY)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrSTBurgP2d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffVecBdrSTBurgP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrSTBurgP2d_sim(&
      rdiscretisationTrial, rdiscretisationTest, rform, nelements,&
      npointsPerElement, Dpoints, ibct, DpointPar, IdofsTrial,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! 1D Burgers equation in space and time.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:), pointer :: DnormalX,DnormalY
    real(DP) :: dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrSTBurgP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,1))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
      
      ! Multiply the velocity vector [0.5*u,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = DnormalX(ipoint,iel)*0.5*Daux(ipoint,iel,1) + DnormalY(ipoint,iel)
          
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale*dnv
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, DnormalX, DnormalY)

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale*BDRC_DIRICHLET_PENALTY
        end do
      end do
      
      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form
      Dcoefficients = 0.0_DP
      
      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrSTBurgP2d_sim')

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,1))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))

      ! Evaluate the solution field in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)

      ! Multiply the velocity vector [0.5*u,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = DnormalX(ipoint,iel)*0.5*Daux(ipoint,iel,1) + DnormalY(ipoint,iel)
          
          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale*dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, DnormalX, DnormalY)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrSTBurgP2d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrSTBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagSTBLevP2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui
    integer :: inode
    
    do inode = 1, nnodes

      ui = DdataAtNode(inode)

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = [a_i,1]*C_{ii}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DmatrixAtNode(1,inode) = dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DcoeffsAtNode(1,inode)&
          +                                          DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = -[a_i,1]*C_{ii}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DmatrixAtNode(1,inode) = -dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DcoeffsAtNode(1,inode)&
          +                                          DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagSTBLevP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalSTBLevP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj
    integer :: iedge
    
    do iedge = 1, nedges

      ui = DdataAtEdge(1,iedge); uj = DdataAtEdge(2,iedge)

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = [a_j,1]*C_{ji}$
      ! where $a_j=u_j/(u_^2+0.5*(1-u_j)^2)$
      DmatrixAtEdge(1,iedge) = dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DcoeffsAtEdge(1,2,iedge)&
          +                                          DcoeffsAtEdge(2,2,iedge))
      
      ! Compute convective coefficient $k_{ji} = [a_i,1]*C_{ij}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DmatrixAtEdge(2,iedge) = dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DcoeffsAtEdge(1,1,iedge)&
          +                                          DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[a_j,1]*C_{ij}$
      ! where $a_j=u_j/(u_j^2+0.5*(1-u_j)^2)$
      DmatrixAtEdge(1,iedge) = -dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DcoeffsAtEdge(1,1,iedge)&
          +                                          DcoeffsAtEdge(2,1,iedge))
      
      ! Compute convective coefficient $k_{ji} = -[a_i,1]*C_{ji}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DcoeffsAtEdge(1,2,iedge)&
          +                                          DcoeffsAtEdge(2,2,iedge))
#endif
    end do
    
  end subroutine transp_calcMatGalSTBLevP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwSTBLevP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,vij
    integer :: iedge
    
    do iedge = 1, nedges

      ui = DdataAtEdge(1,iedge); uj = DdataAtEdge(2,iedge)

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = [a_j,1]*C_{ji}$
      ! where $a_j=u_j/(u_^2+0.5*(1-u_j)^2)$
      DmatrixAtEdge(2,iedge) = dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DcoeffsAtEdge(1,2,iedge)&
          +                                          DcoeffsAtEdge(2,2,iedge))
      
      ! Compute convective coefficient $k_{ji} = [a_i,1]*C_{ij}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DmatrixAtEdge(3,iedge) = dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DcoeffsAtEdge(1,1,iedge)&
          +                                          DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[a_j,1]*C_{ij}$
      ! where $a_j=u_j/(u_j^2+0.5*(1-u_j)^2)$
      DmatrixAtEdge(2,iedge) = -dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DcoeffsAtEdge(1,1,iedge)&
          +                                          DcoeffsAtEdge(2,1,iedge))
      
      ! Compute convective coefficient $k_{ji} = -[a_i,1]*C_{ji}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DmatrixAtEdge(3,iedge) = -dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DcoeffsAtEdge(1,2,iedge)&
          +                                          DcoeffsAtEdge(2,2,iedge))
#endif

      ! Calculate the characteristic speed
      if (abs(ui-uj) .gt. SYS_EPSREAL_DP) then
        vij = (uj*uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))&
              -ui*ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui)))/(uj-ui)
      else
        vij = ui*(1.0_DP-ui)/(ui*ui-0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))**2
      end if

      ! Compute artificial diffusion coefficient
      ! $d_{ij} = abs(0.5*v_{ij}*(c_{ij}-c_{ji}) + max\{cy_{ij}, cy_{ji}\}$
      DmatrixAtEdge(1,iedge) = dscale*&
          abs(0.5_DP*vij*(DcoeffsAtEdge(1,1,iedge)-&
                          DcoeffsAtEdge(1,2,iedge)))+&
          max(-DcoeffsAtEdge(2,1,iedge),-DcoeffsAtEdge(2,2,iedge))
    end do
    
  end subroutine transp_calcMatUpwSTBLevP2d_sim

  !***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrSTBLevP2d_sim(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! 1D Buckley-Leverett equation in space and time.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:), pointer :: DnormalX,DnormalY
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrSTBLevP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))
    
    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP

      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrSTBLevP2d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale*Dcoefficients(1,ipoint,iel)
        end do
      end do


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow [a,1]*u\cdot{\bf n}=[a,1]*g\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale*dval*BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ -([a,1]*u-d\nabla u)\cdot{\bf n} = -([a,1]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! and do not include any boundary integral into the bilinear form at all.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,1))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1). Multiply the velocity vector [a,1] with
      ! the normal in each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,1))
          
          ! Compute the normal velocity and impose Dirichlet boundary condition
          dnv = DnormalX(ipoint,iel)*Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + DnormalY(ipoint,iel)
          Dcoefficients(1,ipoint,iel) = -dscale*dnv*Daux(ipoint,iel,1)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, DnormalX, DnormalY)


    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([a,1]*u-d\nabla u)\cdot{\bf n} = -([a,1]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,2))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))

      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2). Multiply the velocity vector [a,1] with
      ! the normal in each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,2))
          
          ! Compute the normal velocity
          dnv = DnormalX(ipoint,iel)*Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + DnormalY(ipoint,iel)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            ! Compute the prescribed normal velocity
            dnv = DnormalX(ipoint,iel)*Daux(ipoint,iel,2)/(Daux(ipoint,iel,2)**2 +&
                  0.5*(1-Daux(ipoint,iel,2))**2) + DnormalY(ipoint,iel)
            Dcoefficients(1,ipoint,iel) = -dscale*dnv*Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, DnormalX, DnormalY)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrSTBLevP2d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_coeffVecBdrSTBLevP2d_sim

    !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrSTBLevP2d_sim(&
      rdiscretisationTrial, rdiscretisationTest, rform, nelements,&
      npointsPerElement, Dpoints, ibct, DpointPar, IdofsTrial,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! 1D Buckley-Leverett equation in space and time.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:), pointer :: DnormalX,DnormalY
    real(DP) :: dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrSTBLevP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,2))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)

      ! Multiply the velocity vector [a,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = DnormalX(ipoint,iel)*Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + DnormalY(ipoint,iel)
          
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale*dnv
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, DnormalX, DnormalY)

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale*BDRC_DIRICHLET_PENALTY
        end do
      end do
      
      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form
      Dcoefficients = 0.0_DP
      
      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrSTBLevP2d_sim')

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,2))
      allocate(DnormalX(npointsPerElement,nelements))
      allocate(DnormalY(npointsPerElement,nelements))

      ! Evaluate the solution field in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
      
      ! Multiply the velocity vector [a,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = DnormalX(ipoint,iel)*Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + DnormalY(ipoint,iel)
          
          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale*dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, DnormalX, DnormalY)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrSTBLevP2d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrSTBLevP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagBurgP2d_sim(DdataAtNode,&
      DcoeffsAtNode, InodeList, dscale, nnodes,&
      DmatrixAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the primal Burger`s equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: InodeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtNode
!</output>
!</subroutine>
    
    ! local variables
    integer :: inode

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficients $k_{ii} = (0.5*u_i*I)*C_{ii}$
      DmatrixAtNode(1,inode) = dscale*&
          0.5_DP*(DdataAtNode(inode)*DcoeffsAtNode(1,inode)&
                 +DdataAtNode(inode)*DcoeffsAtNode(2,inode))
#else
      ! Compute convective coefficients $k_{ii} = -(0.5*u_i*I)*C_{ii}$
      DmatrixAtNode(1,inode) = -dscale*&
          0.5_DP*(DdataAtNode(inode)*DcoeffsAtNode(1,inode)&
                 +DdataAtNode(inode)*DcoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalBurgP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Burger`s equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
    
    ! local variables
    integer :: iedge

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = (0.5*u_j*I)*C_{ji}$
      DmatrixAtEdge(1,iedge) = dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DcoeffsAtEdge(1,2,iedge)+DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = (0.5*u_i*I)*C_{ij}$
      DmatrixAtEdge(2,iedge) = dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DcoeffsAtEdge(1,1,iedge)+DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -(0.5*u_j*I)*C_{ij}$
      DmatrixAtEdge(1,iedge) = -dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DcoeffsAtEdge(1,1,iedge)+DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -(0.5*u_i*I)*C_{ji}$
      DmatrixAtEdge(2,iedge) = -dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DcoeffsAtEdge(1,2,iedge)+DcoeffsAtEdge(2,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwBurgP2d_sim(DdataAtEdge,&
      DcoeffsAtEdge, IedgeList, dscale, nedges,&
      DmatrixAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Burger`s equation in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DmatrixAtEdge
!</output>
!</subroutine>
    
    ! local variables
    integer :: iedge

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = (0.5*u_j*I)*C_{ji}$
      DmatrixAtEdge(2,iedge) = dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DcoeffsAtEdge(1,2,iedge)+DcoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = (0.5*u_i*I)*C_{ij}$
      DmatrixAtEdge(3,iedge) = dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DcoeffsAtEdge(1,1,iedge)+DcoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -(0.5*u_j*I)*C_{ij}$
      DmatrixAtEdge(2,iedge) = -dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DcoeffsAtEdge(1,1,iedge)+DcoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -(0.5*u_i*I)*C_{ji}$
      DmatrixAtEdge(3,iedge) = -dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DcoeffsAtEdge(1,2,iedge)+DcoeffsAtEdge(2,2,iedge))
#endif

      ! Compute artificial diffusion coefficient
      ! $d_{ij} = abs(v_{ij}*I*0.5*(c_{ij}-c_{ji})$,
      ! where $v_{ij} = 0.5*(f(u_j)-f(u_i))/(u_j-u_i) = 0.5*(u_i+u_j)$
      DmatrixAtEdge(1,iedge) = dscale*&
          abs(0.25_DP*(DdataAtEdge(1,iedge)+DdataAtEdge(2,iedge))*&
          (DcoeffsAtEdge(1,1,iedge)-DcoeffsAtEdge(1,2,iedge)+&
           DcoeffsAtEdge(2,1,iedge)-DcoeffsAtEdge(2,2,iedge)))
    end do

  end subroutine transp_calcMatUpwBurgP2d_sim

  !***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrBurgP2d_sim(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the constant velocities in the primal problem.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: DnormalX,DnormalY,DpointParMirror
    real(DP), dimension(:,:,:), pointer :: Dcoords,Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dlocalData,dnv,dscale,dtime,dval
    real(DP) :: dminParam,dmaxParam,dminParamMirror,dmaxParamMirror
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrBurgP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vectors
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
#endif

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrBurgP2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP

      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrBurgP2d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ D\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.

      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale*Dcoefficients(1,ipoint,iel)
        end do
      end do


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ([0.5*u*I]*u\cdot{\bf n}=([0.5*g*I]*g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,1), DcoeffAtDOF(npoints))
#endif
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      do iel = 1, nelements
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
#endif
        
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,1)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Impose Dirichlet value via penalty method
          DcoeffAtDOF(ipoint) = dscale*dval*BDRC_DIRICHLET_PENALTY
#else
          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale*dval*BDRC_DIRICHLET_PENALTY
#endif
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form
      !
      ! $$ -([0.5*u*I]*u-d\nabla u)\cdot{\bf n} = -([0.5*g*I]*g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Allocate temporal memory
      allocate(DnormalX(npoints,nelements), DnormalY(npoints,nelements))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the physical coordinates of the DOFs on the boundary,
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, DnormalX, DnormalY, 1)
#else
      ! Calculate the normal vectors in cubature on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
#endif
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Evaluate the function parser for the boundary values in the
      ! cubature points or the DOFs on the boundary. Multiply the
      ! velocity vector with the normal in each point to get the
      ! normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, dval)
          
          ! Compute the normal velocity
          dnv = 0.5_DP*(DnormalX(ipoint,iel)+DnormalY(ipoint,iel))*dval

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Store the Robin boundary condition at the DOF on the boundary
          DcoeffAtDOF(ipoint) = -dscale*dnv*dval
#else
          ! Store the Robin boundary condition at the cubature point on the boundary
          Dcoefficients(1,ipoint,iel) = -dscale*dnv*dval
#endif
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(DnormalX, DnormalY)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif


    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u*I]*u-d\nabla u)\cdot{\bf n} = -([0.5*g*I]g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(DnormalX(npoints,nelements),DnormalY(npoints,nelements),&
               Daux(npoints,nelements,1))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the solution values and the physical coordinates of the DOFs
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, DnormalX, DnormalY, 1)
#else
      ! Evaluate the solution vector in the cubature points on the
      ! boundary and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
#endif
      
      ! Initialise values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Evaluate the function parser for the boundary values in the
      ! cubature points or the DOFS on the boundary. Multiply the
      ! velocity vector with the normal in each point to get the
      ! normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Compute the normal velocity
          dnv = 0.5_DP*(DnormalX(ipoint,iel) + DnormalY(ipoint,iel))*Daux(ipoint,iel,1)
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = -dscale*0.5_DP*(DnormalX(ipoint,iel) +&
                                                    DnormalY(ipoint,iel))*dval*dval
#else
            Dcoefficients(1,ipoint,iel) = -dscale*0.5_DP*(DnormalX(ipoint,iel) +&
                                                            DnormalY(ipoint,iel))*dval*dval
#endif
          else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = 0.0_DP
#else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
          end if
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux, DnormalX, DnormalY)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif


    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u*I]*u-d\nabla u)\cdot{\bf n} = -([0.5*g*I]*g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      ! Get mirrored boundary region from collection structure
      p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
          'rboundaryRegionMirror', ssectionName=trim(rcollection%SquickAccess(1)))
      
      ! Get minimum/maximum parameter values from collection structure
      dminParam = rcollection%DquickAccess(3)
      dmaxParam = rcollection%DquickAccess(4)
      dminParamMirror = rcollection%DquickAccess(5)
      dmaxParamMirror = rcollection%DquickAccess(6)
      
      ! Allocate temporal memory
      allocate(DnormalX(npoints,nelements), DnormalY(npoints,nelements))
      allocate(DpointParMirror(npoints,nelements), Daux(npoints,nelements,2))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the solution values and the physical coordinates of the DOFs
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, DnormalX, DnormalY, 1)
#else
      ! Evaluate the solution vector in the cubature points on
      ! the boundary and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
#endif
      
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Rescale parameter values DedgePosition on the boundary
        ! segment where to compute the boundary conditions into
        ! parameter values on the mirror boundary region
        call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
            dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region
        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
      else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Rescale parameter values DedgePosition on the boundary
        ! segment where to compute the boundary conditions into
        ! parameter values on the reversely mirror boundary region
        call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
            dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the reversely mirror boundary region
        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
      end if
      
      ! Evaluate the solution in the cubature points on the mirrored
      ! boundary and store the result in Daux(:,:,2)
      call doEvaluateAtBdr2d(DER_FUNC, npoints*nelements,&
          Daux(:,:,2), p_rsolution%RvectorBlock(1), DpointParMirror,&
          ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)

      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal velocity
          dnv = 0.5_DP*(DnormalX(ipoint,iel) + DnormalY(ipoint,iel))*Daux(ipoint,iel,1)
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = -dscale*0.5_DP*(DnormalX(ipoint,iel) +&
                                  DnormalY(ipoint,iel))*Daux(ipoint,iel,2)**2
#else
            Dcoefficients(1,ipoint,iel) = -dscale*0.5_DP*(DnormalX(ipoint,iel) +&
                                          DnormalY(ipoint,iel))*Daux(ipoint,iel,2)**2
#endif
          else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = 0.0_DP
#else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
          end if
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux, DnormalX, DnormalY)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif
      

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrBurgP2d_sim')
      call sys_halt()
      
    end select

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Deallocate temporal memory
    deallocate(Dbas)
#endif

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate the solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    !***************************************************************************
    
    subroutine doEvaluateAtBdr2d(iderType, n, Dvalues, rvectorScalar,&
        DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n
      real(DP), dimension(n), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdr2d

  end subroutine transp_coeffVecBdrBurgP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrBurgP2d_sim( rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the constant velocities in the primal problem.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! An array accepting the DOF`s on all elements in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: DnormalX,DnormalY,Daux
    real(DP), dimension(:,:,:), pointer :: Dcoords
    real(DP) :: dlocalData,dnv,dscale,dtime
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

    ! REMARK: This subroutine makes use of a lot of preprocessor flags
    ! to distinguish between the classical boundary integal evaluation,
    ! that is, the boundary integrals are approximated by evaluating
    ! the FE-function at the cubature points and computing the flux
    ! based on the FE-function in the cubature points.
    ! This approach is not fully compatible with the group finite
    ! element formulation which approximates the fluxes in the same
    ! way as the approxitame solution in the degrees of freedom. In
    ! this case, the fluxes are evaluated in the degrees of freedom
    ! and interpolated to the cubature points on the boundary afterwards.

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrBurgP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1


#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
#endif
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    ! which will be used to (a) interpolate the FE-functions evaluated
    ! at the degrees of freedom into the cubature points for numerical
    ! integration and to redistribute the function values given at the
    ! cubature points to the degrees of freefom.
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrBurgP2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral

      ! Allocate temporal memory
      allocate(DnormalX(npoints,nelements), DnormalY(npoints,nelements),&
               Daux(npoints,nelements))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the solution values and the physical coordinates of the DOFs
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel) = p_Ddata(ivt)
                
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel) = p_Ddata(ivt)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, DnormalX, DnormalY, 1)
#else
      ! Evaluate the solution vector in the cubature points on the
      ! boundary and store the result in Daux
      call fevl_evaluate_sim(DER_FUNC, Daux,&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
#endif
      
      ! Multiply the velocity vector [0.5*u*I] with the normal in
      ! each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal velocity
          dnv = 0.5_DP*(DnormalX(ipoint,iel) + DnormalY(ipoint,iel))*Daux(ipoint,iel)
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Scale normal velocity by scaling parameter
          DcoeffAtDOF(ipoint) = -dscale*dnv
#else
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale*dnv
#endif
        end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
        
      ! Free temporal memory
      deallocate(Daux, DnormalX, DnormalY)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif
      
      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale*BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form

      Dcoefficients = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrBurgP2d_sim')


    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions

      ! Allocate temporal memory
      allocate(DnormalX(npoints,nelements), DnormalY(npoints,nelements),&
               Daux(npoints,nelements))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
     
      ! We need the solution values and the physical coordinates of the DOFs
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, DnormalX, DnormalY, 1)
#else
      ! Evaluate the velocity field (if any) in the cubature points on
      ! the boundary and store the result in Daux
      call fevl_evaluate_sim(DER_FUNC, Daux,&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, DnormalX, DnormalY, 1)
#endif
      
      ! Multiply the velocity vector [0.5*u*I] with the normal in
      ! each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal velocity
          dnv = 0.5_DP*(DnormalX(ipoint,iel) + DnormalY(ipoint,iel))*Daux(ipoint,iel)

          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = -dscale*dnv
#else
            Dcoefficients(1,ipoint,iel) = -dscale*dnv
#endif
          else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = 0.0_DP
#else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
          end if
        end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Free temporal memory
      deallocate(Daux, DnormalX, DnormalY)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrBurgP2d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_coeffMatBdrBurgP2d_sim

end module transport_callback2d

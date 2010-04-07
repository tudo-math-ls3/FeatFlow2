!##############################################################################
!# ****************************************************************************
!# <name> elasticity_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the elasticity problem that are
!# used during the matrix/vector assembly. Furthermore, the type t_problem is defined
!# here which contains all necessary information for the problem to be computed. Finally,
!# some constants needed for the computation are defined.
!# </purpose>
!##############################################################################

module elasticity_callback

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use triangulation
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use element
  
  implicit none

  integer, parameter :: EQ_POISSON          = 1
  integer, parameter :: EQ_ELASTICITY       = 2

  integer, parameter :: SIMUL_REAL          = 1
  integer, parameter :: SIMUL_ANALYTICAL    = 2

  integer, parameter :: BC_NEUMANN          = 1
  integer, parameter :: BC_DIRICHLET        = 2

  integer, parameter :: SOLVER_DIRECT       = 1
  integer, parameter :: SOLVER_CG           = 2
  integer, parameter :: SOLVER_BICGSTAB     = 3
  integer, parameter :: SOLVER_MG           = 4
  integer, parameter :: SOLVER_CG_MG        = 5
  integer, parameter :: SOLVER_MG_CG        = 6
  integer, parameter :: SOLVER_MG_BICGSTAB  = 7

  integer, parameter :: SMOOTHER_NO         = 0
  integer, parameter :: SMOOTHER_JACOBI     = 1
  integer, parameter :: SMOOTHER_ILU        = 2


  ! type for storing all necessary information about the simulation
  type t_problem
    !   grid file
    character(len=500) :: sgridFileTri
    character(len=500) :: sgridFilePrm

    ! number of boundaries (has to be set manually by the user who has to know
    ! the number of boundaries in the current grid)
    integer :: nboundaries = 1

    ! number of boundary segments (has to be set manually by the user who has to know
    ! the number segments in the current grid)
    integer, dimension(:), pointer :: NboundarySegments

    ! max. number boundary segments over all boundaries
    integer :: nmaxNumBoundSegments
  
    ! kind of equation (possible values: EQ_POISSON, EQ_ELASTICITY)
    integer :: cequation = EQ_ELASTICITY

    ! number of blocks (1 for Poisson equation, 2 for 2D elasticity equation)
    integer :: nblocks
  
    ! material parameters (Poisson ratio nu and shear modulus mu)
    real(DP) :: dnu     = 0.3_DP
    real(DP) :: dmu     = 0.5_DP
    real(DP) :: dlambda = 0.75_DP

    ! definition of boundary conditions (BC_NEUMANN or BC_DIRICHLET)
    ! (dimension  nblocks x max. number segments x nboundaries)
    integer, dimension(:,:,:), pointer :: Cbc

    ! type of simulation (possible values: SIMUL_REAL, SIMUL_ANALYTICAL)
    integer :: csimulation = SIMUL_REAL

    ! given surface forces on Neumann boundary condition segments
    ! (dimension nblocks x max. number segments x nboundaries)
    real(DP), dimension(:,:,:), pointer :: DforceSurface

    ! constant RHS values (only needed in case of csimulation .eq. SIMUL_REAL)
    real(DP) :: dforceVolumeX   = 0.0_DP
    real(DP) :: dforceVolumeY   = 0.0_DP

    ! function IDs (only needed in case of csimulation .eq. SIMUL_ANALYTICAL)
    integer :: cfuncID_u1 = 4
    integer :: cfuncID_u2 = 52

    ! kind of element used (possible values: EL_Q1, EL_Q2)
    integer :: celement

    ! 1D and 2D cubature formulas (they are automatically chosen according to the
    ! selected finite element)
    integer :: ccubature1D, ccubature2D
  
    ! MAX & MIN level where we want to solve.
    integer :: ilevelMax, ilevelMin

    ! kind of solver (possible values: SOLVER_DIRECT,BICGSTAB_SOLVER,SOLVER_MG,SOLVER_CG)
    integer :: csolver = SOLVER_DIRECT

    ! flag whether MG is involved
    logical :: bmgInvolved = .FALSE.
    
    ! max. number of iterations
    integer :: niterations = 5000

    ! tolerance
    real(DP) :: dtolerance = 1e-8_DP

    ! kind of elementary smoother (possible values: SMOOTHER_JACOBI, SMOOTHER_ILU)
    integer :: celementaryPrec = SMOOTHER_JACOBI

    ! MG cycle (0=F-cycle, 1=V-cycle, 2=W-cycle)
    integer :: ccycle = 0

    ! number of smoothing steps
    integer :: nsmoothingSteps = 2

    ! damping parameter
    real(DP) :: ddamp = 0.7_DP

    ! show deformation in gmv(possible values: YES, NO
    integer :: cshowDeformation = NO

    ! number of points where to evaluate the finite element solution
    integer :: nevalPoints = 0
  
    ! number of provided reference solutions in the first
    ! min(nevalPoints, nrefSols) evaluation points
    integer :: nrefSols = 0

    ! given points where to evaluate the FE solution (dimension ndim x nevalPoints)
    real(DP), dimension(:,:), pointer :: DevalPoints

    ! given reference solutions in the evaluation points (dimension ndim x nevalPoints)
    real(DP), dimension(:,:), pointer :: DrefSols

    ! values and derivatives of the FE function in the evaluation points
    ! (dimension ndim x nevalPoints)
    real(DP), dimension(:,:), pointer :: Dvalues, DderivX, DderivY

  end type t_problem

  ! the one and only instance of the t_problem type which is used in this and the main
  ! module to store all necessary information about the simulation
  type(t_problem), save :: rprob


contains


! ****************************************************************************************


!<subroutine>
  subroutine elast_mat_Poisson_2D(rdiscretisationTrial, rdiscretisationTest, rform, &
                                  nelements, npointsPerElement,Dpoints, IdofsTrial, &
                                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the matrix assembly. It computes the coefficients
    ! in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these elements
    ! (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine simultaneously computes
    ! for all these points and all the terms in the bilinear form the corresponding
    ! coefficients.
    ! The routine implements the interface defined in the file
    ! 'intf_coefficientMatrixSc.inc'.
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
    
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(#local DOF in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    Dcoefficients = 1.0_DP

  end subroutine elast_mat_Poisson_2D


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_Poisson_2D_vol(rdiscretisation, rform, nelements, &
                                      npointsPerElement, Dpoints, IdofsTest, &
                                      rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the vector assembly. It computes the function f
    ! in the linear form (usually L(v) = (f,v)_0 + (g,v)_N), i.e. the volumetric 
    ! contributions.
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

    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(\#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
    
!</input>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    real(DP), dimension(:,:), pointer :: Der_u1xx, Der_u1yy
    
    allocate(Der_u1xx(npointsPerElement,nelements), &
             Der_u1yy(npointsPerElement,nelements))
    
    call elast_analFunc_u1(DER_DERIV_XX, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u1xx, rcollection)
    
    call elast_analFunc_u1(DER_DERIV_YY, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u1yy, rcollection)
    
    if (rprob%csimulation .eq. SIMUL_REAL) then
      Dcoefficients(1,:,:) = rprob%dforceVolumeX
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      ! compute Laplace operator
      Dcoefficients(1,:,:) = -Der_u1xx - Der_u1yy
    end if
    
    deallocate(Der_u1xx, Der_u1yy)
  
  end subroutine elast_RHS_Poisson_2D_vol


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_Poisson_2D_bound(rdiscretisation, rform, nelements, &
                                        npointsPerElement, Dpoints, ibct, DpointPar, &
                                        IdofsTest, rdomainIntSubset, Dcoefficients, &
                                        rcollection)
    
    use basicgeometry
    use boundary
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    use fparser
    use spatialdiscretisation

!<description>
    ! This subroutine is called during the vector assembly. It computes the function g
    ! in the linear form (usually L(v) = (f,v)_0 + (g,v)_N), i.e. the contributions
    ! stemming from nonzero Neumann boundary conditions.
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
    
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
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
    real(DP), dimension(:,:,:,:), pointer :: DstressTensor
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny
    integer :: iel,ipoint
    allocate(DstressTensor(2,2,npointsPerElement,nelements))

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      end do
    end do
   

    ! multiply the solution vector with the normal in each point to get the normal velocity
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
      
        dt = DpointPar(ipoint,iel) 
        if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then 
          call elast_stressTensor(rdiscretisation, nelements, npointsPerElement, Dpoints,&
                                  IdofsTest, rdomainIntSubset, DstressTensor, rcollection)
          ! Get the normal vector in the point from the boundary.
          ! Note that the parameter value is in length parametrisation!
          ! When we are at the left or right endpoint of the interval, we
          ! calculate the normal vector based on the current edge.
          ! Without that, the behaviour of the routine may lead to some
          ! confusion if the endpoints of the interval coincide with
          ! the endpoints of a boundary edge. In such a case, the routine
          ! would normally compute the normal vector as a mean on the
          ! normal vectors of the edges adjacent to such a point!
          if (DpointPar(ipoint,iel) .eq. dminPar) then
            ! Start point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
               ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
  
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
           ! End point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
               ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
               ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if

          ! Compute the normal value
!          print *,'1comp', ibct,rcollection%IquickAccess(1),dnx,dny
!          print *,Dpoints(1,ipoint,iel), Dpoints(2,ipoint,iel)
          Dcoefficients(1,ipoint,iel) = dnx * DstressTensor(1,1,ipoint,iel) &
                                      + dny * DstressTensor(1,2,ipoint,iel)
 
        else if (rprob%csimulation .eq. SIMUL_REAL) then
          ! in rcollection%IquickAccess(1) the current segment number is stored
          Dcoefficients(1,ipoint,iel) &
            = rprob%DforceSurface(1,rcollection%IquickAccess(1),ibct)
        end if
      end do
    end do

    deallocate(DstressTensor)

  end subroutine elast_RHS_Poisson_2D_bound


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_volX(rdiscretisation, rform, nelements, npointsPerElement, &
                               Dpoints, IdofsTest, rdomainIntSubset, Dcoefficients, &
                               rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the vector assembly. It computes the function f_1
    ! in the linear form (usually L(v) = (f_1,v_1)_0 + (f_2,v_2)_0 + (g,v)_N), i.e. the
    ! volumetric forces in x-direction.
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

    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(\#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
    
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    real(DP), dimension(:,:), pointer :: Der_u1xx,Der_u1yy,Der_u2xy
    
    allocate(Der_u1xx(npointsPerElement,nelements), &
             Der_u1yy(npointsPerElement,nelements), &
             Der_u2xy(npointsPerElement,nelements))
    
    call elast_analFunc_u1(DER_DERIV_XX, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u1xx, rcollection)
    
    call elast_analFunc_u2(DER_DERIV_XY, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u2xy, rcollection)
    
    call elast_analFunc_u1(DER_DERIV_YY, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u1yy, rcollection)
    
    if (rprob%csimulation .eq. SIMUL_REAL) then
      Dcoefficients(1,:,:) = rprob%dforceVolumeX
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      Dcoefficients(1,:,:) = - (2 * rprob%dmu + rprob%dlambda) * Der_u1xx &
                             -  rprob%dmu                      * Der_u1yy &
                             - (rprob%dmu + rprob%dlambda)     * Der_u2xy
    end if
    
    deallocate(Der_u1xx, Der_u1yy, Der_u2xy)
  
  end subroutine elast_RHS_2D_volX

!BRAL: merge volX and volY by using the rcollection%IquickAccess(2) field


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_volY(rdiscretisation, rform, nelements, npointsPerElement, &
                                Dpoints, IdofsTest, rdomainIntSubset, Dcoefficients, &
                                rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the vector assembly. It computes the function f_2
    ! in the linear form (usually L(v) = (f_1,v_1)_0 + (f_2,v_2)_0 + (g,v)_N), i.e. the
    ! volumetric forces in y-direction.
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

    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(\#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
    
!</input>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
    
!</subroutine>
    real(DP), dimension(:,:), pointer :: Der_u2xx,Der_u2yy,Der_u1yx

    allocate(Der_u2xx(npointsPerElement,nelements), &
             Der_u2yy(npointsPerElement,nelements), &
             Der_u1yx(npointsPerElement,nelements))

    call elast_analFunc_u2(DER_DERIV_XX, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u2xx, rcollection)
    
    call elast_analFunc_u1(DER_DERIV_XY, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u1yx,rcollection)
    
    call elast_analFunc_u2(DER_DERIV_YY, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u2yy, rcollection)
    
    if (rprob%csimulation .eq. SIMUL_REAL) then
	    Dcoefficients(1,:,:) = rprob%dforceVolumeY
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
	    Dcoefficients(1,:,:) = - (rprob%dmu + rprob%dlambda)     * Der_u1yx &
	                          -  rprob%dmu                      * Der_u2xx &
	                          - (2 * rprob%dmu + rprob%dlambda) * Der_u2yy
    end if

    deallocate(Der_u2xx, Der_u2yy, Der_u1yx)

  end subroutine elast_RHS_2D_volY


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_surfX(rdiscretisation, rform, nelements, npointsPerElement, &
                                Dpoints, ibct, DpointPar, IdofsTest, rdomainIntSubset, &
                                Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    use fparser
    use spatialdiscretisation

!<description>
    ! This subroutine is called during the vector assembly. It computes the function g_1
    ! in the linear form (usually L(v) = (f,v)_0 + (g_1,v_1)_N + (g_2,v_2)_N), i.e. the
    ! surface forces in x-direction (nonzero Neumann boundary conditions of the first
    ! block row in the block system).
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
    
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
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
    real(DP), dimension(:,:,:,:), pointer :: DstressTensor
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny
    integer :: iel,ipoint
    allocate(DstressTensor(2,2,npointsPerElement,nelements))

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      end do
    end do
   

    ! Multiply the velocity vector with the normal in each point
    ! to get the normal velocity.
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        dt = DpointPar(ipoint,iel) 
        if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then 
          call elast_stressTensor(rdiscretisation, nelements, npointsPerElement, Dpoints,&
                                  IdofsTest, rdomainIntSubset, DstressTensor, rcollection)
          ! Get the normal vector in the point from the boundary.
          ! Note that the parameter value is in length parametrisation!
          ! When we are at the left or right endpoint of the interval, we
          ! calculate the normal vector based on the current edge.
          ! Without that, the behaviour of the routine may lead to some
          ! confusion if the endpoints of the interval coincide with
          ! the endpoints of a boundary edge. In such a case, the routine
          ! would normally compute the normal vector as a mean on the
          ! normal vectors of the edges adjacent to such a point!
          if (DpointPar(ipoint,iel) .eq. dminPar) then
            ! Start point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
          
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! End point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if

          ! Compute the normal value
!          print *,'1comp', ibct,rcollection%IquickAccess(1),dnx,dny
!          print *,Dpoints(1,ipoint,iel), Dpoints(2,ipoint,iel)
          Dcoefficients(1,ipoint,iel) = dnx * DstressTensor(1,1,ipoint,iel) &
                                      + dny * DstressTensor(1,2,ipoint,iel)
 
        else if (rprob%csimulation .eq. SIMUL_REAL) then
          ! in rcollection%IquickAccess(1) the current segment number is stored
          Dcoefficients(1,ipoint,iel) &
            = rprob%DforceSurface(1,rcollection%IquickAccess(1),ibct)
        end if
      end do
    end do

    deallocate(DstressTensor)

  end subroutine elast_RHS_2D_surfX

!BRAL: merge surfX and surfY by using the rcollection%IquickAccess(2) field
  

! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_surfY(rdiscretisation,rform, nelements, npointsPerElement, &
                                 Dpoints, ibct, DpointPar, IdofsTest, rdomainIntSubset, &
                                 Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    use fparser
    use spatialdiscretisation

!<description>
    ! This subroutine is called during the vector assembly. It computes the function g_2
    ! in the linear form (usually L(v) = (f,v)_0 + (g_1,v_1)_N + (g_2,v_2)_N), i.e. the
    ! surface forces in y-direction (nonzero Neumann boundary conditions of the first
    ! block row in the block system).
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
    
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(#local DOF in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
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
    real(DP), dimension(:,:,:,:), pointer :: DstressTensor
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny
    integer :: iel,ipoint
    allocate(DstressTensor(2,2,npointsPerElement,nelements))
  
    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      end do
    end do
  
     ! Multiply the velocity vector with the normal in each point
     ! to get the normal velocity.
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
    
        dt = DpointPar(ipoint,iel)
  
!BRAL: put this outside the do-loops
        if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then 
          call elast_stressTensor(rdiscretisation, nelements, npointsPerElement, Dpoints,&
                                  IdofsTest, rdomainIntSubset, DstressTensor, rcollection)
          ! Get the normal vector in the point from the boundary.
          ! Note that the parameter value is in length parametrisation!
          ! When we are at the left or right endpoint of the interval, we
          ! calculate the normal vector based on the current edge.
          ! Without that, the behaviour of the routine may lead to some
          ! confusion if the endpoints of the interval coincide with
          ! the endpoints of a boundary edge. In such a case, the routine
          ! would normally compute the normal vector as a mean on the
          ! normal vectors of the edges adjacent to such a point!
          if (DpointPar(ipoint,iel) .eq. dminPar) then
            ! Start point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
    
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! End point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if
  
         ! Compute the normal value
  !       print *,'2comp', ibct,rcollection%IquickAccess(1),dnx,dny
  !       print *,Dpoints(1,ipoint,iel), Dpoints(2,ipoint,iel)
          Dcoefficients(1,ipoint,iel) = dnx * DstressTensor(2,1,ipoint,iel)&
                                      + dny * DstressTensor(2,2,ipoint,iel)
     
        else if (rprob%csimulation .eq. SIMUL_REAL) then
           ! in rcollection%IquickAccess(1) the current segment number is stored
           Dcoefficients(1,ipoint,iel) = &
             rprob%DforceSurface(2,rcollection%IquickAccess(1),ibct)
        end if
      end do
    end do

    deallocate(DstressTensor)

  end subroutine elast_RHS_2D_surfY


! ****************************************************************************************

  
!<subroutine>
  subroutine elast_stressTensor(rdiscretisation, nelements, npointsPerElement, Dpoints, &
                                IdofsTest, rdomainIntSubset, Dvalues, rcollection)
  
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
  
!<description>
!</description>
  
!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(\#local DOF in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest
  
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
  
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    real(DP), dimension(:,:), pointer :: Der_u1x,Der_u2x,Der_u1y,Der_u2y
    allocate(Der_u1x(npointsPerElement,nelements), &
             Der_u2x(npointsPerElement,nelements), &
             Der_u1y(npointsPerElement,nelements), &
             Der_u2y(npointsPerElement,nelements))
    
    
    call elast_analFunc_u1(DER_DERIV_X,rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u1x, rcollection)
    
    call elast_analFunc_u1(DER_DERIV_Y,rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest,rdomainIntSubset, Der_u1y, rcollection)
    
    call elast_analFunc_u2(DER_DERIV_X, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u2x, rcollection)
    
    call elast_analFunc_u2(DER_DERIV_Y, rdiscretisation, nelements, npointsPerElement, &
                           Dpoints, IdofsTest, rdomainIntSubset, Der_u2y, rcollection)                       
    
    Dvalues(1,1,:,:) =   2 * rprob%dmu * Der_u1x(:,:) &
                       + rprob%dlambda * (Der_u1x(:,:) + Der_u2y(:,:))
    Dvalues(1,2,:,:) = rprob%dmu * (Der_u1y(:,:) + Der_u2x(:,:))
    Dvalues(2,1,:,:) = rprob%dmu * (Der_u2x(:,:) + Der_u1y(:,:))
    Dvalues(2,2,:,:) =   2 * rprob%dmu * Der_u2y(:,:) &
                       + rprob%dlambda * (Der_u1x(:,:) + Der_u2y(:,:))
    
    deallocate(Der_u1x, Der_u2x, Der_u1y, Der_u2y)

  end subroutine elast_stressTensor


! ****************************************************************************************


!<subroutine>
  subroutine elast_analFunc_u1(cderivative, rdiscretisation, nelements, npointsPerElement,&
                               Dpoints, IdofsTest, rdomainIntSubset, Dvalues, rcollection)
  
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine computes the (analytical) values of a function in a couple of points
    ! on a couple of elements.
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
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(\#local DOF in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest
  
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

     Dvalues(:,:) = elast_danalyticFunction(Dpoints, nelements, npointsPerElement, &
                                            cderivative, rprob%cfuncID_u1)

  end subroutine elast_analFunc_u1


!BRAL: merge refFunc_u1 and refFunc_u2 by using the rcollection%IquickAccess(2) field
!      and using an array instead of rprob%cfuncID_u1


! ****************************************************************************************


!<subroutine>
  subroutine elast_analFunc_u2(cderivative, rdiscretisation, nelements, npointsPerElement,&
                               Dpoints, IdofsTest, rdomainIntSubset, Dvalues, rcollection)
  
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
  
!<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the (analytical) values of a function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
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
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
    ! An array accepting the DOF on all elements trial in the trial space.
    ! DIMENSION(\#local DOF in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest
  
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    Dvalues(:,:) = elast_danalyticFunction(Dpoints, nelements, npointsPerElement, &
                                           cderivative, rprob%cfuncID_u2)

  end subroutine elast_analFunc_u2


! ****************************************************************************************


!<subroutine>
  subroutine elast_boundValue_2D(Icomponents, rdiscretisation, rboundaryRegion, &
                                 ielement, cinfoNeeded, iwhere, dwhere, Dvalues, &
                                 rcollection)
  
    use collection
    use spatialdiscretisation
    use discretebc
    
!<description>
    ! This subroutine is called during the discretisation of boundary conditions. It
    ! calculates a special quantity on the boundary, which is then used by the
    ! discretisation routines to generate a discrete 'snapshot' of the (actually analytic)
    ! boundary conditions.
!</description>
  
!<input>
    ! Component specifierFor Dirichlet boundary: 
    !   Icomponents(1) defines the number of the boundary component, the value
    !   is to be calculated for (e.g. 1=1st solution component, i.e. x-displacement, 
    !   2=2nd solution component, i.e. y-displacement,...)
    integer, dimension(:), intent(in) :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in) :: rboundaryRegion
    
    ! The element number on the boundary which is currently being processed
    integer, intent(in) :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in) :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    integer, intent(in) :: iwhere
  
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    real(DP), intent(in) :: dwhere
  
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
  
!</input>

!<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    real(DP), dimension(:), intent(out) :: Dvalues
!</output>
  
!</subroutine>

    real(DP), dimension(2,1,1) :: Dpoints
    real(DP), dimension(1,1) :: Daux

    ! To get the X/Y-coordinates of the boundary point, use:
    real(DP) :: dx,dy
    
    call boundary_getCoords(rdiscretisation%p_rboundary, &
                            rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! set coordinates of the current point
    Dpoints(1,1,1) = dx
    Dpoints(2,1,1) = dy
    ! calculate value
    if (rprob%csimulation .eq. SIMUL_REAL) then
      ! currently, only zero Dirichlet boundary values are supported for real simulations
! BRAL: das muss verbessert werden!
      Dvalues(1) = 0.0_DP
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      select case (Icomponents(1))
      case(1) ! x-displacement
        Daux = elast_danalyticFunction(Dpoints,1,1, DER_FUNC, rprob%cfuncID_u1)
        Dvalues(1) = Daux(1,1)
      case(2) ! y-displacement
        Daux = elast_danalyticFunction(Dpoints,1,1, DER_FUNC, rprob%cfuncID_u2)
        Dvalues(1) = Daux(1,1)
      end select
    end if
  end subroutine elast_boundValue_2D


! ****************************************************************************************


!<function>
  function elast_danalyticFunction(Dpts, nelements, npointsPerElement, cderiv, cselect, &
                                   dparam) result(Dval)

!<description>
    ! This function provides some analytic functions, which can be used for
    ! validating your FE code.
!</description>

!<input>
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpts

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
  
    ! derivative of the function to be calculated
    integer(I32), intent(in) :: cderiv
  
    ! selector for the desired function
    integer(I32), intent(in) :: cselect
  
    ! optional parameter to influence the solution
    real(DP), intent(in), optional :: dparam
!</input>

!<!--
    !
    ! cselect   u(x,y)
    !   0        0.0
    !   1        0.1 * x
    !   2        0.1 * x^2
    !   3        0.1 * y
    !   4        0.1 * y^2
    !   5        4 * x * (1 - x)
    !   6        4 * y * (1 - y)
    !   7        x * (1 - x) * y * (1 - y)    (zero on the boundary of the unit square;
    !                                          used in FBENCHMARK)
    !   8        -(x * x + y * y) * x         (used in FBENCHMARK)
    !   9        y^2 * (1 - y)^2 * x
    !  10        y^2 * (1 - y)^2 * (x - 1)
    !  11        0.25 *  (1/sqrt(2) + x - y) * (1/sqrt(2) - x + y)
    !                 * ((1-sqrt(2))/sqrt(2) + x + y) * ((1+sqrt(2))/sqrt(2) - x - y)
    !                                         (zero on the boundary of the unit square
    !                                          rotated by Pi/4)
    !  12        sin(x) * sin(y)
    !  13        0.05 * sin(4 * PI * x) * sin(4 * PI * y)
    !            (zero on the boundary of the unit square)
    !  14        cos(x) * cos(y)
    !  15        cos(PI/2 * (x + y))     (zero divergence)
    !  16        -cos(PI/2 * (x + y))    (zero divergence)
    !  17        2 * cos(x) * sin(y) - 2 * (1 - cos(1)) * sin(1)
    !                        (has zero integral over the boundary of the unit square)
    !  18        8 * (1 - x) (Stokes: pressure for parabolic inflow (6) on unit square)
    !  19        sin(PI/2 * (x - y))  (Stokes: to use as pressure with (15)/(16))
    !  20        -(x * x + y * y) * x * x (slight modification of (8); used in Poisson app.)
    !  21        cos(Pi/2)*x - sin(Pi/2)*y - x ((21) and (22) simulate a rigid body
    !  22        sin(Pi/2)*x + cos(Pi/2)*y - y  rotation of Pi/2 around the origin)
    !  23        1.0
    !  24        -(2x - 1)(2y^3 - 3y^2) / 6    (zero divergence together with 6.0_DP x (7))
    !  25        sin(x) cos(y)
    !  26        -cos(x) sin(y)
    !  27        4 - 8x          (Stokes, unit square, zero mean pressure: with /6/ + /0/)
    !  28        2 cos(x) sin(y) - 2 sin(1) + 2 sin(1) cos(1)
    !  29        xy - 1/4        (Stokes, unit square, zero mean pressure: with /7/ + /24/)
    !
    !  Function triple taken from Bochev/ Gunzburger/ Lehoucq, On stabilized finite
    !  element methods for the Stokes problem in the small time limit (preprint)
    !  30        sin(PI*x - 7/10) * sin(PI*y + 1/5)
    !  31        cos(PI*x - 7/10) * cos(PI*y + 1/5)
    !  32        sin(x) * cos(y) + (cos(1) - 1) * sin(1)
    !
    !  33        -sin(gamma x) * sin(gamma y)
    !  34        -cos(gamma x) * cos(gamma y)
    !            (gamma should be set to a multiple of pi. gamma=pi is just fine.
    !             A setting of gamma=3*Pi gives in combination with stabilised Q1/Q1
    !             approaches serious problems for the pressure solution of coarse grids.)
    !
    !  Function triple taken from Bochev/ Gunzburger/ Lehoucq, On stabilized finite
    !  element methods for the Stokes problem in the small time limit (journal version)
    !  (together with 32 for p). The difference to the preprint version (30,31) is that
    !  the velocities are not only divergence free but also zero on the boundary of the
    !  unitsquare.
    !  35        x^2*(1-x)^2 * 2*PI*sin(PI*y)*cos(PI*y)
    !  36        -(2*x*(1-x)^2 - 2*x^2*(1-x))*sin^2(Pi*y)
    !            (these are the first two components of curl(0,0,psi) where
    !             psi(x,y) = x^2(1-x)^2 * sin^2(PI*y) )
    !
    !  37        16 * x * (1 - x) * y * (1 - y)
    !  38        42*x^2 * (2 - y) + 42*sin(2 - 5*y*x) * cos(x + y + 1)
    !  39        (x^2 - 1)^2 (y^2 - 1)y / 4      (zero divergence together with 40)
    !  40        (y^2 - 1)^2 (1 - x^2)x / 4
    !  41        5 x^3 (y-1) + y^3
    !  42        x*(y-0.01) (for solid beam configuration; use together with 43)
    !  43        -0.5*x^2 (for solid beam configuration; use together with 42)
    !  44        sin(c1*x+c2) * sin(c1*y+c2)  (variant of 12/30)
    !  45        cos(c1*x+c2) * cos(c1*y+c2)  (variant of 14/31)
    !  46        -sin(c1*x+c2) * cos(c1*y+c2)  (variant of 25/32)
    !  47        sin(PI*x+0.4) cos(PI*y-0.3) (variant of 25/32)
    !  48        c*x^2*y*sin(c*(x-0.5*y^2))
    !  49        -2*x*cos(c*(x-0.5*y^2)) + c*x^2*sin(c*(x-0.5*y^2))
    !  50        0.05 * sin(2*PI*x)*sin(2*PI*y)
    !            (same as 13, only different factors)
    !  51        sin(PI/2 (x-1)) sin(PI/2 (y-1))
    !            (Circular harmonic function on x^2 + y^2 - 1 = 0)
    !  52        0.05 * cos(2*PI*x)*cos(2*PI*y)
    !            (same as 50, only sin replaced by cos in order to get nonzero values on
    !             the boundary of the unitsquare)


    ! Stokes pairs (for unit square with zero mean pressure):
    ! / 0,  0,  0/ , /23, 23,  0/ , / 6,  0, 27/ , /15, 16, 19/,
    ! /12, 14, 28/ , / 7, 24, 29/ , /25, 26, 19/, /30, 31, 32/, /33, 34, 32/, /35, 36, 32/
    !
    ! Stokes pairs (for square box [-1,1]x[-1,1], zero mean pressure):
    ! /39, 40, 41/
    !
! -->

!<result>
    ! (The result of the function calculation)
    real(DP), dimension(npointsPerElement,nelements) :: Dval
!</result>

!<errors>
    ! none
!</errors>
!</function>

    real(DP) :: daux, daux1

    real(DP) :: dgamma
#ifdef SOLUTION_CAUSING_SERIOUS_PROBLEMS_FOR_STABILISED_Q1_Q1_STOKES
    ! gamma should be set to a multiple of pi. gamma=pi is just fine.
    ! A setting of gamma=3*Pi gives in combination with stabilised Q1/Q1 approaches
    ! serious problems for the pressure solution of coarse grids.
    dgamma = 3.0_DP * SYS_PI
#else
    dgamma = 1.0_DP * SYS_PI
#endif

    ! avoid misleading warnings about uninitialised variables
    Dval(:,:) = 0.0_DP

    select case (cselect)

    case (0) ! u(x,y) = 0.0
      Dval(:,:) = 0.0_DP

    case (1) ! u(x,y) = 0.1 * x
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.1_DP
      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (2) ! u(x,y) = 0.1 * x**2
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(1,:,:) * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.2_DP * Dpts(1,:,:)
      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dval(:,:) = 0.2_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (3) ! u(x,y) = 0.1 * y
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(2,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_Y);  Dval(:,:) = 0.1_DP
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (4) ! u(x,y) = 0.1 * y**2
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 0.1_DP * Dpts(2,:,:) * Dpts(2,:,:)
      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_Y);  Dval(:,:) = 0.2_DP * Dpts(2,:,:)
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.2_DP
      end select

    case (5) ! u(x,y) = 4 * x * (1 - x)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 4.0_DP * Dpts(1,:,:) * (1.0_DP - Dpts(1,:,:))
      case (DER_DERIV_X);  Dval(:,:) = 4.0_DP * (1.0_DP - 2.0_DP * Dpts(1,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_XX); Dval(:,:) = -8.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
      end select

    case (6) ! u(x,y) = 4 * y * (1 - y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = 4.0_DP * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
      case (DER_DERIV_Y);  Dval(:,:) = 4.0_DP * (1.0_DP - 2.0_DP * Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
      case (DER_DERIV_YY); Dval(:,:) = -8.0_DP
      end select

    case (7) ! u(x,y) = x * (1 - x) * y * (1 - y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        Dpts(1,:,:) * (1.0_DP - Dpts(1,:,:)) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = &
        (-1.0_DP + 2.0_DP * Dpts(1,:,:)) * Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = &
        Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:)) * (-1.0_DP + 2.0_DP * Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) = (2.0_DP * Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:)))
      case (DER_DERIV_XY); Dval(:,:) = &
        (-1.0_DP + 2.0_DP * Dpts(1,:,:)) * (-1.0_DP + 2.0_DP * Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) = 2.0_DP * Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:))
      end select

    case (8) ! u(x,y) = -(x * x + y * y) * x
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        -(Dpts(1,:,:) * Dpts(1,:,:) + Dpts(2,:,:) * Dpts(2,:,:)) * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = &
        -(3.0_DP * Dpts(1,:,:) * Dpts(1,:,:) + Dpts(2,:,:) * Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) = -2.0_DP * Dpts(1,:,:) * Dpts(2,:,:)
      case (DER_DERIV_XX); Dval(:,:) = -6.0_DP * Dpts(1,:,:)
      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP * Dpts(2,:,:)
      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP * Dpts(1,:,:)
      end select

    case (9) ! u(x,y) = y^2 * (1 - y)^2 * x
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 * Dpts(1,:,:)
      case (DER_DERIV_X);  Dval(:,:) = &
        Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2
      case (DER_DERIV_Y);  Dval(:,:) = &
          2.0_DP * Dpts(1,:,:) * (Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 - Dpts(2,:,:) &
        * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:)))
      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
      case (DER_DERIV_XY); Dval(:,:) = &
          2.0_DP * (Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 &
        - Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:)))
      case (DER_DERIV_YY); Dval(:,:) = &
        2.0_DP * Dpts(1,:,:) * ((1.0_DP - Dpts(2,:,:))**2 &
        - 4.0_DP * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:)) + Dpts(2,:,:) * Dpts(2,:,:))
      end select

!    case (10) ! u(x,y) = y^2 * (1 - y)^2 * (x - 1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2 * (Dpts(1,:,:) - 1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) = Dpts(2,:,:) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))**2
!      case (DER_DERIV_Y);  Dval(:,:) =   2.0_DP * Dpts(2,:,:) * (2.0_DP * Dpts(2,:,:) - 1.0_DP) &
!                              * (Dpts(2,:,:) - 1.0_DP) * (Dpts(1,:,:) - 1.0_DP)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 2.0_DP * Dpts(2,:,:) * (2.0_DP * Dpts(2,:,:) - 1.0_DP) * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_YY); Dval(:,:) =   2.0_DP * (1.0_DP - 6.0_DP * Dpts(2,:,:) + 6.0_DP * Dpts(2,:,:) * Dpts(2,:,:)) &
!                              * (Dpts(2,:,:) - 1.0_DP)
!      end select
!
!    case (11) ! u1(x,y) = 0.25 * (1/sqrt(2) + x - y) * (1/sqrt(2) - x + y)
!              !            * ((1-sqrt(2))/sqrt(2) + x + y) * ((1+sqrt(2))/sqrt(2) - x - y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   0.25_DP * (1.0_DP / sqrt(2.0_DP) + Dpts(1,:,:) - Dpts(2,:,:)) &
!                              * (1.0_DP / sqrt(2.0_DP) - Dpts(1,:,:) + Dpts(2,:,:)) &
!                              * ((1.0_DP - sqrt(2.0_DP)) / sqrt(2.0_DP) + Dpts(1,:,:) + Dpts(2,:,:)) &
!                              * ((1.0_DP + sqrt(2.0_DP)) / sqrt(2.0_DP) - Dpts(1,:,:) - Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  -1.5_DP * Dpts(1,:,:) * Dpts(1,:,:) - 1.0_DP * Dpts(2,:,:) * Dpts(2,:,:) * Dpts(1,:,:) &
!                              + 0.5_DP * Dpts(2,:,:) * Dpts(2,:,:) - 0.5_DP * Dpts(2,:,:) + 0.25_DP &
!                              + 1.0_DP * Dpts(1,:,:) * Dpts(2,:,:) + 1.0_DP * Dpts(1,:,:)**3
!      case (DER_DERIV_Y);  Dval(:,:) =  -1.0_DP * Dpts(1,:,:) * Dpts(1,:,:) * Dpts(2,:,:) + 0.5_DP*Dpts(1,:,:) * Dpts(1,:,:) &
!                              - 0.5_DP*Dpts(1,:,:) - 1.5_DP * Dpts(2,:,:) * Dpts(2,:,:) + 0.25_DP &
!                              + 1.0_DP * Dpts(1,:,:) * Dpts(2,:,:) + 1.0_DP * Dpts(2,:,:)**3
!      case (DER_DERIV_XX); Dval(:,:) = -1.0_DP * Dpts(2,:,:) * Dpts(2,:,:) + 1.0_DP * Dpts(2,:,:) + 3.0_DP * Dpts(1,:,:) * Dpts(1,:,:)&
!                              - 3.0_DP * Dpts(1,:,:)
!      case (DER_DERIV_XY); Dval(:,:) = -0.5_DP + 1.0_DP * Dpts(1,:,:) + 1.0_DP * Dpts(2,:,:) - 2.0_DP * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_YY); Dval(:,:) =  3.0_DP * Dpts(2,:,:) * Dpts(2,:,:) - 3.0_DP * Dpts(2,:,:) - 1.0_DP * Dpts(1,:,:) * Dpts(1,:,:) &
!                               + 1.0_DP * Dpts(1,:,:)
!      end select
!
!    case (12) ! u(x,y) = sin(x) * sin(y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      end select
!
!    case (13) ! u(x,y) = 0.05 * sin(4*PI*x)*sin(4*PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  0.05_DP*sin(4.0_DP * SYS_PI * Dpts(1,:,:)) * sin(4.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =   0.2_DP * SYS_PI &
!                              * cos(4.0_DP * SYS_PI * Dpts(1,:,:)) * sin(4.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =   0.2_DP * SYS_PI &
!                              * sin(4.0_DP * SYS_PI * Dpts(1,:,:)) * cos(4.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  -0.8_DP * SYS_PI * SYS_PI &
!                              * sin(4.0_DP * SYS_PI * Dpts(1,:,:)) * sin(4.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =   0.8_DP * SYS_PI * SYS_PI &
!                              * cos(4.0_DP * SYS_PI * Dpts(1,:,:)) * cos(4.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  -0.8_DP * SYS_PI * SYS_PI &
!                              * sin(4.0_DP * SYS_PI * Dpts(1,:,:)) * sin(4.0_DP * SYS_PI * Dpts(2,:,:))
!      end select
!
!    case (14) ! u(x,y) = cos(x) * cos(y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =  sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      end select
!
!    case (15) ! u(x,y) = cos(PI/2 * (x + y))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                              cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_X);  Dval(:,:) =           -0.5_DP * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_Y);  Dval(:,:) =           -0.5_DP * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XX); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_YY); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      end select
!
!    case (16) ! u(x,y) = -cos(PI/2 * (x + y))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                            -cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_X);  Dval(:,:) =           0.5_DP * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_Y);  Dval(:,:) =           0.5_DP * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XX); Dval(:,:) = 0.25_DP * SYS_PI * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) = 0.25_DP * SYS_PI * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      case (DER_DERIV_YY); Dval(:,:) = 0.25_DP * SYS_PI * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) + Dpts(2,:,:)))
!      end select
!
!    case (17) ! u(x,y) = 2 * cos(x) * sin(y) - 2 * (1 - cos(1)) * sin(1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  2.0_DP * cos(Dpts(1,:,:)) * sin(Dpts(2,:,:)) &
!                              -2.0_DP * (1.0_DP - cos(1.0_DP)) * sin(1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) = -2.0_DP * sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  2.0_DP * cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -2.0_DP * cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP * sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP * cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      end select
!
!    case (18) ! u(x,y) = 8 * (1 - x)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  8.0_DP*(1.0_DP - Dpts(1,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = -8.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (19) ! u(x,y) = sin(PI/2 * (x - y))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                              sin(0.5_DP * SYS_PI * (Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_X);  Dval(:,:) =            0.5_DP * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_Y);  Dval(:,:) =           -0.5_DP * SYS_PI * cos(0.5_DP * SYS_PI * (Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_XX); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) =  0.25_DP * SYS_PI * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) - Dpts(2,:,:)))
!      case (DER_DERIV_YY); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI * sin(0.5_DP * SYS_PI * (Dpts(1,:,:) - Dpts(2,:,:)))
!      end select
!
!    case (20) ! u(x,y) = -(x * x + y * y) * x * x
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = -(Dpts(1,:,:) * Dpts(1,:,:) + Dpts(2,:,:) * Dpts(2,:,:)) * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = -4.0_DP * Dpts(1,:,:)**3 + 2.0_DP * Dpts(1,:,:) * Dpts(2,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) = -2.0_DP * Dpts(1,:,:) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XX); Dval(:,:) = -12.0_DP * Dpts(1,:,:) * Dpts(1,:,:) - 2.0_DP * Dpts(2,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XY); Dval(:,:) = -4.0_DP * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP * Dpts(1,:,:) * Dpts(1,:,:)
!      end select
!
!    case (21) ! u(x,y) = cos(Pi/2)*x - sin(Pi/2)*y - x
!      if (present(dparam)) then
!        daux = 0.125_DP * dparam * SYS_PI
!      else
!        daux = 0.5_DP * SYS_PI
!      endif
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = cos(daux) * Dpts(1,:,:) - sin(daux) * Dpts(2,:,:) - Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = cos(daux) - 1.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) = - sin(daux)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (22) ! u(x,y) = sin(Pi/2)*x + cos(Pi/2)*y - y
!      if (present(dparam)) then
!        daux = 0.125_DP * dparam * SYS_PI
!      else
!        daux = 0.5_DP * SYS_PI
!      endif
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = sin(daux) * Dpts(1,:,:) + cos(daux) * Dpts(2,:,:) - Dpts(2,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = sin(daux)
!      case (DER_DERIV_Y);  Dval(:,:) = cos(daux) - 1.0_DP
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (23) ! u(x,y) = 1.0
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 1.0_DP
!      case (DER_DERIV_X);  Dval(:,:) = 0.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (24) ! u(x,y) = -(2x - 1)(2y^3 - 3y^2) / 6
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = - (2.0_DP * Dpts(1,:,:) - 1.0_DP) * (2.0_DP * Dpts(2,:,:)**3 - 3.0_DP * Dpts(2,:,:)**2) / 6.0_DP
!      case (DER_DERIV_X);  Dval(:,:) = - Dpts(2,:,:)**2 * (2.0_DP * Dpts(2,:,:) - 3.0_DP) / 3.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) = -Dpts(2,:,:) * (Dpts(2,:,:) - 1.0_DP) * (2.0_DP * Dpts(1,:,:) - 1.0_DP)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP * Dpts(2,:,:) * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_YY); Dval(:,:) = -(2.0_DP * Dpts(1,:,:) - 1.0_DP) * (2.0_DP * Dpts(2,:,:) - 1.0_DP)
!      end select
!
!    case (25) ! u(x,y) = sin(x) cos(y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      end select
!
!    case (26) ! u(x,y) = -cos(x) sin(y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      end select
!
!    case (27) ! u(x,y) = 4 - 8x
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  4.0_DP - 8.0_DP * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = -8.0_DP
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (28) ! u(x,y) = 2 cos(x) sin(y) - 2 sin(1) + 2 sin(1) cos(1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   2.0_DP * cos(Dpts(1,:,:)) * sin(Dpts(2,:,:)) &
!                                - 2.0_DP * sin(1.0_DP) &
!                                + 2.0_DP * sin(1.0_DP) * cos(1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) = -2.0_DP * sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  2.0_DP * cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -2.0_DP * cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP * sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -2.0_DP * cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      end select
!
!    case (29) ! u(x,y) = xy - 1/4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(1,:,:) * Dpts(2,:,:) - 0.25_DP
!      case (DER_DERIV_X);  Dval(:,:) = Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) = Dpts(1,:,:)
!      case (DER_DERIV_XX); Dval(:,:) = 0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) = 1.0_DP
!      case (DER_DERIV_YY); Dval(:,:) = 0.0_DP
!      end select
!
!    case (30) ! u(x,y) = sin(PI*x - 7/10) * sin(PI*y + 1/5)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(SYS_PI * Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI * Dpts(2,:,:) + 0.2_DP)
!      case (DER_DERIV_X);  Dval(:,:) =  cos(SYS_PI * Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_Y);  Dval(:,:) =  sin(SYS_PI * Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -sin(SYS_PI * Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI * SYS_PI
!      case (DER_DERIV_XY); Dval(:,:) =  cos(SYS_PI * Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI * SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) = -sin(SYS_PI * Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI * SYS_PI
!      end select
!
!    case (31) ! u(x,y) = cos(PI*x - 7/10) cos(PI*y + 1/5)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  cos(SYS_PI * Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI * Dpts(2,:,:) + 0.2_DP)
!      case (DER_DERIV_X);  Dval(:,:) = -sin(SYS_PI * Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_Y);  Dval(:,:) = -cos(SYS_PI * Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -cos(SYS_PI * Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI * SYS_PI
!      case (DER_DERIV_XY); Dval(:,:) =  sin(SYS_PI * Dpts(1,:,:) - 0.7_DP) * sin(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI * SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) = -cos(SYS_PI * Dpts(1,:,:) - 0.7_DP) * cos(SYS_PI * Dpts(2,:,:) + 0.2_DP) * SYS_PI * SYS_PI
!      end select
!
!    case (32) ! u(x,y) = sin(x) cos(y) + (cos(1) - 1) sin(1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(Dpts(1,:,:)) * cos(Dpts(2,:,:)) + (cos(1.0_DP) - 1.0_DP) * sin(1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) =  cos(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -sin(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -cos(Dpts(1,:,:)) * sin(Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = -sin(Dpts(1,:,:)) * cos(Dpts(2,:,:))
!      end select
!
!    case (33) ! u(x,y) = -sin(gamma x) * sin(gamma y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                   -sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = -dgamma *          cos(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = -dgamma *          sin(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  dgamma * dgamma * sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -dgamma * dgamma * cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  dgamma * dgamma * sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      end select
!
!    case (34) ! u(x,y) = -cos(gamma x) * cos(gamma y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                   -cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =  dgamma *          sin(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  dgamma *          cos(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  dgamma * dgamma * cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = -dgamma * dgamma * sin(dgamma * Dpts(1,:,:)) * sin(dgamma * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  dgamma * dgamma * cos(dgamma * Dpts(1,:,:)) * cos(dgamma * Dpts(2,:,:))
!      end select
!
!    case (35) ! u(x,y) = x^2*(1-x)^2 * 2*PI*sin(PI*y)*cos(PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))**2*2.0_DP*SYS_PI*sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = 4.0_DP*(Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 - Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) &
!                                * SYS_PI*sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = 2.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))**2 &
!                                * (SYS_PI**2*cos(SYS_PI*Dpts(2,:,:))**2 - SYS_PI**2*sin(SYS_PI*Dpts(2,:,:))**2)
!      case (DER_DERIV_XX); Dval(:,:) = (4.0_DP*(1.0_DP-Dpts(1,:,:))**2 - 16.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:)) + 4.0_DP*Dpts(1,:,:)**2) &
!                                * SYS_PI*sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) = 4.0_DP*SYS_PI**2 * ((Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 &
!                                - Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) * cos(SYS_PI*Dpts(2,:,:))**2 &
!                                - (Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 + Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) * sin(SYS_PI*Dpts(2,:,:))**2)
!      case (DER_DERIV_YY); Dval(:,:) =-8.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))**2 &
!                                * SYS_PI**3*cos(SYS_PI*Dpts(2,:,:))*sin(SYS_PI*Dpts(2,:,:))
!      end select
!
!    case (36) ! u(x,y) = -(2*x*(1-x)^2 - 2*x^2*(1-x))*sin^2(PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = -(2.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 &
!                                - 2.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:)))*sin(SYS_PI*Dpts(2,:,:))**2
!      case (DER_DERIV_X);  Dval(:,:) = -(2.0_DP*(1.0_DP-Dpts(1,:,:))**2 - 8.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:)) &
!                                + 2.0_DP*Dpts(1,:,:)**2) * sin(SYS_PI*Dpts(2,:,:))**2
!      case (DER_DERIV_Y);  Dval(:,:) = -2.0_DP*(2.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2-2.0_DP*Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:))) &
!                                * sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))*SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -(-12.0_DP + 24.0_DP*Dpts(1,:,:))*sin(SYS_PI*Dpts(2,:,:))**2
!      case (DER_DERIV_XY); Dval(:,:) = -2.0_DP*(2.0_DP*(1.0_DP-Dpts(1,:,:))**2 -8.0_DP*Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))&
!                                + 2.0_DP*Dpts(1,:,:)**2) * sin(SYS_PI*Dpts(2,:,:))*cos(SYS_PI*Dpts(2,:,:))*SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) =   4.0_DP*(Dpts(1,:,:)*(1.0_DP-Dpts(1,:,:))**2 - Dpts(1,:,:)**2*(1.0_DP-Dpts(1,:,:)))*SYS_PI**2 &
!                                * (-cos(SYS_PI*Dpts(2,:,:))**2 + sin(SYS_PI*Dpts(2,:,:))**2)
!      end select
!
!    case (37) ! u(x,y) = 16 * x * (1 - x) * y * (1 - y)
!      Dval(:,:) = 0_DP
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = Dpts(1,:,:) * (1.0_DP - Dpts(1,:,:)) * Dpts(2,:,:) * (1.0_DP - Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) = (-1.0_DP + 2.0_DP * Dpts(1,:,:)) * Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) = Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:)) * (-1.0_DP + 2.0_DP * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) = (2.0_DP * Dpts(2,:,:) * (-1.0_DP + Dpts(2,:,:)))
!      case (DER_DERIV_XY); Dval(:,:) = (-1.0_DP + 2.0_DP * Dpts(1,:,:)) * (-1.0_DP + 2.0_DP * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) = 2.0_DP * Dpts(1,:,:) * (-1.0_DP + Dpts(1,:,:))
!      end select
!      Dval(:,:) = Dval(:,:)*16.0_DP
!
!    case (38) ! u(x,y) = 42*x^2*(2 - y) + 42*sin(2 - 5*y*x)*cos(x + y + 1)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 42.0_DP * Dpts(1,:,:)** 2 * (2.0_DP - Dpts(2,:,:)) &
!                                + 42.0_DP * sin(2.0_DP - 5.0_DP * Dpts(2,:,:)* Dpts(1,:,:)) * cos(Dpts(1,:,:) + Dpts(2,:,:) + 1.0_DP)
!      case (DER_DERIV_X);  Dval(:,:) =   84.0_DP * Dpts(1,:,:)* (2.0_DP - Dpts(2,:,:)) &
!                                - 210.0_DP * cos(-2.0_DP + 5.0_DP * Dpts(2,:,:)* Dpts(1,:,:)) * Dpts(2,:,:) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                - 42.0_DP * sin(2.0_DP - 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =  -42.0_DP * Dpts(1,:,:) ** 2 &
!                                - 210.0_DP * cos(-2.0_DP + 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(1,:,:) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                - 42.0_DP * sin(2.0_DP - 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =   168.0_DP - 84.0_DP * Dpts(2,:,:) &
!                                + 1050.0_DP * sin(-2.0_DP + 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) ** 2 * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 420.0_DP * cos(-2.0_DP + 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                -  42.0_DP * sin(2.0_DP - 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =   -84.0_DP * Dpts(2,:,:) + 1050.0_DP * sin(-2.0_DP &
!                                + 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) * Dpts(1,:,:) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                - 210.0_DP * cos(-2.0_DP + 5.0_DP* Dpts(2,:,:)* Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 210.0_DP * cos(-2.0_DP + 5.0_DP* Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(1,:,:) * sin(1.0_DP +Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 210.0_DP * cos(-2.0_DP + 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(2,:,:) * sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                -  42.0_DP * sin(2.0_DP - 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  1050.0_DP * sin(-2.0_DP &
!                                + 5.0_DP * Dpts(2,:,:)* Dpts(1,:,:)) * Dpts(1,:,:) ** 2 * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                + 420.0_DP * cos(-2.0_DP + 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * Dpts(1,:,:)* sin(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:)) &
!                                -  42.0_DP * sin(2.0_DP - 5.0_DP * Dpts(2,:,:) * Dpts(1,:,:)) * cos(1.0_DP + Dpts(1,:,:) + Dpts(2,:,:))
!      end select
!
!    case (39) ! u(x,y) = (x^2 - 1)^2 (y^2 - 1)y / 4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 0.25_DP * (Dpts(1,:,:)**2 - 1.0_DP)**2 * (Dpts(2,:,:)**2 - 1.0_DP) * Dpts(2,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = (Dpts(1,:,:)**2 - 1.0_DP) * (Dpts(2,:,:)**2 - 1.0_DP) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) = 0.25_DP * (Dpts(1,:,:)**2 - 1.0_DP)**2 * (3.0_DP * Dpts(2,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_XX); Dval(:,:) = Dpts(2,:,:) * (Dpts(2,:,:)**2 - 1.0_DP) * (3.0_DP * Dpts(1,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_XY); Dval(:,:) = Dpts(1,:,:) * (Dpts(1,:,:)**2 - 1.0_DP) * (3.0_DP * Dpts(2,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_YY); Dval(:,:) = 1.50_DP * (Dpts(1,:,:)**2 - 1)**2 * Dpts(2,:,:)
!      end select
!
!    case (40) ! u(x,y) = (y^2 - 1)^2 (1 - x^2)x / 4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 0.25_DP * (Dpts(2,:,:)**2 - 1.0_DP)**2 * (1.0_DP - Dpts(1,:,:)**2) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) = -0.25_DP * (Dpts(2,:,:)**2 - 1.0_DP)**2 * (3.0_DP * Dpts(1,:,:)**2 - 1.0_DP)
!      case (DER_DERIV_Y);  Dval(:,:) = (Dpts(2,:,:)**2 - 1.0_DP) * (1.0_DP - Dpts(1,:,:)**2) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XX); Dval(:,:) = -1.50_DP * (Dpts(2,:,:)**2 - 1)**2 * Dpts(1,:,:)
!      case (DER_DERIV_XY); Dval(:,:) = Dpts(2,:,:) * (Dpts(2,:,:)**2 - 1.0_DP) * (1.0_DP - 3.0_DP * Dpts(1,:,:)**2)
!      case (DER_DERIV_YY); Dval(:,:) = Dpts(1,:,:) * (Dpts(1,:,:)**2 - 1.0_DP) * (1.0_DP - 3.0_DP * Dpts(2,:,:)**2)
!      end select
!
!    case (41) ! u(x,y) = 5 x^3 (y-1) + y^3
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) = 5.0_DP * Dpts(1,:,:)**3 * (Dpts(2,:,:) - 1.0_DP) + Dpts(2,:,:)**3
!      case (DER_DERIV_X);  Dval(:,:) = 15.0_DP * Dpts(1,:,:)**2 * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_Y);  Dval(:,:) = 5.0_DP * Dpts(1,:,:)**3 + 3.0_DP * Dpts(2,:,:)**2
!      case (DER_DERIV_XX); Dval(:,:) = 30.0_DP * Dpts(1,:,:) * (Dpts(2,:,:) - 1.0_DP)
!      case (DER_DERIV_XY); Dval(:,:) = 15.0_DP * Dpts(1,:,:)**2
!      case (DER_DERIV_YY); Dval(:,:) = 6.0 * Dpts(2,:,:)
!      end select
!
!    case (42) ! u(x,y) = x*(y-0.01)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  Dpts(1,:,:) * (Dpts(2,:,:)-0.01_DP)
!      case (DER_DERIV_X);  Dval(:,:) =      (Dpts(2,:,:)-0.01_DP)
!      case (DER_DERIV_Y);  Dval(:,:) =  Dpts(1,:,:)
!      case (DER_DERIV_XX); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  1.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (43) ! u(x,y) = -0.5*x^2
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  -0.5_DP * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) =  -Dpts(1,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  -1.0_DP
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (44) ! u(x,y) = sin(c1*x + c2) * sin(c1*y + c2)
!      daux = 2.0_DP*SYS_PI
!      daux1 = 0.02_DP*daux
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_X);  Dval(:,:) =  daux *        cos(daux * Dpts(1,:,:)+ daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_Y);  Dval(:,:) =  daux *        sin(daux *Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XX); Dval(:,:) = -daux * daux * sin(daux * Dpts(1,:,:)+ daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XY); Dval(:,:) =  daux * daux * cos(daux * Dpts(1,:,:)+ daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_YY); Dval(:,:) = -daux * daux * sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      end select
!
!    case (45) ! u(x,y) = cos(c1*x + c2) * cos(c1*y + c2)
!      daux = 2.0_DP*SYS_PI
!      daux1 = 0.02_DP*daux
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =                cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_X);  Dval(:,:) = -daux *        sin(daux *Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_Y);  Dval(:,:) = -daux *        cos(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XX); Dval(:,:) = -daux * daux * cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:)+ daux1)
!      case (DER_DERIV_XY); Dval(:,:) =  daux * daux * sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_YY); Dval(:,:) = -daux * daux * cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      end select
!
!    case (46) ! u(x,y) = -sin(c1*x + c2) * cos(c1*y + c2)
!      daux = (2.0_DP*SYS_PI)**2
!      daux1 = 0.02_DP*daux
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =               -sin(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_X);  Dval(:,:) = -daux *        cos(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_Y);  Dval(:,:) =  daux *        sin(daux * Dpts(1,:,:) + daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XX); Dval(:,:) =  daux * daux * sin(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_XY); Dval(:,:) =  daux * daux * cos(daux * Dpts(1,:,:)+ daux1) * sin(daux * Dpts(2,:,:) + daux1)
!      case (DER_DERIV_YY); Dval(:,:) =  daux * daux * sin(daux * Dpts(1,:,:) + daux1) * cos(daux * Dpts(2,:,:) + daux1)
!      end select
!
!    case (47) ! u(x,y) = sin(PI*x+0.4) cos(PI*y-0.3)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP)
!      case (DER_DERIV_X);  Dval(:,:) =  cos(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI
!      case (DER_DERIV_Y);  Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * sin(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI
!      case (DER_DERIV_XX); Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI * SYS_PI
!      case (DER_DERIV_XY); Dval(:,:) = -cos(SYS_PI*Dpts(1,:,:) + 0.4_DP) * sin(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI * SYS_PI
!      case (DER_DERIV_YY); Dval(:,:) = -sin(SYS_PI*Dpts(1,:,:) + 0.4_DP) * cos(SYS_PI*Dpts(2,:,:) - 0.3_DP) * SYS_PI * SYS_PI
!      end select
!
!
!    case (48) ! u(x,y) = c*x^2*y*sin(c*(x-0.5*y^2))
!      daux = 20.0_DP*SYS_PI
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   daux*Dpts(1,:,:)**2*Dpts(2,:,:) * sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_X);  Dval(:,:) =   2*daux*Dpts(1,:,:)*Dpts(2,:,:)*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_Y);  Dval(:,:) =  -daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)**2*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + daux*Dpts(1,:,:)**2*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XX); Dval(:,:) =   (2*daux*Dpts(2,:,:) - daux**3*Dpts(1,:,:)**2*Dpts(2,:,:))*sin(daux*(Dpts(1,:,:) &
!                                -0.5_DP*Dpts(2,:,:)**2)) + 4*daux**2*Dpts(1,:,:)*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XY); Dval(:,:) =   (daux**2*Dpts(1,:,:)**2 - 2*daux**2*Dpts(1,:,:)*Dpts(2,:,:)**2)*cos(daux*(Dpts(1,:,:) &
!                                -0.5_DP*Dpts(2,:,:)**2)) + (2*daux*Dpts(1,:,:) + daux**3*Dpts(1,:,:)**2*Dpts(2,:,:)**2)*sin(daux*(Dpts(1,:,:)&
!                                -0.5_DP*Dpts(2,:,:)**2)) 
!      case (DER_DERIV_YY); Dval(:,:) =  -daux**3*Dpts(1,:,:)**2*Dpts(2,:,:)**3*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                - 3*daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      end select
!
!    case (49) ! u(x,y) = -2*x*cos(c*(x-0.5*y^2)) + c*x^2*sin(c*(x-0.5*y^2))
!      daux = 20.0_DP*SYS_PI
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =   -2*Dpts(1,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + daux*Dpts(1,:,:)**2*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_X);  Dval(:,:) =   (daux**2*Dpts(1,:,:)**2 - 2)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + 4*daux*Dpts(1,:,:)*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) 
!      case (DER_DERIV_Y);  Dval(:,:) =  -2*daux*Dpts(1,:,:)*Dpts(2,:,:)*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                - daux**2*Dpts(1,:,:)**2*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XX); Dval(:,:) =   (6*daux - daux**3*Dpts(1,:,:)**2) * sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2)) &
!                                + 6*daux**2*Dpts(1,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_XY); Dval(:,:) =   (daux**3*Dpts(1,:,:)**2*Dpts(2,:,:) - 2*daux*Dpts(2,:,:))*sin(daux*(Dpts(1,:,:)& 
!                                -0.5_DP*Dpts(2,:,:)**2)) - 4*daux**2*Dpts(1,:,:)*Dpts(2,:,:)*cos(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      case (DER_DERIV_YY); Dval(:,:) =   (2*daux**2*Dpts(1,:,:)*Dpts(2,:,:)**2  - daux**2*Dpts(1,:,:)**2)*cos(daux*(Dpts(1,:,:)&
!                                -0.5_DP*Dpts(2,:,:)**2)) - (daux**3*Dpts(1,:,:)**2*Dpts(2,:,:)**2 &
!                                + 2*daux*Dpts(1,:,:))*sin(daux*(Dpts(1,:,:)-0.5_DP*Dpts(2,:,:)**2))
!      end select
!      
!    case (50) ! u(x,y) = 0.05 * sin(2*PI*x)*sin(2*PI*y)
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  0.05_DP*sin(2.0_DP * SYS_PI * Dpts(1,:,:)) * sin(2.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_X);  Dval(:,:) =   0.1_DP * SYS_PI &
!                              * cos(2.0_DP * SYS_PI * Dpts(1,:,:)) * sin(2.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_Y);  Dval(:,:) =   0.1_DP * SYS_PI &
!                              * sin(2.0_DP * SYS_PI * Dpts(1,:,:)) * cos(2.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_XX); Dval(:,:) =  -0.2_DP * SYS_PI * SYS_PI &
!                              * sin(2.0_DP * SYS_PI * Dpts(1,:,:)) * sin(2.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_XY); Dval(:,:) =   0.2_DP * SYS_PI * SYS_PI &
!                              * cos(2.0_DP * SYS_PI * Dpts(1,:,:)) * cos(2.0_DP * SYS_PI * Dpts(2,:,:))
!      case (DER_DERIV_YY); Dval(:,:) =  -0.2_DP * SYS_PI * SYS_PI &
!                              * sin(2.0_DP * SYS_PI * Dpts(1,:,:)) * sin(2.0_DP * SYS_PI * Dpts(2,:,:))
!      end select
!      
!    case (51) ! u(x,y) = sin(PI/2 (x-1)) sin(PI/2 (y-1))
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  sin(0.5_DP * SYS_PI * (Dpts(1,:,:)-1)) * sin(0.5_DP * SYS_PI * (Dpts(2,:,:)-1))
!      case (DER_DERIV_X);  Dval(:,:) =  0.5_DP * SYS_PI &
!                               * cos(0.5_DP * SYS_PI * (Dpts(1,:,:)-1)) * sin(0.5_DP * SYS_PI * (Dpts(2,:,:)-1))
!      case (DER_DERIV_Y);  Dval(:,:) =  0.5_DP * SYS_PI &
!                               * sin(0.5_DP * SYS_PI * (Dpts(1,:,:)-1)) * cos(0.5_DP * SYS_PI * (Dpts(2,:,:)-1))
!      case (DER_DERIV_XX); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI &
!                               * sin(0.5_DP * SYS_PI * (Dpts(1,:,:)-1)) * sin(0.5_DP * SYS_PI * (Dpts(2,:,:)-1))
!      case (DER_DERIV_XY); Dval(:,:) =  0.25_DP * SYS_PI * SYS_PI &
!                               * cos(0.5_DP * SYS_PI * (Dpts(1,:,:)-1)) * cos(0.5_DP * SYS_PI * (Dpts(2,:,:)-1))
!      case (DER_DERIV_YY); Dval(:,:) = -0.25_DP * SYS_PI * SYS_PI &
!                               * sin(0.5_DP * SYS_PI * (Dpts(1,:,:)-1)) * sin(0.5_DP * SYS_PI * (Dpts(2,:,:)-1))
!      end select
!
    case (52) ! u(x,y) = 0.05 * cos(2*PI*x)*cos(2*PI*y)
      select case (cderiv)
      case (DER_FUNC);     Dval(:,:) = &
        0.05_DP*cos(2.0_DP * SYS_PI *Dpts(1,:,:)) * cos(2.0_DP * SYS_PI * Dpts(2,:,:))
      case (DER_DERIV_X);  Dval(:,:) = -0.1_DP * SYS_PI &
        * sin(2.0_DP * SYS_PI * Dpts(1,:,:)) * cos(2.0_DP * SYS_PI * Dpts(2,:,:))
      case (DER_DERIV_Y);  Dval(:,:) =   -0.1_DP * SYS_PI &
        * cos(2.0_DP * SYS_PI * Dpts(1,:,:)) * sin(2.0_DP * SYS_PI * Dpts(2,:,:))
      case (DER_DERIV_XX); Dval(:,:) =  -0.2_DP * SYS_PI * SYS_PI &
        * cos(2.0_DP * SYS_PI * Dpts(1,:,:)) * cos(2.0_DP * SYS_PI * Dpts(2,:,:))
      case (DER_DERIV_XY); Dval(:,:) =   0.2_DP * SYS_PI * SYS_PI &
        * sin(2.0_DP * SYS_PI * Dpts(1,:,:)) * sin(2.0_DP * SYS_PI * Dpts(2,:,:))
      case (DER_DERIV_YY); Dval(:,:) =  -0.2_DP * SYS_PI * SYS_PI &
        * cos(2.0_DP * SYS_PI *Dpts(1,:,:)) * cos(2.0_DP * SYS_PI * Dpts(2,:,:))
      end select
!
!    case (53) ! u(x,y) = -x^3*y
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  -Dpts(1,:,:) * Dpts(1,:,:)*Dpts(1,:,:)*Dpts(2,:,:)
!      case (DER_DERIV_X);  Dval(:,:) =  -3.0_DP * Dpts(1,:,:) * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) =  -Dpts(1,:,:) * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_XX); Dval(:,:) =  -6.0_DP * Dpts(1,:,:) * Dpts(2,:,:)
!      case (DER_DERIV_XY); Dval(:,:) =  -3.0_DP * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!
!    case (54) ! u(x,y) = 1/3*x^4
!      select case (cderiv)
!      case (DER_FUNC);     Dval(:,:) =  (1.0_DP/4.0_DP) * Dpts(1,:,:) * Dpts(1,:,:) * &
!                                            Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_X);  Dval(:,:) =  Dpts(1,:,:) * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_Y);  Dval(:,:) =  0.0_DP
!      case (DER_DERIV_XX); Dval(:,:) =  3.0_DP * Dpts(1,:,:) * Dpts(1,:,:)
!      case (DER_DERIV_XY); Dval(:,:) =  0.0_DP
!      case (DER_DERIV_YY); Dval(:,:) =  0.0_DP
!      end select
!      
    end select

  end function elast_danalyticFunction

end module elasticity_callback


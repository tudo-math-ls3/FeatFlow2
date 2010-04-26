!#########################################################################################
!# ***************************************************************************************
!# <name> elasticity_callback </name>
!# ***************************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the elasticity problem that are
!# used during the matrix/vector assembly.
!# </purpose>
!#########################################################################################

module elasticity_callback

  implicit none

contains


! ****************************************************************************************


!<subroutine>
  subroutine elast_mat_Poisson_2D(rdiscretisationTrial, rdiscretisationTest, rform, &
                                  nelements, npointsPerElement,Dpoints, IdofsTrial, &
                                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use fsystem
    use spatialdiscretisation
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
    
    use fsystem
    use spatialdiscretisation
    use derivatives
    use collection
    use scalarpde
    use domainintegration
    use elasticity_basic
    
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

    real(DP), dimension(:,:), pointer :: Duxx, Duyy
    
    allocate(Duxx(npointsPerElement,nelements), Duyy(npointsPerElement,nelements))
    call elast_analFunc(1, DER_DERIV_XX, rdiscretisation, nelements, npointsPerElement, &
                        Dpoints, rdomainIntSubset, Duxx, rcollection)
    call elast_analFunc(1, DER_DERIV_YY, rdiscretisation, nelements, npointsPerElement, &
                        Dpoints, rdomainIntSubset, Duyy, rcollection)
    
    if (rprob%csimulation .eq. SIMUL_REAL) then
      Dcoefficients(1,:,:) = rprob%dforceVolumeX
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      ! compute Laplace operator
      Dcoefficients(1,:,:) = -Duxx - Duyy
    endif
    
    deallocate(Duxx, Duyy)
  
  end subroutine elast_RHS_Poisson_2D_vol


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_Poisson_2D_bound(rdiscretisation, rform, nelements, &
                                        npointsPerElement, Dpoints, ibct, DpointPar, &
                                        IdofsTest, rdomainIntSubset, Dcoefficients, &
                                        rcollection)
    
    use fsystem
    use boundary
    use collection
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use derivatives
    use elasticity_basic

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
    real(DP) :: dminPar, dmaxPar, dt, dnx, dny
    integer :: iel, ipoint, iseg
    real(DP), dimension(:,:), pointer :: Dux, Duy

    ! get the minimum and maximum parameter value, corresponding to start and end point
    ! of the interval
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      enddo
    enddo

    if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then 
      ! in case of an analytical test function, the Neumann contributions are given by
      ! grad(u) * n, where n is the normal vector

      ! compute grad(u)
      allocate(Dux(npointsPerElement,nelements), Duy(npointsPerElement,nelements))
      call elast_analFunc(1, DER_DERIV_X, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Dux, rcollection)
      call elast_analFunc(1, DER_DERIV_Y, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Duy, rcollection)
  
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          ! get the parameter value corresponding to the current boundary point
          dt = DpointPar(ipoint,iel) 
          ! Compute the normal vector in the current boundary point. At the endpoints of
          ! the interval the normal vector is based on the current edge.
          if (DpointPar(ipoint,iel) .eq. dminPar) then
            ! start point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
               ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
  
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! end point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
               ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
               ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          endif
  
          ! compute the normal value grad(u)*n
          Dcoefficients(1,ipoint,iel) = dnx * Dux(ipoint,iel) + dny * Duy(ipoint,iel)
   
        enddo
      enddo
      deallocate(Dux, Duy)
    else if (rprob%csimulation .eq. SIMUL_REAL) then
      ! in case of a real simulation, constant nonzero Neumann contributions are provided
      ! by the user for each boundary segment

      ! in rcollection%IquickAccess(2) the current segment number is stored
      iseg = rcollection%IquickAccess(2)
      Dcoefficients(1,ipoint,iel) = rprob%DbcValue(1, iseg, ibct)
    endif


  end subroutine elast_RHS_Poisson_2D_bound


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_vol(rdiscretisation, rform, nelements, npointsPerElement, &
                              Dpoints, IdofsTest, rdomainIntSubset, Dcoefficients, &
                              rcollection)
    
    use fsystem
    use collection
    use scalarpde
    use domainintegration
    use derivatives
    use elasticity_basic

!<description>
    ! This subroutine is called during the vector assembly. It computes the function f
    ! in the linear form (usually L(v) = (f,v)_0 + (g,v)_N), i.e. the volumetric forces
    ! When rcollection%IquickAccess(1) .eq. 1, then forces in x-direction are computed,
    ! when rcollection%IquickAccess(1) .eq. 2, then forces in y-direction.
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

    real(DP), dimension(:,:), pointer :: Du1xx, Du1yy, Du2xy, Du2xx, Du2yy, Du1yx
    
    if (rcollection%IquickAccess(1) .eq. 1) then
      ! forces in x-direction
      allocate(Du1xx(npointsPerElement,nelements), &
               Du1yy(npointsPerElement,nelements), &
               Du2xy(npointsPerElement,nelements))
      
      call elast_analFunc(1, DER_DERIV_XX, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Du1xx, rcollection)
      
      call elast_analFunc(2, DER_DERIV_XY, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Du2xy, rcollection)
      
      call elast_analFunc(1, DER_DERIV_YY, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Du1yy, rcollection)
      if (rprob%csimulation .eq. SIMUL_REAL) then
        Dcoefficients(1,:,:) = rprob%dforceVolumeX
      else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
        Dcoefficients(1,:,:) = - (2 * rprob%dmu + rprob%dlambda) * Du1xx &
                               -  rprob%dmu                      * Du1yy &
                               - (rprob%dmu + rprob%dlambda)     * Du2xy
      endif

      deallocate(Du1xx, Du1yy, Du2xy)
  
    else if (rcollection%IquickAccess(1) .eq. 2) then
      ! forces in y-direction
      allocate(Du2xx(npointsPerElement,nelements), &
               Du2yy(npointsPerElement,nelements), &
               Du1yx(npointsPerElement,nelements))
  
      call elast_analFunc(2, DER_DERIV_XX, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Du2xx, rcollection)
      
      call elast_analFunc(1, DER_DERIV_XY, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Du1yx,rcollection)
      
      call elast_analFunc(2, DER_DERIV_YY, rdiscretisation, nelements, npointsPerElement, &
                          Dpoints, rdomainIntSubset, Du2yy, rcollection)
      
      if (rprob%csimulation .eq. SIMUL_REAL) then
        Dcoefficients(1,:,:) = rprob%dforceVolumeY
      else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
        Dcoefficients(1,:,:) = - (rprob%dmu + rprob%dlambda)     * Du1yx &
                               -  rprob%dmu                      * Du2xx &
                               - (2 * rprob%dmu + rprob%dlambda) * Du2yy
      endif
  
      deallocate(Du2xx, Du2yy, Du1yx)
      
    endif
    
  end subroutine elast_RHS_2D_vol


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_surf(rdiscretisation, rform, nelements, npointsPerElement, &
                               Dpoints, ibct, DpointPar, IdofsTest, rdomainIntSubset, &
                               Dcoefficients, rcollection)
    
    use fsystem
    use boundary
    use collection
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use elasticity_basic

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
    real(DP) :: dminPar, dmaxPar, dt, dnx, dny
    integer :: iel, ipoint, icomp, iseg

    ! get the minimum and maximum parameter value, corresponding to start and end point
    ! of the interval
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      enddo
    enddo

    ! in rcollection%IquickAccess(1) the current component is stored
    icomp = rcollection%IquickAccess(1)

    if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      ! in case of an analytical test function, the Neumann contributions are given by
      ! sigma * n, where sigma is the Cauchy stress tensor and n the normal vector
      
      ! compute the stress tensor sigma
      allocate(DstressTensor(2,2,npointsPerElement,nelements))
      call elast_stressTensor(rdiscretisation, nelements, npointsPerElement, Dpoints,&
                              rdomainIntSubset, DstressTensor, rcollection)
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          ! get the parameter value corresponding to the current boundary point
          dt = DpointPar(ipoint,iel) 
          ! Compute the normal vector in the current boundary point. At the endpoints of
          ! the interval the normal vector is based on the current edge.
          if (DpointPar(ipoint,iel) .eq. dminPar) then
            ! start point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
          
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! end point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, &
                                         dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          endif

          ! Now compute sigma * n
          Dcoefficients(1,ipoint,iel) = dnx * DstressTensor(icomp,1,ipoint,iel) &
                                      + dny * DstressTensor(icomp,2,ipoint,iel)
        enddo
      enddo
      deallocate(DstressTensor)

    else if (rprob%csimulation .eq. SIMUL_REAL) then
      ! in case of a real simulation, the Neumann contributions are given by constant
      ! surface forces on the current segment, which are provided by the user

      ! in rcollection%IquickAccess(2) the current segment number is stored
      iseg = rcollection%IquickAccess(2)
      Dcoefficients(1,:,:) = rprob%DbcValue(icomp, iseg, ibct)
    endif


  end subroutine elast_RHS_2D_surf


! ****************************************************************************************

  
!<subroutine>
  subroutine elast_stressTensor(rdiscretisation, nelements, npointsPerElement, Dpoints, &
                                rdomainIntSubset, Dvalues, rcollection)

!BRAL: nur eine Zeile des Spannungstensors berechnen (also icomp uebergeben)!!
!BRAL: an mixed formulation anpassen

    use fsystem
    use collection
    use scalarpde
    use domainintegration
    use derivatives
    use elasticity_basic

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

    real(DP), dimension(:,:), pointer :: Du1x, Du2x, Du1y, Du2y
    allocate(Du1x(npointsPerElement,nelements), Du2x(npointsPerElement,nelements), &
             Du1y(npointsPerElement,nelements), Du2y(npointsPerElement,nelements))
    
    
    call elast_analFunc(1, DER_DERIV_X,rdiscretisation, nelements, npointsPerElement, &
                        Dpoints, rdomainIntSubset, Du1x, rcollection)
    
    call elast_analFunc(1, DER_DERIV_Y,rdiscretisation, nelements, npointsPerElement, &
                        Dpoints, rdomainIntSubset, Du1y, rcollection)
    
    call elast_analFunc(2, DER_DERIV_X, rdiscretisation, nelements, npointsPerElement, &
                        Dpoints, rdomainIntSubset, Du2x, rcollection)
    
    call elast_analFunc(2, DER_DERIV_Y, rdiscretisation, nelements, npointsPerElement, &
                        Dpoints, rdomainIntSubset, Du2y, rcollection)                       
    
    Dvalues(1,1,:,:) =   2*rprob%dmu * Du1x(:,:) + rprob%dlambda * (Du1x(:,:) + Du2y(:,:))
    Dvalues(1,2,:,:) = rprob%dmu * (Du1y(:,:) + Du2x(:,:))
    Dvalues(2,1,:,:) = rprob%dmu * (Du2x(:,:) + Du1y(:,:))
    Dvalues(2,2,:,:) =   2*rprob%dmu * Du2y(:,:) + rprob%dlambda * (Du1x(:,:) + Du2y(:,:))
    
    deallocate(Du1x, Du2x, Du1y, Du2y)

  end subroutine elast_stressTensor


! ****************************************************************************************


!<subroutine>
  subroutine elast_boundValue_2D(Icomponents, rdiscretisation, rboundaryRegion, &
                                 ielement, cinfoNeeded, iwhere, dwhere, Dvalues, &
                                 rcollection)
  
    use fsystem
    use collection
    use spatialdiscretisation
    use boundary
    use derivatives
    use elasticity_basic

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

    ! coordinates of the boundary point
    real(DP) :: dx, dy
    
    ! calculate value
    if (rprob%csimulation .eq. SIMUL_REAL) then
      ! currently, only zero Dirichlet boundary values are supported for real simulations
! BRAL: das muss verbessert werden!
      Dvalues(1) = 0.0_DP
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      ! set coordinates of the current point
      call boundary_getCoords(rdiscretisation%p_rboundary, &
                              rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
      Dpoints(1,1,1) = dx
      Dpoints(2,1,1) = dy
      Daux = elast_danalyticFunction(Dpoints,1,1, DER_FUNC, rprob%CfuncID(Icomponents(1)))
      Dvalues(1) = Daux(1,1)
    endif
  end subroutine elast_boundValue_2D

end module elasticity_callback


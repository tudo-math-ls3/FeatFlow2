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
                                  nel, nptsPerEl,Dpoints, IdofsTrial, IdofsTest, &
                                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use fsystem
    use spatialdiscretisation
    use collection
    use scalarpde
    use domainintegration
    
!<description>
    ! This subroutine is called during the matrix assembly. It computes the coefficients
    ! in front of the terms of the bilinear form.
    ! The routine accepts a set of elements and a set of points on these elements
    ! (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine simultaneously computes
    ! for all these points and all the terms in the bilinear form the corresponding
    ! coefficients.
    ! The routine implements the interface defined in the file
    ! 'intf_coefficientMatrixSc.inc'.
!</description>
    
!<input>
    ! discretisation structure of the trial space
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial
    
    ! discretisation structure of the test space
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! bilinear form to be evaluated
    type(t_bilinearForm), intent(in) :: rform
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! DIMENSION(spatial dimension, nptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    
    ! array of trial space DOF in all elements
    ! DIMENSION(#local DOF in trial space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! array of test space DOF in all elements
    ! DIMENSION(#local DOF in test space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
  
!<output>
    ! array of coefficients of the terms in the bilinear form for all given points on
    ! all given elements, DIMENSION(itermCount, nptsPerEl, nel) with itermCount the
    ! number of terms in the bilinear form
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    Dcoefficients = 1.0_DP

  end subroutine elast_mat_Poisson_2D


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_Poisson_2D_vol(rdiscretisation, rform, nel, nptsPerEl, Dpoints, &
                                      IdofsTest, rdomainIntSubset, Dcoefficients, &
                                      rcollection)
    
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
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! linear form to be evaluated
    type(t_linearForm), intent(in) :: rform
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! DIMENSION(spatial dimension, nptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! array of test space DOF in all elements
    ! DIMENSION(#local DOF in test space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
    
  
!<output>
    ! array of coefficients of the terms in the linear form for all given points on
    ! all given elements, DIMENSION(itermCount, nptsPerEl, nel) with itermCount the
    ! number of terms in the linear form
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    real(DP), dimension(:,:), pointer :: Duxx, Duyy
    
    if (rprob%csimulation .eq. SIMUL_REAL) then
      Dcoefficients(1,:,:) = rprob%dforceVolumeX
    else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
      ! compute Laplace operator
      allocate(Duxx(nptsPerEl,nel), Duyy(nptsPerEl,nel))
      call elast_analFunc(1, DER_DERIV_XX, rdiscretisation, nel, nptsPerEl, Dpoints, &
                          rdomainIntSubset, Duxx, rcollection)
      call elast_analFunc(1, DER_DERIV_YY, rdiscretisation, nel, nptsPerEl, Dpoints, &
                          rdomainIntSubset, Duyy, rcollection)
      Dcoefficients(1,:,:) = -Duxx - Duyy
      deallocate(Duxx, Duyy)
    endif

  end subroutine elast_RHS_Poisson_2D_vol


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_Poisson_2D_bound(rdiscretisation, rform, nel, nptsPerEl, Dpoints, &
                                        ibct, DpointPar, IdofsTest, rdomainIntSubset, &
                                        Dcoefficients, rcollection)
    
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
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! linear form to be evaluated
    type(t_linearForm), intent(in) :: rform
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! DIMENSION(spatial dimension, nptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    
    ! index of the boundary component that contains the points
    integer, intent(in) :: ibct
    
    ! parameter values (LENGTH PARAMETRISATION) of the points on the boundary component
    ! DIMENSION(nptsPerEl, nel)
    real(DP), dimension(:,:), intent(in) :: DpointPar
    
    ! array of test space DOF in all elements
    ! DIMENSION(#local DOF in test space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! array of coefficients of the terms in the linear form for all given points on
    ! all given elements, DIMENSION(itermCount, nptsPerEl, nel) with itermCount the
    ! number of terms in the linear form
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
    do iel = 1, nel
      do ipoint = 1, nptsPerEl
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      enddo
    enddo

    if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then 
      ! in case of an analytical test function, the Neumann contributions are given by
      ! grad(u) * n, where n is the normal vector

      ! compute grad(u)
      allocate(Dux(nptsPerEl,nel), Duy(nptsPerEl,nel))
      call elast_analFunc(1, DER_DERIV_X, rdiscretisation, nel, nptsPerEl, Dpoints, &
                          rdomainIntSubset, Dux, rcollection)
      call elast_analFunc(1, DER_DERIV_Y, rdiscretisation, nel, nptsPerEl, Dpoints, &
                          rdomainIntSubset, Duy, rcollection)
  
      do iel = 1, nel
        do ipoint = 1, nptsPerEl
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
  subroutine elast_RHS_2D_vol(rdiscretisation, rform, nel, nptsPerEl, Dpoints, &
                              IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
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
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! linear form to be evaluated
    type(t_linearForm), intent(in) :: rform
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! DIMENSION(spatial dimension, nptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! array of test space DOF in all elements
    ! DIMENSION(#local DOF in test space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
    

!<output>
    ! array of coefficients of the terms in the linear form for all given points on
    ! all given elements, DIMENSION(itermCount, nptsPerEl, nel) with itermCount the
    ! number of terms in the linear form
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    real(DP), dimension(:,:), pointer :: Du1xx, Du1yx, Du1yy, Du2xx, Du2xy, Du2yy, &
                                         Dpx, Dpy
    
    if (rcollection%IquickAccess(1) .eq. 1) then
      ! forces in x-direction

      if (rprob%csimulation .eq. SIMUL_REAL) then
        Dcoefficients(1,:,:) = rprob%dforceVolumeX

      else ! rprob%csimulation .eq. SIMUL_ANALYTICAL
        allocate(Du1xx(nptsPerEl,nel), Du1yy(nptsPerEl,nel), Du2xy(nptsPerEl,nel))
        call elast_analFunc(1, DER_DERIV_XX, rdiscretisation, nel, nptsPerEl, Dpoints, &
                            rdomainIntSubset, Du1xx, rcollection)
        if (rprob%cformulation .ne. FORMULATION_STOKES) then
          ! for Stokes the mixed derivatives are not required
          call elast_analFunc(2, DER_DERIV_XY, rdiscretisation, nel, nptsPerEl, Dpoints, &
                              rdomainIntSubset, Du2xy, rcollection)
        endif
        call elast_analFunc(1, DER_DERIV_YY, rdiscretisation, nel, nptsPerEl, Dpoints, &
                            rdomainIntSubset, Du1yy, rcollection)
        if (rprob%cformulation .eq. FORMULATION_DISPL) then
          Dcoefficients(1,:,:) = - (2*rprob%dmu + rprob%dlambda) * Du1xx &
                                 -  rprob%dmu                    * Du1yy &
                                 - (rprob%dmu + rprob%dlambda)   * Du2xy
        else
          ! mixed formulation or Stokes
          allocate(Dpx(nptsPerEl,nel))
          call elast_analFunc(3, DER_DERIV_X, rdiscretisation, nel, nptsPerEl, Dpoints, &
                              rdomainIntSubset, Dpx, rcollection)
          if (rprob%cformulation .eq. FORMULATION_MIXED) then
            Dcoefficients(1,:,:) = - 2*rprob%dmu * Du1xx - rprob%dmu * (Du1yy + Du2xy) &
                                   + 2*rprob%dmu * Dpx
            ! Why the factor 2 * mu in front of Dpx?
            ! This code merely provides analytic solutions to test the application.
            ! Typically, mu has order O(1e+6) or even larger. Using a typical analytic
            ! solution provided by elast_danalyticFunction(...) (which is in the range
            ! of O(1)), e.g. some trigonomic function, this would result in a massive
            ! imbalance between u and p:
            !      -2mu * div(eps(u))  +  p    = f
            !          O(1e+6)         + O(1)  = O(1e+6)
            ! Such a problem is not well posed. Results for p would be way off for such
            ! imbalanced equation system. Remedy:
            ! * Use a set of analytic solutions for u and p with u O(1) and p O(mu), or
            ! * Scale the pressure solution here, not in elast_danalyticFunction(...).
          else ! FORMULATION_STOKES
            Dcoefficients(1,:,:) = - 2*rprob%dmu * (Du1xx + Du1yy) + Dpx
          endif
          deallocate(Dpx)
        endif
        deallocate(Du1xx, Du1yy, Du2xy)
      endif
  
    else if (rcollection%IquickAccess(1) .eq. 2) then
      ! forces in y-direction

      if (rprob%csimulation .eq. SIMUL_REAL) then
        Dcoefficients(1,:,:) = rprob%dforceVolumeY

      else if (rprob%csimulation .eq. SIMUL_ANALYTICAL) then
        allocate(Du2xx(nptsPerEl,nel), Du2yy(nptsPerEl,nel), Du1yx(nptsPerEl,nel))
        call elast_analFunc(2, DER_DERIV_XX, rdiscretisation, nel, nptsPerEl, Dpoints, &
                            rdomainIntSubset, Du2xx, rcollection)
        if (rprob%cformulation .ne. FORMULATION_STOKES) then
          ! for Stokes the mixed derivatives are not required
          call elast_analFunc(1, DER_DERIV_XY, rdiscretisation, nel, nptsPerEl, Dpoints, &
                              rdomainIntSubset, Du1yx, rcollection)
        endif
        call elast_analFunc(2, DER_DERIV_YY, rdiscretisation, nel, nptsPerEl, Dpoints, &
                            rdomainIntSubset, Du2yy, rcollection)
        if (rprob%cformulation .eq. FORMULATION_DISPL) then
          Dcoefficients(1,:,:) = - (rprob%dmu + rprob%dlambda)     * Du1yx &
                                 -  rprob%dmu                      * Du2xx &
                                 - (2 * rprob%dmu + rprob%dlambda) * Du2yy
        else
          ! mixed formulation or Stokes
          allocate(Dpy(nptsPerEl,nel))
          call elast_analFunc(3, DER_DERIV_Y, rdiscretisation, nel, nptsPerEl, Dpoints, &
                              rdomainIntSubset, Dpy, rcollection)
          if (rprob%cformulation .eq. FORMULATION_MIXED) then
            Dcoefficients(1,:,:) = - rprob%dmu * (Du2xx + Du1yx) - 2*rprob%dmu * Du2yy &
                                   + 2*rprob%dmu * Dpy
            ! Why the factor 2 * mu in front of Dpx? See comment above!

          else ! FORMULATION_STOKES
            Dcoefficients(1,:,:) = - 2*rprob%dmu * (Du2xx + Du2yy) + Dpy
          endif
          deallocate(Dpy)
        endif
        deallocate(Du2xx, Du2yy, Du1yx)
      endif
    endif
    
  end subroutine elast_RHS_2D_vol


! ****************************************************************************************


!<subroutine>
  subroutine elast_RHS_2D_surf(rdiscretisation, rform, nel, nptsPerEl, &
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
    ! This subroutine is called during the vector assembly. It computes the function g
    ! in the linear form (usually L(v) = (f,v)_0 + (g,v)_N), i.e. the
    ! surface forces in x- or y-direction (nonzero Neumann boundary conditions).
!</description>

!<input>
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! linear form to be evaluated
    type(t_linearForm), intent(in) :: rform
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! DIMENSION(spatial dimension, nptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    
    ! index of the boundary component that contains the points
    integer, intent(in) :: ibct
    
    ! parameter values (LENGTH PARAMETRISATION) of the points on the boundary component
    ! DIMENSION(nptsPerEl, nel)
    real(DP), dimension(:,:), intent(in) :: DpointPar
    
    ! array of test space DOF in all elements
    ! DIMENSION(#local DOF in test space, nel)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! array of coefficients of the terms in the linear form for all given points on
    ! all given elements, DIMENSION(itermCount, nptsPerEl, nel)
    ! with itermCount the number of terms in the linear form
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), pointer :: DstressTensor
    real(DP) :: dminPar, dmaxPar, dt, dnx, dny
    integer :: iel, ipoint, icomp, iseg

    ! get the minimum and maximum parameter value, corresponding to start and end point
    ! of the interval
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, nel
      do ipoint = 1, nptsPerEl
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      enddo
    enddo

    ! in rcollection%IquickAccess(1) the current component is stored
    icomp = rcollection%IquickAccess(1)

    if (rprob%csimulation .eq. SIMUL_ANALYTICAL .and. icomp .lt. 3) then
      ! In case of an analytical test function, the Neumann contributions are given by
      ! sigma * n, where sigma is the Cauchy stress tensor and n the normal vector.
      ! For the third component (pressure) in a mixed formulation or Stokes, there is
      ! no such contribution.
      
      allocate(DstressTensor(2,nptsPerEl,nel))
      ! compute the required row icomp of the Cauchy stress tensor
      call elast_stressTensor(icomp, rdiscretisation, nel, nptsPerEl, Dpoints, &
                              rdomainIntSubset, DstressTensor, rcollection)
      do iel = 1, nel
        do ipoint = 1, nptsPerEl
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
          Dcoefficients(1,ipoint,iel) = dnx * DstressTensor(1,ipoint,iel) &
                                      + dny * DstressTensor(2,ipoint,iel)
        enddo
      enddo
      deallocate(DstressTensor)

    else if (rprob%csimulation .eq. SIMUL_REAL) then
      ! in case of a real simulation, the Neumann contributions are given by constant
      ! surface forces on the current segment, which are provided by the user

      ! in rcollection%IquickAccess(2) the current segment number is stored
      iseg = rcollection%IquickAccess(2)
      Dcoefficients(1,:,:) = rprob%DbcValue(icomp, iseg, ibct)
    else
      ! in every other case set the values to zero
      Dcoefficients(1,:,:) = 0.0_DP
    endif

  end subroutine elast_RHS_2D_surf


! ****************************************************************************************

  
!<subroutine>
  subroutine elast_stressTensor(irow, rdiscretisation, nel, nptsPerEl, Dpoints, &
                                rdomainIntSubset, Dstresses, rcollection)

    use fsystem
    use genoutput
    use collection
    use scalarpde
    use domainintegration
    use derivatives
    use elasticity_basic

!<description>
    ! This routine computes one row of the 2x2 Cauchy stress tensor in a given set
    ! of points/elements. The stress tensor is given by:
    !
    ! pure displacement formulation:
    !   2*mu*eps + lambda*div(u)*I
    !      / 2*mu*eps11 + lambda*div(u)       2*mu*eps12      \
    !   =  |                                                  |
    !      \      2*mu*eps12       2*mu*eps22 + lambda*div(u) /
    !
    !     / (2*mu + lambda)*u1x + lambda*u2y      mu*(u1y + u2x)       \
    !   = |                                                            |
    !     \        mu*(u1y + u2x)     (2*mu + lambda)*u2y + lambda*u1x /
    !
    ! mixed formulation:
    !   2*mu*eps - p*I
    !     / 2*mu*eps11 - p     2*mu*eps12    \   / 2*mu*u1x - p    mu*(u1y + u2x) \
    !   = |                                  | = |                                |
    !     \   2*mu*eps12      2*mu*eps22 - p /   \ mu*(u1y + u2x)   2*mu*u2y - p  /
    !
    ! Stokes:
    !   2*mu*grad(u) - p*I
    !     /  2*mu*grad(u)11 - p     2*mu*grad(u)12      \   / 2*mu*u1x - p   2*mu*u1y   \
    !   = |                                             | = |                           |
    !     \  2*mu*grad(u)21         2*mu*grad(u)22 - p  /   \   2*mu*u2x   2*mu*u2y - p /
!</description>
  
!<input>
    ! specifies which row of the stress tensor is to be computed
    integer, intent(in) :: irow
  
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! number of elements where the coefficients are to be computed
    integer, intent(in) :: nel
    
    ! number of points per element where the coefficients are to be computed
    integer, intent(in) :: nptsPerEl
    
    ! array of all points on all elements where coefficients are to be computed
    ! (usually coincides with rdomainSubset%p_DcubPtsReal)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
    ! structure providing more detailed information about the current element set
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
  
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! array for storing the desired stresses in all points on all elements
    !   DIMENSION(2, snptsPerEl, nel)
    real(DP), dimension(:,:,:), intent(out) :: Dstresses
!</output>
!</subroutine>

    real(DP), dimension(:,:), pointer :: Du1x, Du2x, Du1y, Du2y, Dpfunc
    allocate(Du1x(nptsPerEl,nel), Du2x(nptsPerEl,nel), &
             Du1y(nptsPerEl,nel), Du2y(nptsPerEl,nel))
    if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
        rprob%cformulation .eq. FORMULATION_STOKES) then
      allocate(Dpfunc(nptsPerEl,nel))
    endif
    call elast_analFunc(1, DER_DERIV_X,rdiscretisation, nel, nptsPerEl, Dpoints, &
                        rdomainIntSubset, Du1x, rcollection)
    call elast_analFunc(1, DER_DERIV_Y,rdiscretisation, nel, nptsPerEl, Dpoints, &
                        rdomainIntSubset, Du1y, rcollection)
    call elast_analFunc(2, DER_DERIV_X, rdiscretisation, nel, nptsPerEl, Dpoints, &
                        rdomainIntSubset, Du2x, rcollection)
    call elast_analFunc(2, DER_DERIV_Y, rdiscretisation, nel, nptsPerEl, Dpoints, &
                        rdomainIntSubset, Du2y, rcollection)                       
    if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
        rprob%cformulation .eq. FORMULATION_STOKES) then
      call elast_analFunc(3, DER_FUNC, rdiscretisation, nel, nptsPerEl, Dpoints, &
                          rdomainIntSubset, Dpfunc, rcollection)                       
    endif
    
    if (irow .eq. 1) then
      if (rprob%cformulation .eq. FORMULATION_DISPL) then
        Dstresses(1,:,:) = (2*rprob%dmu + rprob%dlambda) * Du1x +  rprob%dlambda * Du2y
        Dstresses(2,:,:) = rprob%dmu * (Du1y + Du2x)
      else if (rprob%cformulation .eq. FORMULATION_MIXED) then
        Dstresses(1,:,:) = 2*rprob%dmu * Du1x - Dpfunc
        Dstresses(2,:,:) = rprob%dmu * (Du1y + Du2x)
      else
        Dstresses(1,:,:) = 2*rprob%dmu * Du1x - Dpfunc
        Dstresses(2,:,:) = 2*rprob%dmu * Du1y
      endif
    else if (irow .eq. 2) then
      if (rprob%cformulation .eq. FORMULATION_DISPL) then
        Dstresses(1,:,:) = rprob%dmu * (Du2x + Du1y)
        Dstresses(2,:,:) = (2*rprob%dmu + rprob%dlambda) * Du2y +  rprob%dlambda * Du1x
      else if (rprob%cformulation .eq. FORMULATION_MIXED) then
        Dstresses(1,:,:) = rprob%dmu * (Du2x + Du1y)
        Dstresses(2,:,:) = 2*rprob%dmu * Du2y - Dpfunc
      else
        Dstresses(1,:,:) = 2*rprob%dmu * Du2x
        Dstresses(2,:,:) = 2*rprob%dmu * Du2y - Dpfunc
      endif
    else
      call output_line('Invalid row number ' // trim(sys_siL(irow,4)) // '!', &
                       OU_CLASS_ERROR, OU_MODE_STD, 'elast_stressTensor')
      call sys_halt()
    endif    

    deallocate(Du1x, Du2x, Du1y, Du2y)
    if (rprob%cformulation .eq. FORMULATION_MIXED .or. &
        rprob%cformulation .eq. FORMULATION_STOKES) then
      deallocate(Dpfunc)
    endif

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
    ! specifier for Dirichlet boundary
    ! Icomponents(1) defines the index of the boundary component the value is to be
    ! calculated for (e.g., 1 = 1st solution component = x-displacement, 2 = 2nd solution
    ! component = y-displacement)
    integer, dimension(:), intent(in) :: Icomponents
  
    ! discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! current boundary region
    type(t_boundaryRegion), intent(in) :: rboundaryRegion
    
    ! number of the current element on the boundary
    integer, intent(in) :: ielement
    
    ! The type of information, the routine should calculate. One of the DISCBC_NEEDxxxx
    ! constants. Depending on the constant, the routine has to return one or more
    ! information values in the result array.
    integer, intent(in) :: cinfoNeeded
    
    ! Reference to a geometric object where information is to be computed.
    ! cinfoNeeded = DISCBC_NEEDFUNC / DISCBC_NEEDDERIV:
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded = DISCBC_NEEDINTMEAN:
    !   iwhere = number of the edge where the integral mean value is to be computed
    integer, intent(in) :: iwhere
  
    ! Reference to a geometric object where information is to be computed.
    ! cinfoNeeded = DISCBC_NEEDFUNC / DISCBC_NEEDDERIV:
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded = DISCBC_NEEDINTMEAN:
    !   dwhere = 0 (not used)
    real(DP), intent(in) :: dwhere
!</input>
  
!<inputoutput>
    ! optional collection structure for additional information provided by the user
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
  
!<output>
    ! Array for storing the calculated values. If only one value is needed it is stored
    ! in Dvalues(1). If multiple values are needed, then, e.g., for DISCBC_NEEDDERIV:
    ! Dvalues(1) = x-derivative, Dvalues(2) = y-derivative, ...
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
      if (Icomponents(1) .lt. 3) then
        ! boundary values for x- or y-displacements
        Daux = elast_danalyticFunction(Dpoints, 1, 1, DER_FUNC, &
                                       rprob%CfuncID(Icomponents(1)))
      else
        ! boundary values for pressure solution
  
        if (rprob%dnu .eq. 0.5) then
          ! incompressible case
          Daux = elast_danalyticFunction(Dpoints, 1, 1, DER_FUNC, rprob%CfuncID(3))
        else
          ! in the compressible case the pressure is related to the displacements
          ! via p = - dlambda * div(u) = -dlambda * (u1_x + u2_y)
          Daux = -rprob%dlambda * &
            (elast_danalyticFunction(Dpoints, 1, 1, DER_DERIV_X, rprob%CfuncID(1)) &
           + elast_danalyticFunction(Dpoints, 1, 1, DER_DERIV_Y, rprob%CfuncID(2)))  
        endif
      endif
      Dvalues(1) = Daux(1,1)
    endif
  end subroutine elast_boundValue_2D

end module elasticity_callback


!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the Poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# --- 2D version ---
!#
!# 1.) coeff_Laplace_2D
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS_2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 2D case, Q2 bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) coeff_RHS_Sin2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 2D case, sinus bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) getBoundaryValues_2D
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 5.) getBoundaryValuesFBC_2D
!#     -> Returns analytic values in the inner of the domain on
!#        fictitious boundary objects
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcfassembly.inc'
!#
!# 6.) getBoundaryValuesMR_2D
!#     -> Returns discrete values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_discretebc.inc'
!#
!# 7.) getReferenceFunction_2D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_2D, Q2 bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 8.) getReferenceFunction_Sin2D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_Sin2D, sinus bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 9.) gethadaptMonitorFunction_2D
!#     -> Controls the grid adaption strategy in poisson2d_method1_hadapt.
!#
!# </purpose>
!##############################################################################

module poisson2d_callback

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

contains

! ***************************************************************************
  !<subroutine>

  subroutine coeff_Laplace_2D (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

   subroutine coeff_RHS_X_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form corresponding
    ! to the X-velocity.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients (1,:,:) = 0.0_DP
!     Dcoefficients (1,:,:) = -0.06_DP - &
!                      0.16_DP * SYS_PI**2 * sin(2.0_DP * SYS_PI * Dpoints(1,:,:)) &
!                      * sin(2.0_DP * SYS_PI * Dpoints(2,:,:))
    
!     Dcoefficients (1,:,:) = 2.2_DP * &
!                     Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
!                     0.6_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) - &
!                     0.16_DP * SYS_PI**2 * cos(2.0_DP * SYS_PI * Dpoints(1,:,:)) &
!                     * cos(2.0_DP * SYS_PI * Dpoints(2,:,:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

   subroutine coeff_RHS_Y_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form corresponding
    ! to the X-velocity.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients (1,:,:) = 0.0_DP
!     Dcoefficients (1,:,:) = 0.28_DP * &
!                      SYS_PI**2 * cos(2.0_DP * SYS_PI * Dpoints(1,:,:)) &
!                      * cos(2.0_DP * SYS_PI * Dpoints(2,:,:))

!    Dcoefficients (1,:,:) = - 0.8_DP * &
!                     (1.0_DP-2.0_DP*Dpoints(1,:,:)) * (1.0_DP-2.0_DP*Dpoints(2,:,:)) + &
!                     0.06_DP * SYS_PI**2 * sin(2.0_DP * SYS_PI * Dpoints(1,:,:)) * &
!                      sin(2.0_DP * SYS_PI * Dpoints(2,:,:)) + 0.22_DP * &
!                      SYS_PI**2 * sin(2.0_DP * SYS_PI * Dpoints(1,:,:)) * &
!                      sin(2.0_DP * SYS_PI * Dpoints(2,:,:))

  end subroutine

! ***************************************************************************

!<subroutine>

   subroutine coeff_RHS_neumBdr_X_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  ibct, DpointPar, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
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

   ! An array accepting the DOF's on all elements trial in the trial space.
   ! DIMENSION(#local DOF's in test space,nelements)
   integer, dimension(:,:), intent(in) :: IdofsTest

   ! This is a t_domainIntSubset structure specifying more detailed information
   ! about the element set that is currently being integrated.
   ! It's usually used in more complex situations (e.g. nonlinear matrices).
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
   real(DP) :: dminPar,dmaxPar,dt,dnx,dny,dnv
   integer :: icomp,iel,ipoint,ndim
   real(DP)                         :: fx,fy


   ! Get the minimum and maximum parameter value. The point with the minimal
   ! parameter value is the start point of the interval, the point with the
   ! maximum parameter value the endpoint.
   dminPar = DpointPar(1,1)
   dmaxPar = DpointPar(1,1)
   do iel = 1, size(rdomainIntSubset%p_Ielements)
     do ipoint = 1, ubound(Dpoints,2)
       dminPar = min(DpointPar(ipoint,iel), dminPar)
       dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
     end do
   end do

   ! Multiply the velocity vector with the normal in each point
   ! to get the normal velocity.
   do iel = 1, size(rdomainIntSubset%p_Ielements)
     do ipoint = 1, ubound(Dpoints,2)

       dt = DpointPar(ipoint,iel)

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
!     call f_x(ipoint, iel, fx)
!     call f_y(ipoint, iel, fy)
!     print *, dnx, dny, ibct, dt
! 
!        dnv = dnx * fx + dny * fy

       Dcoefficients(1,ipoint,iel) = 0.0_DP

     end do
   end do

 end subroutine coeff_RHS_neumBdr_X_2D


 ! ***************************************************************************

!<subroutine>

   subroutine coeff_RHS_neumBdr_Y_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  ibct, DpointPar, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
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

   ! An array accepting the DOF's on all elements trial in the trial space.
   ! DIMENSION(#local DOF's in test space,nelements)
   integer, dimension(:,:), intent(in) :: IdofsTest

   ! This is a t_domainIntSubset structure specifying more detailed information
   ! about the element set that is currently being integrated.
   ! It's usually used in more complex situations (e.g. nonlinear matrices).
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
   real(DP) :: dminPar,dmaxPar,dt,dnx,dny,dnv
   integer :: icomp,iel,ipoint,ndim
   real(DP)                         :: fx,fy


   ! Get the minimum and maximum parameter value. The point with the minimal
   ! parameter value is the start point of the interval, the point with the
   ! maximum parameter value the endpoint.
   dminPar = DpointPar(1,1)
   dmaxPar = DpointPar(1,1)
   do iel = 1, size(rdomainIntSubset%p_Ielements)
     do ipoint = 1, ubound(Dpoints,2)
       dminPar = min(DpointPar(ipoint,iel), dminPar)
       dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
     end do
   end do

   ! Multiply the velocity vector with the normal in each point
   ! to get the normal velocity.
   do iel = 1, size(rdomainIntSubset%p_Ielements)
     do ipoint = 1, ubound(Dpoints,2)

       dt = DpointPar(ipoint,iel)

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
!     call f_x(ipoint, iel, fx)
!     call f_y(ipoint, iel, fy)
    print *, dnx, dny, ibct, dt
! 
!        dnv = dnx * fx + dny * fy

       Dcoefficients(1,ipoint,iel) = -3.0e-2_DP

     end do
   end do

 end subroutine coeff_RHS_neumBdr_Y_2D


 ! ***************************************************************************

!<subroutine>

  subroutine  f_x(ipoint, iel, Dvalues)
  
!   use spatialdiscretisation
!   use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
!   integer, dimension(:), intent(IN)                           :: Icomponents

   ! The element number on the boundary which is currently being processed
!   integer, intent(IN)                                         :: ielement
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
!   real(DP), intent(IN)                                        :: dwhere

     integer, intent(IN) :: iel,ipoint

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), intent(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
!     integer :: icomponent
!      REAL(DP) :: dx,dy
!     
!      CALL boundary_getCoords(rdiscretisation%p_rboundary, &
!          rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Get from the current component of the PDE we are discretising:
!     icomponent = Icomponents(1)
    

     ! Return zero Dirichlet boundary values for all situations by default.
  Dvalues = 0.0_DP

  end subroutine

  ! ***************************************************************************
! <subroutine>

  subroutine  f_y(ipoint, iel, Dvalues)
  
!   use spatialdiscretisation
!   use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
!   integer, dimension(:), intent(IN)                           :: Icomponents

   ! The element number on the boundary which is currently being processed
!   integer, intent(IN)                                         :: ielement
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
!   real(DP), intent(IN)                                        :: dwhere

    integer, intent(IN) :: iel,ipoint

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), intent(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
!     integer :: icomponent
!      REAL(DP) :: dx,dy
!     
!      CALL boundary_getCoords(rdiscretisation%p_rboundary, &
!          rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Get from the current component of the PDE we are discretising:
!     icomponent = Icomponents(1)
    

     ! Return zero Dirichlet boundary values for all situations by default.
  Dvalues = -3.0_DP

  end subroutine
  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_X_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
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
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

!     select case (cderivative)
!     case (DER_FUNC)
! !     u(x,y) = x*(1-x)*y*(1-y)
!     Dvalues (:,:) = Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) * &
!                              Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
!   case (DER_DERIV_X)
! !        u(x,y)   = x*(1-x)*y*(1-y)
! !     => u_x(x,y) =  y*(1-x)*(1-y)-x*y*(1-y) 
!    Dvalues (:,:) = Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
!        Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(2,:,:))
!   case (DER_DERIV_Y)
! !        u(x,y)   = x*(1-x)*y*(1-y)
! !     => u_y(x,y) = x*(1-x)*(1-y)-x*y*(1-x)
!    Dvalues (:,:) = Dpoints(1,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
!        Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:))
!  case DEFAULT
! !     Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = 0.1*y**2
    Dvalues (:,:) = 0.1_DP * Dpoints(2,:,:) * Dpoints(2,:,:)
  case (DER_DERIV_X)
    !    u(x,y) = 0.1 * y**2
    ! => u_x(x,y) =  0 
    Dvalues (:,:) = 0.0_DP
  case (DER_DERIV_Y)
    !    u(x,y) = 0.1 * y**2
    ! => u_y(x,y) = 0.2 * y
    Dvalues (:,:) = 0.2_DP * Dpoints(2,:,:)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
!   

  end subroutine


  ! ***************************************************************************
!<subroutine>

  subroutine getReferenceFunction_Y_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
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
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = 0.05*COS(2*PI * x) * COS(2*PI * y)
    Dvalues (:,:) = 0.05_DP * cos(2.0_DP * SYS_PI*Dpoints(1,:,:)) * cos(2.0_DP * SYS_PI*Dpoints(2,:,:))
  case (DER_DERIV_X)
    !    u(x,y)   = 0.05*COS(2*PI * x) * COS(2*PI * y)
    ! => u_x(x,y) = -0.1 * PI * SIN(2* PI * x) * COS(2* PI * y)
    Dvalues (:,:) = -0.1_DP * SYS_PI * sin(2.0_DP * SYS_PI*Dpoints(1,:,:)) * cos(2.0_DP * SYS_PI*Dpoints(2,:,:))
  case (DER_DERIV_Y)
    !    u(x,y)   = 0.05*COS(2*PI * x) * COS(2*PI * y)
    ! => u_y(x,y) = -0.1 * PI * COS(2 * PI * x) * SIN(2 * PI * y)
    Dvalues (:,:) = -0.1_DP * SYS_PI * cos(2.0_DP * SYS_PI*Dpoints(1,:,:)) * sin(2.0_DP * SYS_PI*Dpoints(2,:,:))
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select

 
!   select case (cderivative)
!   case (DER_FUNC)
!     ! u(x,y) = 0.05*SIN(2*PI * x) * SIN(2*PI * y)
!     Dvalues (:,:) = 0.05_DP * sin(2.0_DP * SYS_PI*Dpoints(1,:,:)) * sin(2.0_DP * SYS_PI*Dpoints(2,:,:))
!   case (DER_DERIV_X)
!     !    u(x,y)   = 0.05*SIN(2 * PI * x) * SIN(2 * PI * y)
!     ! => u_x(x,y) = 0.1 * PI * COS(2* PI * x) * SIN(2* PI * y)
!     Dvalues (:,:) = 0.1_DP * SYS_PI * cos(2.0_DP * SYS_PI*Dpoints(1,:,:)) * sin(2.0_DP * SYS_PI*Dpoints(2,:,:))
!   case (DER_DERIV_Y)
!     !    u(x,y)   = 0.05*SIN(2 * PI * x) * SIN(2 * PI * y)
!     ! => u_y(x,y) = 0.1 * PI * SIN(2 * PI * x) * COS(2 * PI * y)
!     Dvalues (:,:) = 0.1_DP * SYS_PI * sin(2.0_DP * SYS_PI*Dpoints(1,:,:)) * cos(2.0_DP * SYS_PI*Dpoints(2,:,:))
!   case DEFAULT
!     ! Unknown. Set the result to 0.0.
!     Dvalues = 0.0_DP
!   end select
  

  end subroutine



  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere, dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(IN)                                         :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(IN)                                         :: cinfoNeeded
  
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
  integer, intent(IN)                                          :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  real(DP), intent(IN)                                        :: dwhere

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
!     integer :: icomponent
!      REAL(DP) :: dx,dy
!     
!      CALL boundary_getCoords(rdiscretisation%p_rboundary, &
!          rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Get from the current component of the PDE we are discretising:
!     icomponent = Icomponents(1)
    

     ! Return zero Dirichlet boundary values for all situations by default.
  Dvalues(1) = 0.0_DP

  ! Now, depending on the problem, calculate the actual velocity value.

!   select case (icomponent)
!   case (1) ! X-velocity
!     Dvalues(1) = 0.1_DP * dy * dy
! 
!    case (2) ! Y-velocity
!     Dvalues(1) = 0.05_DP *cos(2.0_DP*SYS_PI*dx)* cos(2.0_DP*SYS_PI*dy)
! 
!    end select


  end subroutine


end module

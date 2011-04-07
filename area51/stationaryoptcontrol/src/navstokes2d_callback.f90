!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_callback </name>
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

module navstokes2d_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use feevaluation
  use derivatives
  use spatialdiscretisation
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use feevaluation
  
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
    integer, intent(IN)                        :: nelements
    
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

  subroutine coeff_RHS_primal_2D (rdiscretisation,rform, &
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
    ! the coefficients in front of the terms of the linear form.
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
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    !    u(x,y) = 16*x*(1-x)*y*(1-y)
    ! => f(x,y) = 32 * (y*(1-y)+x*(1-x))
    !Dcoefficients (1,:,:) = 32.0_DP * &
    !                ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
    !                  Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )
    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSx_dual_2D (rdiscretisation,rform, &
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
    ! the coefficients in front of the terms of the linear form.
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
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer, dimension(2) :: Isize
    integer :: iel
    real(DP), dimension(:,:), allocatable :: Dvalues
    real(DP), dimension(:), allocatable :: Dvalues2

    ! RHS of the dual equation is the negative target function:
    !
    !    u(x,y) = 16*x*(1-x)*y*(1-y)
    ! => z(x,y) = -u(x,y)
    if (rcollection%IquickAccess(1) .eq. 0) then
      ! Analytically given target flow
      Isize(1) = Ubound(Dcoefficients,2)
      Isize(2) = Ubound(Dcoefficients,3)
      allocate (Dvalues(Isize(1),Isize(2)))
      call getReferenceFunctionX_2D (DER_FUNC,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,Dvalues,rcollection)
      Dcoefficients (1,:,:) = -Dvalues(:,:)
      deallocate (Dvalues)
    else if (rcollection%IquickAccess(1) .eq. 1) then
      ! Target flow read from file
      allocate (Dvalues2(Ubound(Dcoefficients,2)))
      do iel = 1,nelements
        call fevl_evaluate (DER_FUNC, Dvalues2, &
            rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints(:,:,iel))
        Dcoefficients (1,:,iel) = -Dvalues2(:)
      end do
      deallocate(Dvalues2)
    else if (rcollection%IquickAccess(1) .eq. 2) then
      Dcoefficients (1,:,:) = 0.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSy_dual_2D (rdiscretisation,rform, &
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
    ! the coefficients in front of the terms of the linear form.
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
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer, dimension(2) :: Isize
    integer :: iel
    real(DP), dimension(:,:), allocatable :: Dvalues
    real(DP), dimension(:), allocatable :: Dvalues2

    ! RHS of the dual equation is the negative target function:
    !
    !    u(x,y) = 16*x*(1-x)*y*(1-y)
    ! => z(x,y) = -u(x,y)
    if (rcollection%IquickAccess(1) .eq. 0) then
      ! Analytically given target flow
      Isize(1) = Ubound(Dcoefficients,2)
      Isize(2) = Ubound(Dcoefficients,3)
      allocate (Dvalues(Isize(1),Isize(2)))
      call getReferenceFunctionY_2D (DER_FUNC,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,Dvalues,rcollection)
      Dcoefficients (1,:,:) = -Dvalues(:,:)
      deallocate (Dvalues)
    else if (rcollection%IquickAccess(1) .eq. 1) then
      ! Target flow read from file
      allocate (Dvalues2(Ubound(Dcoefficients,2)))
      do iel = 1,nelements
        call fevl_evaluate (DER_FUNC, Dvalues2, &
            rcollection%p_rvectorQuickAccess1%RvectorBlock(2), Dpoints(:,:,iel))
        Dcoefficients (1,:,iel) = -Dvalues2(:)
      end do
      deallocate(Dvalues2)
    else if (rcollection%IquickAccess(1) .eq. 2) then
      Dcoefficients (1,:,:) = 0.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunctionX_2D (cderivative,rdiscretisation, &
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

    integer :: iel

    if (rcollection%IquickAccess(1) .eq. 0) then
      ! Analytical function
      select case (cderivative)
      case (DER_FUNC)
        ! u(x,y) = 16*x*(1-x)*y*(1-y)
        Dvalues (:,:) = 0.3_DP * 4.0_DP * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
        !Dvalues (:,:) = 0.0_DP
        !Dvalues (:,:) = 2.0_DP/3.0_DP
        !Dvalues (:,:) = 4.0_DP * Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) * (1.0_DP-Dpoints(1,:,:)) + &
        !  2.0_DP/3.0_DP * Dpoints(1,:,:)
      case (DER_DERIV_X)
        !    u(x,y)   = 16*x*(1-x)*y*(1-y)
        ! => u_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
        Dvalues (:,:) = 16.0_DP * ( &
            Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
            Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(2,:,:)) )
      case (DER_DERIV_Y)
        !    u(x,y)   = 16*x*(1-x)*y*(1-y)
        ! => u_y(x,y) = 16 * ( x*(1-x)*(1-y)-x*y*(1-x) )
        Dvalues (:,:) = 16.0_DP * ( &
            Dpoints(1,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
            Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) )
      case DEFAULT
        ! Unknown. Set the result to 0.0.
        Dvalues = 0.0_DP
      end select
    else if (rcollection%IquickAccess(1) .eq. 1) then
      ! Function given by vector
      ! Target flow read from file
      do iel = 1,nelements
        call fevl_evaluate (DER_FUNC, Dvalues(:,iel), &
            rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints(:,:,iel))
      end do
    else if (rcollection%IquickAccess(1) .eq. 2) then
      Dvalues(:,:) = 0.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunctionY_2D (cderivative,rdiscretisation, &
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

    integer :: iel

    if (rcollection%IquickAccess(1) .eq. 0) then
      select case (cderivative)
      case (DER_FUNC)
        ! u(x,y) = 16*x*(1-x)*y*(1-y)
        Dvalues (:,:) = 0.0_DP
      case (DER_DERIV_X)
        !    u(x,y)   = 16*x*(1-x)*y*(1-y)
        ! => u_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
        Dvalues (:,:) = 16.0_DP * ( &
            Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
            Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(2,:,:)) )
      case (DER_DERIV_Y)
        !    u(x,y)   = 16*x*(1-x)*y*(1-y)
        ! => u_y(x,y) = 16 * ( x*(1-x)*(1-y)-x*y*(1-x) )
        Dvalues (:,:) = 16.0_DP * ( &
            Dpoints(1,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
            Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) )
      case DEFAULT
        ! Unknown. Set the result to 0.0.
        Dvalues = 0.0_DP
      end select
    else if (rcollection%IquickAccess(1) .eq. 1) then
      ! Function given by vector
      ! Target flow read from file
      do iel = 1,nelements
        call fevl_evaluate (DER_FUNC, Dvalues(:,iel), &
            rcollection%p_rvectorQuickAccess1%RvectorBlock(2), Dpoints(:,:,iel))
      end do
    else if (rcollection%IquickAccess(1) .eq. 2) then
      Dvalues(:,:) = 0.0_DP
    end if

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine getErrorX_2D (cderivative,rdiscretisation, &
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
  ! the error y-z in a couple of points on a couple
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

    integer, dimension(3) :: Isize
    real(dp), dimension(:,:,:), allocatable :: Ddata

    ! Calculate z
    Isize(1:2) = ubound(Dvalues)
    Isize(3) = 2
    allocate(Ddata(Isize(1),Isize(2),Isize(3)))
    
    call getReferenceFunctionX_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Ddata(:,:,1),rcollection)
    call fevl_evaluate_sim1 (DER_FUNC, Ddata(:,:,2), &
        rcollection%p_rvectorQuickAccess2%RvectorBlock(1), Dpoints, &
        rdomainIntSubset%p_Ielements)
    Dvalues = Ddata(:,:,2)-Ddata(:,:,1)
    deallocate(Ddata)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getErrorY_2D (cderivative,rdiscretisation, &
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
  ! the error y-z in a couple of points on a couple
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

    integer, dimension(3) :: Isize
    real(dp), dimension(:,:,:), allocatable :: Ddata

    ! Calculate z
    Isize(1:2) = ubound(Dvalues)
    Isize(3) = 2
    allocate(Ddata(Isize(1),Isize(2),Isize(3)))
    
    call getReferenceFunctionY_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Ddata(:,:,1),rcollection)
    call fevl_evaluate_sim1 (DER_FUNC, Ddata(:,:,2), &
        rcollection%p_rvectorQuickAccess2%RvectorBlock(2), Dpoints, &
        rdomainIntSubset%p_Ielements)
    Dvalues = Ddata(:,:,2)-Ddata(:,:,1)
    deallocate(Ddata)

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValuesX_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
    
    use fsystem
    use boundary
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
    !   Icomponents(1) defines the number of the solution component, the value
    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
    !   2=2nd solution component, e.g. Y-velocity,...,
    !   3=3rd solution component, e.g. pressure)
    ! For pressure drop boundary / normal stress:
    !   Velocity components that are affected by the normal stress
    !   (usually "1 2" for x- and y-velocity while returned value musr specify
    !   the pressure at the boundary)
    integer, dimension(:), intent(in)                           :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
    
    ! The element number on the boundary which is currently being processed
    integer, intent(in)                                         :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in)                                         :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDFUNCMID : 
    !   iwhere = number of the edge in which midpoint the value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   iwhere = Number of the edge where the normal stress should be computed.
    integer, intent(in)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   dwhere = parameter value of the point on edge iwhere where the normal
    !            stress should be computed.
    real(DP), intent(in)                                        :: dwhere
     
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    !
    ! The function may return SYS_INFINITY_DP as a value. This indicates the
    ! framework to ignore the node and treat it as 'natural boundary condition'
    ! node.
    real(DP), dimension(:), intent(out)                         :: Dvalues
  !</output>
    
  !</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    REAL(DP) :: dx,dy
    
    CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.3_DP*4.0_DP*dy*(1.0_DP-dy)
  
  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues1_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
    
    use fsystem
    use boundary
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
    !   Icomponents(1) defines the number of the solution component, the value
    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
    !   2=2nd solution component, e.g. Y-velocity,...,
    !   3=3rd solution component, e.g. pressure)
    ! For pressure drop boundary / normal stress:
    !   Velocity components that are affected by the normal stress
    !   (usually "1 2" for x- and y-velocity while returned value musr specify
    !   the pressure at the boundary)
    integer, dimension(:), intent(in)                           :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
    
    ! The element number on the boundary which is currently being processed
    integer, intent(in)                                         :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in)                                         :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDFUNCMID : 
    !   iwhere = number of the edge in which midpoint the value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   iwhere = Number of the edge where the normal stress should be computed.
    integer, intent(in)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   dwhere = parameter value of the point on edge iwhere where the normal
    !            stress should be computed.
    real(DP), intent(in)                                        :: dwhere
     
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    !
    ! The function may return SYS_INFINITY_DP as a value. This indicates the
    ! framework to ignore the node and treat it as 'natural boundary condition'
    ! node.
    real(DP), dimension(:), intent(out)                         :: Dvalues
  !</output>
    
  !</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! REAL(DP) :: dx,dy
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 1.0_DP
  
  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValuesY_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
    
    use fsystem
    use boundary
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
    !   Icomponents(1) defines the number of the solution component, the value
    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
    !   2=2nd solution component, e.g. Y-velocity,...,
    !   3=3rd solution component, e.g. pressure)
    ! For pressure drop boundary / normal stress:
    !   Velocity components that are affected by the normal stress
    !   (usually "1 2" for x- and y-velocity while returned value musr specify
    !   the pressure at the boundary)
    integer, dimension(:), intent(in)                           :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
    
    ! The element number on the boundary which is currently being processed
    integer, intent(in)                                         :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in)                                         :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDFUNCMID : 
    !   iwhere = number of the edge in which midpoint the value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   iwhere = Number of the edge where the normal stress should be computed.
    integer, intent(in)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   dwhere = parameter value of the point on edge iwhere where the normal
    !            stress should be computed.
    real(DP), intent(in)                                        :: dwhere
     
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    !
    ! The function may return SYS_INFINITY_DP as a value. This indicates the
    ! framework to ignore the node and treat it as 'natural boundary condition'
    ! node.
    real(DP), dimension(:), intent(out)                         :: Dvalues
  !</output>
    
  !</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    REAL(DP) :: dx,dy
    
    CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.1_DP*4.0_DP*dy*(1.0_DP-dy)
    !Dvalues(1) = 0.0_DP
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValuesZero_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
    
    use fsystem
    use boundary
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
    !   Icomponents(1) defines the number of the solution component, the value
    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry, 
    !   2=2nd solution component, e.g. Y-velocity,...,
    !   3=3rd solution component, e.g. pressure)
    ! For pressure drop boundary / normal stress:
    !   Velocity components that are affected by the normal stress
    !   (usually "1 2" for x- and y-velocity while returned value musr specify
    !   the pressure at the boundary)
    integer, dimension(:), intent(in)                           :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
    
    ! The element number on the boundary which is currently being processed
    integer, intent(in)                                         :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in)                                         :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDFUNCMID : 
    !   iwhere = number of the edge in which midpoint the value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   iwhere = Number of the edge where the normal stress should be computed.
    integer, intent(in)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS : 
    !   dwhere = parameter value of the point on edge iwhere where the normal
    !            stress should be computed.
    real(DP), intent(in)                                        :: dwhere
     
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    !
    ! The function may return SYS_INFINITY_DP as a value. This indicates the
    ! framework to ignore the node and treat it as 'natural boundary condition'
    ! node.
    real(DP), dimension(:), intent(out)                         :: Dvalues
  !</output>
    
  !</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! REAL(DP) :: dx,dy
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP
  
  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMass (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the projective mass matrix.
    !
    ! The coefficients is c=c(lambda_1) or =c(lambda_2) with
    ! c=1/alpha if a < -1/alpha lambda_i < b and c=0 otherwise. This is the derivative
    ! of the projection operator "-P[a,b](-1/alpha lambda_i)" on the left hand
    ! side of the equation.
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
    integer, intent(IN)                        :: nelements
    
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
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    type(t_vectorScalar), pointer :: p_rsubvector
    real(dp), dimension(:,:), allocatable :: Dfunc
    integer(I32) :: celement
    real(DP) :: da, db, dalpha, dp1, dval
    integer :: ipt, iel
    
    ! Get the bounds and the multiplier from the collection
    dalpha = rcollection%DquickAccess(3)
    dp1 = 1.0_DP/dalpha
    
    ! Scale the bounds by -alpha as we analyse lambda
    ! and not u. Change the role of min/max because of the "-"
    ! sign!
    da = rcollection%DquickAccess(1)
    db = rcollection%DquickAccess(2)
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Do we have to analyse lambda_1 or lambda_2?
    if (rcollection%IquickAccess(1) .eq. 1) then
      p_rsubvector => p_rvector%RvectorBlock(4)
    else
      p_rsubvector => p_rvector%RvectorBlock(5)
    end if
  
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    call fevl_evaluate_sim (p_rsubvector, &
        rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
    
    ! Now check the function values lambda.
    ! If a < -1/alpha lambda < b, return 1/alpha.
    ! Otherwise, return 0.
    do iel = 1,ubound(Dcoefficients,3)
!      dp1 = 0.0_dp
!      do ipt = 1,ubound(dcoefficients,2)
!        ! check if the dual variable is in the bounds for the control.
!        dval = -dcoefficients(1,ipt,iel)/dalpha
!        if ((dval .gt. da) .and. (dval .lt. db)) then
!          dp1 = dp1 + 1.0_dp/dalpha
!        end if
!      end do
!      dp1 = dp1 * 0.25_dp
      
      do ipt = 1,ubound(Dcoefficients,2)
        ! Check if the dual variable is in the bounds for the control.
        dval = -Dcoefficients(1,ipt,iel)/dalpha
        if ((dval .gt. da) .and. (dval .lt. db)) then
          Dcoefficients(1,ipt,iel) = dp1
        else
          Dcoefficients(1,ipt,iel) = 0.0_DP
        end if
      end do
    end do
    
  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> sse_callback_poisson </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the SSE problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in "intf_xxxx.inc" files.
!#
!#
!# 1.) coeff_MatrixA_Poisson
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 2.) coeff_MatrixA_Bdr_Poisson
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 3.) coeff_RHS_Poisson
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 4.) coeff_RHS_Bdr_Poisson
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 5.) getBoundaryValues_Poisson
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#
!# 6.) getReferenceFunction_Poisson
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 7.) getReferenceDerivX_Poisson
!#     -> Returns the values of the analytic derivatives in x-direction.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the recovered FE gradient in comparison to the
!#        analytic derivative
!#
!# 8.) getReferenceDerivY_Poisson
!#      -> Returns the values of the analytic derivatives in x-direction.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#         $H_1$-error of the recovered FE gradient in comparison to the
!#         analytic derivative
!#
!# 9.) getReferenceDerivXX_Poisson
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#         $H_1$-error of the FE function in comparison to the analytic
!#         function
!#
!# 10.) getReferenceDerivXY_Poisson
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 11.) getReferenceDerivYX_Poisson
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 12.) getReferenceDerivYY_Poisson
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 13.) getAnalyticValues_Poisson
!#      -> Returns the values of the analytic function and its derivatives.
!#
!# </purpose>
!##############################################################################

module sse_callback_poisson

  use fsystem
  use storage
  use genoutput
  use derivatives
  use boundary
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use element
  use cubature
  use spatialdiscretisation
  use scalarpde
  use domainintegration
  use collection
  use discretebc
  use discretefbc
  use pprocgradients
  use pprocerror
  
  use sse_base

  implicit none

  private

  public :: coeff_MatrixA_Poisson
  public :: coeff_MatrixA_Bdr_Poisson
  public :: coeff_RHS_Poisson
  public :: coeff_RHS_Bdr_Poisson
  public :: getBoundaryValues_Poisson
  public :: getReferenceFunction_Poisson
  public :: getReferenceDerivX_Poisson
  public :: getReferenceDerivY_Poisson
  public :: getReferenceDerivXX_Poisson
  public :: getReferenceDerivXY_Poisson
  public :: getReferenceDerivYX_Poisson
  public :: getReferenceDerivYY_Poisson
  public :: getAnalyticValues_Poisson

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Poisson(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,IdofsTrial,IdofsTest,&
      rdomainIntSubset,Dcoefficients,rcollection)
    
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
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
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

    Dcoefficients(:,:,:) = dpoisson

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Bdr_Poisson(rdiscretisationTrial,rdiscretisationTest,&
      rform,nelements,npointsPerElement,Dpoints,ibct,DpointPar,IdofsTrial,&
      IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use triangulation
    use spatialdiscretisation, only: t_spatialDiscretisation

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
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>

  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>

  !</subroutine>

    Dcoefficients = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Poisson(rdiscretisation,rform, &
      nelements,npointsPerElement,Dpoints, &
      IdofsTest,rdomainIntSubset,&
      Dcoefficients,rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

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

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
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

    !    u(x,y) = SIN(PI * x) * SIN(PI * y)
    ! => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)
    Dcoefficients (1,:,:) = 2.0_DP * SYS_PI * SYS_PI &
                          * sin(SYS_PI * Dpoints(1,:,:)) &
                          * sin(SYS_PI * Dpoints(2,:,:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_Poisson(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

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

#if defined(CASE_POISSON_DIRICHLET)

    call output_line("There is no boundary integral in this benc", &
        OU_CLASS_ERROR,OU_MODE_STD,"coeff_RHS_Bdr_Poisson")
    call sys_halt()

#elif defined(CASE_POISSON_NEUMANN)

    Dcoefficients (1,:,:) = -ddirichlet

#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_Poisson(Icomponents,rdiscretisation,rboundaryRegion,&
      ielement,cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! "snapshot" of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary:
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
  !   2=2nd solution component, e.g. Y-velocity,...)
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

#if defined(CASE_POISSON_DIRICHLET)

    select case(rcollection%IquickAccess(1))
    case(1)
      
      ! Return Dirichlet boundary values for all situations.
      Dvalues(1) = ddirichlet
      
    case default
      call output_line("There is no boundary integral in this benchmark", &
          OU_CLASS_ERROR,OU_MODE_STD,"getBoundaryValues_Poisson")
      call sys_halt()
    end select

#elif defined(CASE_POISSON_NEUMANN)

    ! local variables
    real(DP) :: dnx,dny

    select case(rcollection%IquickAccess(1))
    case(1)

      ! Return Dirichlet boundary values for all situations.
      Dvalues(1) = ddirichlet

    case(2)

      ! Compute the normal vector in the point on the boundary
      call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

      ! Return Neumann boundary values for all situations.
      Dvalues(1) = dneumann*dnx

    case(3)

      ! Compute the normal vector in the point on the boundary
      call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

      ! Return Neumann boundary values for all situations.
      Dvalues(1) = dneumann*dny

    case default
      call output_line("There is no boundary integral in this benchmark", &
          OU_CLASS_ERROR,OU_MODE_STD,"getBoundaryValues_Poisson")
      call sys_halt()
    end select

#else
#error 'Test case is undefined.'
#endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivX_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xx(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xy(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivY_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yy(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yx(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXX_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xx(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xxx(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xxy(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXY_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_xy(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xyx(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_xyy(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYX_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yx(x,y)= PI * PI * COS(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yxx(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yxy(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYY_Poisson(cderivative,rdiscretisation, &
      nelements,npointsPerElement,Dpoints,IdofsTest,rdomainIntSubset,&
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

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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

  select case (cderivative)
  case (DER_FUNC)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_yy(x,y)= -PI * PI * SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_X)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yyx(x,y)= -PI * PI * PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))

  case (DER_DERIV_Y)
    !    u(x,y)    = SIN(PI * x) * SIN(PI * y)
    ! => u_yyy(x,y)= -PI * PI * PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = -SYS_PI*SYS_PI*SYS_PI * &
        sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  elemental subroutine getAnalyticValues_Poisson(dx,dy,&
      du,du_x,du_y,du_xx,du_xy,du_yy)

!<description>
    ! This function computes the values of the analytic function and
    ! its derivatives.
!</description>

!<input>
    ! Coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<output>
    ! Solution value
    real(DP), intent(out) :: du

    ! Values of the first derivatives
    real(DP), intent(out) :: du_x,du_y

    ! Values of the second derivatives
    real(DP), intent(out) :: du_xx,du_xy,du_yy
!</output>
!</subroutine>

    ! Solution values
    du = sin(SYS_PI*dx) * sin(SYS_PI*dy)

    ! Values of first derivatives
    du_x = SYS_PI * cos(SYS_PI*dx) * sin(SYS_PI*dy)
    du_y = SYS_PI * sin(SYS_PI*dx) * cos(SYS_PI*dy)

    ! Values of second derivatives
    du_xx = -SYS_PI*SYS_PI * sin(SYS_PI*dx) * sin(SYS_PI*dy)
    du_xy =  SYS_PI*SYS_PI * cos(SYS_PI*dx) * cos(SYS_PI*dy)
    du_yy = -SYS_PI*SYS_PI * sin(SYS_PI*dx) * sin(SYS_PI*dy)

  end subroutine getAnalyticValues_Poisson

end module sse_callback_poisson

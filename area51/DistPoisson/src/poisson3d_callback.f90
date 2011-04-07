!##############################################################################
!# ****************************************************************************
!# <name> poisson3d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the Poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# --- 3D version ---
!#
!# 1.) coeff_Laplace_3D
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS_3D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 3D case, Q2 bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) coeff_RHS_Sin3D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 3D case, sinus bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) getBoundaryValuesMR_3D
!#     -> Returns discrete values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_discretebc.inc'
!#
!# 5.) getReferenceFunction_3D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_3D, Q2 bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 6.) getReferenceFunction_Sin3D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_3D, sinus bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# </purpose>
!##############################################################################

module poisson3d_callback

  use fsystem
  use storage
  use linearsolver
  use boundary
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  
  implicit none

contains

! ***************************************************************************
  !<subroutine>

  subroutine coeff_Laplace_3D(rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)    :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)    :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)             :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                          :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                          :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)       :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in)          :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in)          :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)          :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional  :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)      :: Dcoefficients
  !</output>
    
  !</subroutine>

      Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_3D (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)    :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)               :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                          :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                          :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)       :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in)          :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)          :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional  :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

      !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
      ! => f(x,y,z) = 128 * ( y*(1-y)*z*(1-z)
      !             + x*(1-x)*z*(1-z) + x*(1-x)*y*(1-y))
      Dcoefficients(1,:,:) = 1.0_dp
          

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_3D (cderivative,rdiscretisation, &
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
  integer, intent(in)                            :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)      :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                            :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                            :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(in)         :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in)            :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)            :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional    :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)          :: Dvalues
!</output>
  
!</subroutine>

    select case (cderivative)
    case (DER_FUNC3D)
      !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
      Dvalues(:,:) = 64.0_DP * Dpoints(1,:,:) * (1.0_DP - Dpoints(1,:,:))&
                             * Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))&
                             * Dpoints(3,:,:) * (1.0_DP - Dpoints(3,:,:))
    case (DER_DERIV3D_X)
      !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
      ! => u_x(x,y,z) = 64*(1-2*x)*y*(1-y)*z*(1-z)
      Dvalues(:,:) = 64.0_DP * (1.0_DP - 2.0_DP * Dpoints(1,:,:))&
                             * Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))&
                             * Dpoints(3,:,:) * (1.0_DP - Dpoints(3,:,:))
    case (DER_DERIV3D_Y)
      !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
      ! => u_y(x,y,z) = 64*(1-2*y)*x*(1-x)*z*(1-z)
      Dvalues(:,:) = 64.0_DP * (1.0_DP - 2.0_DP * Dpoints(2,:,:))&
                             * Dpoints(1,:,:) * (1.0_DP - Dpoints(1,:,:))&
                             * Dpoints(3,:,:) * (1.0_DP - Dpoints(3,:,:))
    case (DER_DERIV3D_Z)
      !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
      ! => u_y(x,y,z) = 64*(1-2*z)*x*(1-x)*y*(1-y)
      Dvalues(:,:) = 64.0_DP * (1.0_DP - 2.0_DP * Dpoints(3,:,:))&
                             * Dpoints(1,:,:) * (1.0_DP - Dpoints(1,:,:))&
                             * Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))
    case default
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Sin3D (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)    :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)               :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                          :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                          :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)       :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in)          :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)          :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional  :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

      !    u(x,y,z) = SIN(PI * x) * SIN(PI * y) * SIN(PI * z)
      ! => f(x,y,z) = 3 * PI^2 * SIN(PI * x) * SIN(PI * y) * SIN(PI * z)
      Dcoefficients(1,:,:) = 3.0_DP * SYS_PI**2 * sin(SYS_PI*Dpoints(1,:,:))*&
                       sin(SYS_PI*Dpoints(2,:,:))*sin(SYS_PI*Dpoints(3,:,:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_Sin3D (cderivative,rdiscretisation, &
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
  integer, intent(in)                            :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)      :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                            :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                            :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(in)         :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in)            :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)            :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional    :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)          :: Dvalues
!</output>
  
!</subroutine>

    select case (cderivative)
    case (DER_FUNC3D)
      !    u(x,y,z) = SIN(PI * x) * SIN(PI * y) * SIN(PI * z)
      Dvalues(:,:) = sin(SYS_PI*Dpoints(1,:,:))* &
                     sin(SYS_PI*Dpoints(2,:,:))*sin(SYS_PI*Dpoints(3,:,:))
    case (DER_DERIV3D_X)
      !    u(x,y,z)   = SIN(PI * x) * SIN(PI * y) * SIN(PI * z)
      ! => u_x(x,y,z) = PI * COS(PI * x) * SIN(PI * y) * SIN(PI * z)
      Dvalues(:,:) = SYS_PI * cos(SYS_PI*Dpoints(1,:,:))* &
                     sin(SYS_PI*Dpoints(2,:,:))*sin(SYS_PI*Dpoints(3,:,:))
    case (DER_DERIV3D_Y)
      !    u(x,y,z)   = SIN(PI * x) * SIN(PI * y) * SIN(PI * z)
      ! => u_y(x,y,z) = PI * SIN(PI * x) * COS(PI * y) * SIN(PI * z)
      Dvalues(:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:))* &
                     cos(SYS_PI*Dpoints(2,:,:))*sin(SYS_PI*Dpoints(3,:,:))
    case (DER_DERIV3D_Z)
      !    u(x,y,z)   = SIN(PI * x) * SIN(PI * y) * SIN(PI * z)
      ! => u_z(x,y,z) = PI * SIN(PI * x) * SIN(PI * y) * COS(PI * z)
      Dvalues(:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:))* &
                     sin(SYS_PI*Dpoints(2,:,:))*cos(SYS_PI*Dpoints(3,:,:))
    case default
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select
  
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine getBoundaryValuesMR_3D (Icomponents,rdiscretisation,rmeshRegion,&
                                      cinfoNeeded,Iwhere,Dwhere,Dcoords,Dvalues,&
                                      rcollection)
  
  use collection
  use spatialdiscretisation
  use meshregion
  
!<description>
  ! This subroutine is called during the assembly of boundary conditions which
  ! are defined on mesh regions.
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
  integer, dimension(:), intent(in)              :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)      :: rdiscretisation
  
  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in)                 :: rmeshRegion

  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                            :: cinfoNeeded
  
  ! An array holding information about what type of DOF is currently processed.
  ! The information is build up as follows:
  ! Iwhere(1) = vertice number of the DOF, if the DOF is vertice-based, otherwise 0
  ! Iwhere(2) = edge number of the DOF, if the DOF is edge-based, otherwise 0
  ! Iwhere(3) = face number of the DOF, if the DOF is face-based, otherwise 0
  ! Iwhere(4) = currently processed element number.
  ! If Iwhere(1) = Iwhere(2) = Iwhere(3) = 0, then the DOF is element based.
  integer, dimension(4), intent(in)              :: Iwhere
  
  ! The coordinates of the point which is currently processed, given in
  ! reference coordinates of the currently processed cell type (edge,face,element).
  ! If the DOF is vertice-based, then Dwhere is undefined.
  ! If the DOF is edge-based or element-based in 1D, then Dwhere has dimension 1.
  ! If the DOF is face-based or element-based in 2D, then Dwhere has dimension 2.
  ! IF the DOF is element-based in 3D, then Dwhere has dimension 3.
  real(DP), dimension(:), intent(in)             :: Dwhere

  ! The coordinates of the point for which the boundary values are to be
  ! calculated.
  real(DP), dimension(:), intent(in)             :: Dcoords

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional    :: rcollection

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
  real(DP), dimension(:), intent(out)            :: Dvalues
!</output>
  
!</subroutine>

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValuesFBC_3D(Icomponents,rdiscretisation,&
                                     Revaluation, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretefbc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions on fictitious boundary components. It calculates a special quantity 
  ! on the boundary, which is then used by the discretisation routines to 
  ! generate a discrete 'snapshot' of the (actually analytic) boundary conditions.
  !
  ! The routine must calculate the values on all elements of the element
  ! list Ielements simultaneously. Iwhere is a list with vertex or edge numbers
  ! where information is to be retrieved. Dvalues is filled with function values
  ! while Binside is set to TRUE for every vertex/edge that is inside of the
  ! corresponding fictitious boundary region (identified by rbcRegion).
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary: 
  !   Icomponents(1..SIZE(Icomponents)) defines the number of the solution component,
  !   the value should be calculated for 
  !   (e.g. 1=1st solution component, e.g. X-velocity, 
  !         2=2nd solution component, e.g. Y-velocity,...,
  !         3=3rd solution component, e.g. pressure)
  !   Example: Icomponents(:) = [1,2] -> Compute velues for X- and Y-velocity
  !     (1=x, 2=y component)
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_blockDiscretisation), intent(in)                     :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<inputoutput>
  ! A t_discreteFBCevaluation structure array that defines what to evaluate, 
  ! where to evaluate and which accepts the return values.
  ! This callback routine must check out the cinfoNeeded-entry in this structure
  ! to find out what to evaluate.
  ! The other entries in this structure describe where to evaluate.
  ! The result of the evaluation must be written into the p_Dvalues array entry
  ! in this structure.
  !
  ! The number of structures in this array depend on what to evaluate:
  !
  ! For Dirichlet boundary:
  !   revaluation contains as many entries as Icomponents; every entry in
  !   Icomponent corresponds to one entry in revaluation
  !   (so Icomponent(1)=1 defines to evaluate the X-velocity while the 
  !    values for the X-velocity are written to revaluation(1)\%p_Dvalues;
  !    Icomponent(2)=2 defines to evaluate the Y-velocity while the values 
  !    for the Y-velocity are written to revaluation(2)\%p_Dvalues, etc).
  !
  type(t_discreteFBCevaluation), dimension(:), intent(inout) :: Revaluation
!</inputoutput>
  
!</subroutine>

      ! local variables
      real(DP) :: ddistance, dxcenter, dycenter,dzcenter, dradius, dx, dy,eps,dz
      real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
      type(t_triangulation), pointer :: p_rtriangulation
      integer :: ipoint,idx
      
      eps=0.016_dp
      
      ! Just make sure we are evaluating in the corners.
      if (Revaluation(1)%cinfoNeeded .ne. DISCFBC_NEEDFUNC) then
        print *,'FBC: only corner evaluation supported at the moment!'
        stop
      end if
      
      ! Get the triangulation array for the point coordinates
      p_rtriangulation => rdiscretisation%RspatialDiscr(1)%p_rtriangulation
      call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                     p_DvertexCoordinates)

      ! Definition of the circle
      dxcenter = 0.5
      dycenter = 0.5
      dzcenter = 0.5
      dradius  = 0.2
      
      ! Loop through the points where to evaluate:
      do idx = 1,Revaluation(1)%nvalues
      
        ! Get the number of the point to process
        ipoint = Revaluation(1)%p_Iwhere(idx)
        
        ! Get x- and y-coordinate
        dx = p_DvertexCoordinates(1,ipoint)
        dy = p_DvertexCoordinates(2,ipoint)
        dz = p_DvertexCoordinates(3,ipoint)        
        
        ! Get the distance to the center
        ddistance = sqrt( (dx-dxcenter)**2 + (dy-dycenter)**2 + (dz-dzcenter)**2 )
        ddistance=ddistance-dradius
        ! Point inside?
        if(abs(ddistance) .le. eps)then
        !if(ddistance .le. dradius)then
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here.
          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
        
        end if
        
      end do

    
  end subroutine

end module

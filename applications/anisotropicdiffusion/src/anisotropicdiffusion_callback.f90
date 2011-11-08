!##############################################################################
!# ****************************************************************************
!# <name> anisotropicdiffusion_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) coeff_Laplace
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getBoundaryValues
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 4.) getBoundaryValuesFBC
!#
!#     -> Returns analytic values in the inner of the domain on
!#        fictitious boundary objects
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcfassembly.inc'
!#
!# 5.) getReferenceFunction
!#
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 6.) getMonitorFunction
!#
!#     -> Returns the values of the monitor function which is used to
!#        perform h-adaptivity
!# </purpose>
!##############################################################################

module anisotropicdiffusion_callback

  use fsystem
  use storage
  use genoutput
  use mprimitives
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

  subroutine coeff_Laplace (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer :: isolution
    real(DP), dimension(2,2) :: A

    ! Get information about the solution from the collection --
    ! as configured in the main program.
    isolution = rcollection%IquickAccess(1)
    A(1,1) = rcollection%DquickAccess(1)
    A(1,2) = rcollection%DquickAccess(2)
    A(2,1) = rcollection%DquickAccess(3)
    A(2,2) = rcollection%DquickAccess(4)

    select case (isolution)
    case (0,2)

      !    u(x,y) = 16*x*(1-x)*y*(1-y)
      ! => f(x,y) = 32 * (y*(1-y)+x*(1-x))
      !    + Rotation, given by the matrix A
      
      Dcoefficients (1,:,:) = -(&
        -32.0_DP*A(1,1)*DPoints(2,:,:)*(1 - DPoints(2,:,:)) + &
          A(1,2)*(16.0_DP*(1.0_DP - DPoints(1,:,:))*(1.0_DP - DPoints(2,:,:)) &
        -16.0_DP*(1.0_DP - DPoints(1,:,:))*DPoints(2,:,:) &
        -16.0_DP*DPoints(1,:,:)*(1.0_DP - DPoints(2,:,:)) &
        +16.0_DP*DPoints(1,:,:)*DPoints(2,:,:)) &
        +A(2,1)*(16.0_DP*(1.0_DP - DPoints(1,:,:))*(1 - DPoints(2,:,:)) &
        -16.0_DP*(1.0_DP - DPoints(1,:,:))*DPoints(2,:,:) &
        -16.0_DP*DPoints(1,:,:)*(1 - DPoints(2,:,:)) &
        +16.0_DP*DPoints(1,:,:)*DPoints(2,:,:)) &
        -32.0_DP*A(2,2)*DPoints(1,:,:)*(1 - DPoints(1,:,:)) )

    case (3)
      Dcoefficients (1,:,:) = 32.0_DP * &
                      ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
                        Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )
                        
    case DEFAULT
    
      Dcoefficients (1,:,:) = 0.0_DP
                        
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction (cderivative,rdiscretisation, &
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
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)                      :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    integer :: isolution
    real(DP), dimension(2,2) :: A

    ! Get information about the solution from the collection --
    ! as configured in the main program.
    isolution = rcollection%IquickAccess(1)
    A(1,1) = rcollection%DquickAccess(1)
    A(1,2) = rcollection%DquickAccess(2)
    A(2,1) = rcollection%DquickAccess(3)
    A(2,2) = rcollection%DquickAccess(4)

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = 16*x*(1-x)*y*(1-y)
    Dvalues (:,:) = 16.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) * &
                              Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
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
  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
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
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  integer, intent(in)                                          :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   dwhere = 0 (not used)
  real(DP), intent(in)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    real(DP) :: dx,dy
    integer :: isolution
    real(DP), dimension(2,2) :: A

    ! Get information about the solution from the collection --
    ! as configured in the main program.
    isolution = rcollection%IquickAccess(1)
    A(1,1) = rcollection%DquickAccess(1)
    A(1,2) = rcollection%DquickAccess(2)
    A(2,1) = rcollection%DquickAccess(3)
    A(2,2) = rcollection%DquickAccess(4)

    ! Get the X/Y-coordinates of the boundary point:
    !
    call boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

!    ! Return zero Dirichlet boundary values for all situations.
!    SELECT CASE (isolution)
!    CASE (0)
!      Dvalues(1) = 0.0_DP
!
!    CASE (1)
!      Dvalues(1) = 0.0_DP
!      IF ((dwhere .GE. 0.0_DP) .AND. (dwhere .LE. 0.2_DP)) Dvalues(1) = 1.0_DP
!      IF ((dwhere .GE. 0.2_DP) .AND. (dwhere .LE. 0.3_DP)) &
!        CALL mprim_linearRescale(dwhere,0.2_DP,0.3_DP,1.0_DP,0.5_DP,Dvalues(1))
!      IF ((dwhere .GE. 0.3_DP) .AND. (dwhere .LE. 1.7_DP)) Dvalues(1) = 0.5_DP
!      IF ((dwhere .GE. 1.7_DP) .AND. (dwhere .LE. 1.8_DP)) &
!        CALL mprim_linearRescale(dwhere,1.7_DP,1.8_DP,0.5_DP,0.0_DP,Dvalues(1))
!      IF ((dwhere .GE. 1.8_DP) .AND. (dwhere .LE. 2.2_DP)) Dvalues(1) = 0.0_DP
!      IF ((dwhere .GE. 2.2_DP) .AND. (dwhere .LE. 2.3_DP)) &
!        CALL mprim_linearRescale(dwhere,2.2_DP,2.3_DP,0.0_DP,0.5_DP,Dvalues(1))
!      IF ((dwhere .GE. 2.3_DP) .AND. (dwhere .LE. 3.7_DP)) Dvalues(1) = 0.5_DP
!      IF ((dwhere .GE. 3.7_DP) .AND. (dwhere .LE. 3.8_DP)) &
!        CALL mprim_linearRescale(dwhere,3.7_DP,3.8_DP,0.5_DP,1.0_DP,Dvalues(1))
!      IF ((dwhere .GE. 3.8_DP) .AND. (dwhere .LE. 4.0_DP)) Dvalues(1) = 1.0_DP
!
!    CASE DEFAULT
!      Dvalues(1) = 0.0_DP
!
!    END SELECT
    Dvalues(1) = analyticalFunction (isolution,dx,dy)
  
  end subroutine

  ! ***************************************************************************

!<function>

  real(DP) function analyticalFunction (ifunction,dx,dy)
  
!<description>
  ! Calculates the value of an analytical function in a point (dx,dy).
  ! ifunction is an identifier for the function to be evaluated.
!</description>
  
!<input>
  ! Identifier for the function to evaluate.
  integer, intent(in) :: ifunction

  ! Point where to evaluate
  real(DP), intent(in) :: dx,dy
!</input>

!<return>
  ! Value of the function.
!</return>
  
!</subroutine>

    real(DP) :: dvalue

    ! Return zero Dirichlet boundary values for all situations.
    select case (ifunction)
    case (0)
      dvalue = 0.0_DP
      
    case (1)
      dvalue = 0.0_DP
      if (dy .eq. 0.0_DP) then
        if ((dx .ge. 0.0_DP) .and. (dx .le. 0.2_DP)) dvalue = 1.0_DP
        if ((dx .ge. 0.2_DP) .and. (dx .le. 0.3_DP)) &
          call mprim_linearRescale(dx,0.2_DP,0.3_DP,1.0_DP,0.5_DP,dvalue)
        if ((dx .ge. 0.3_DP) .and. (dx .le. 1.0_DP)) dvalue = 0.5_DP
      else if (dx .eq. 1.0_DP) then
        if ((dy .ge. 0.0_DP) .and. (dy .le. 0.7_DP)) dvalue = 0.5_DP
        if ((dy .ge. 0.7_DP) .and. (dy .le. 0.8_DP)) &
          call mprim_linearRescale(dy,0.7_DP,0.8_DP,0.5_DP,0.0_DP,dvalue)
        if ((dy .ge. 0.8_DP) .and. (dy .le. 1.0_DP)) dvalue = 0.0_DP
      else if (dy .eq. 1.0_DP) then
        if ((dx .ge. 0.8_DP) .and. (dx .le. 1.0_DP)) dvalue = 0.0_DP
        if ((dx .ge. 0.7_DP) .and. (dx .le. 0.8_DP)) &
          call mprim_linearRescale(dx,0.7_DP,0.8_DP,0.5_DP,0.0_DP,dvalue)
        if ((dx .ge. 0.0_DP) .and. (dx .le. 0.7_DP)) dvalue = 0.5_DP
      else if (dx .eq. 0.0_DP) then
        if ((dy .ge. 0.3_DP) .and. (dy .le. 1.0_DP)) dvalue = 0.5_DP
        if ((dy .ge. 0.2_DP) .and. (dy .le. 0.3_DP)) &
          call mprim_linearRescale(dy,0.2_DP,0.3_DP,1.0_DP,0.5_DP,dvalue)
        if ((dy .ge. 0.0_DP) .and. (dy .le. 0.2_DP)) dvalue = 1.0_DP
      end if
    case (2)
      dvalue = 16.0_DP*dx*(1-dx)*dy*(1-dy)
    case (4)
      if ((dx .eq. 0.0_DP) .or. (dy .eq. 0.0_DP) .or. &
          (dx .eq. 1.0_DP) .or. (dy .eq. 1.0_DP)) then
        dvalue = 0.0_DP
      else
        dvalue = 2.0_DP
      end if
      
    case DEFAULT
      dvalue = 0.0_DP
      
    end select
    
    analyticalFunction = dvalue
  
  end function

  ! ***************************************************************************

!<subroutine>

  subroutine getMonitorFunction(rtriangulation,rsolution,ieltype,ierrorestimator,rindicator)
  
    use pprocgradients
    use pprocerror

!<description>
  ! This routine defines a 'monitor function' for the adaptive grid refinement
  ! with the h-adaptivity refinement strategy. rindicator is a vector with
  ! NEL entries for all the elements in the triangulation. The routine must
  ! fill each entry with a value that tells the h-adaptivity routines whether
  ! to refine that element or not.
!</descrition>
  
!<input>
    ! The triangulation structure of the underlying mesh which is to be refined
    type(t_triangulation), intent(in) :: rtriangulation

    ! The solution vector
    type(t_vectorScalar), intent(inout)  :: rsolution

    ! The type of element used for the FE solution
    integer, intent(in) :: ieltype

    ! The type of error estimator
    integer, intent(in) :: ierrorestimator
!</input>
    
!</inputoutput>
    ! An indicator vector. Entry i in the vector rindicatir that tells the
    ! mesh adaption routines whether to refine element i or to do coarsening
    ! with it. A value > 1.0 will refine element i, a value < 0.01 will result
    ! in coarsening -- as specified during the initialisation of the
    ! mesh refinement in the main program.
    type(t_vectorScalar), intent(inout) :: rindicator

!</subroutine>

    ! local variables
    type(t_vectorBlock)         :: rgradient,rgradientRef
    type(t_blockDiscretisation) :: rdiscrBlock,rdiscrBlockRef
    real(DP)                    :: dsolutionError,dgradientError,daux

    ! Initialise block discretisations
    call spdiscr_initBlockDiscr (rdiscrBlock,2,&
        rtriangulation, rsolution%p_rspatialdiscr%p_rboundary)
    call spdiscr_initBlockDiscr (rdiscrBlockRef,2,&
        rtriangulation, rsolution%p_rspatialdiscr%p_rboundary)

    ! What kind of element type is used for the FE solution
    select case(ieltype)
    case(1)
      ! Initialise spatial discretisations for gradient with P0-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E000, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E000, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))
      
    case(2)
      ! Initialise spatial discretisations for gradient with P1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P2-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E002, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E002, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))

    case(11)
      ! Initialise spatial discretisations for gradient with Q0-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E010, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E010, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with Q1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))

    case(13)
      ! Initialise spatial discretisations for gradient with Q1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with Q2-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E013, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialdiscr,&
          EL_E013, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))

    case(-1)
      ! Initialise spatial discretisations for gradient with P0/Q0-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P1/Q1-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlockRef%RspatialDiscr(2))
      
    case(-2)
      ! Initialise spatial discretisations for gradient with P1/Q1-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P2/Q2-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E002, EL_E013, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialdiscr,&
          EL_E002, EL_E013, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlockRef%RspatialDiscr(2))

    case DEFAULT
      call output_line('Unsupproted element type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'getMonitorFunction')
      call sys_halt()
    end select
    
    ! Create block vector for gradient values
    call lsysbl_createVecBlockByDiscr (rdiscrBlock,   rgradient,    .true.)
    call lsysbl_createVecBlockByDiscr (rdiscrBlockRef,rgradientRef, .true.)

    ! Recover consistent gradient
    call ppgrd_calcGradient (rsolution, rgradient)

    ! Recover smoothed gradient
    select case(ierrorestimator)
    case (1)
      call ppgrd_calcGradient (rsolution, rgradientRef)

    case (2)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_NODEPATCH)

    case (3)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_ELEMPATCH)

    case (4)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_FACEPATCH)

    case DEFAULT
      call output_line('Invalid type of error estimator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'getMonitorFunction')
      call sys_halt()
    end select

    ! Compute gradient error
    call pperr_blockErrorEstimate(rgradient,rgradientRef,&
        PPERR_L2ERROR, dgradientError,relementError=rindicator)

    ! rindicator is currently the absolute error. If the relative error
    ! is to be used, the following lines must be commented in:
    !
    ! ! Compute L2-norm of solution
    ! CALL pperr_scalar(rsolution,PPERR_L2ERROR,dsolutionError)
    !
    ! ! Prepare indicator for grid refinement/coarsening
    ! daux=SQRT((dsolutionError**2+dgradientError**2)/REAL(rindicator%NEQ,DP))
    ! CALL lsyssc_scaleVector(rindicator,1._DP/daux)

    call output_line('!!gradient error!! = '//trim(sys_sdEL(dgradientError,10)))
    
    ! Release temporal discretisation structure
    call spdiscr_releaseBlockDiscr(rdiscrBlock)
    call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
    
    ! Release vectors
    call lsysbl_releaseVector(rgradient)
    call lsysbl_releaseVector(rgradientRef)
  end subroutine

end module

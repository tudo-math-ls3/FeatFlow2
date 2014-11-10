!##############################################################################
!# ****************************************************************************
!# <name> elemdbg3d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module elemdbg3d_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use spatialdiscretisation
  use discretebc
  use bcassembly
  use derivatives
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS3D (rdiscretisation,rform, &
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
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    integer :: itest, isolution
    real(DP) :: dnu, dbeta1, dbeta2, dbeta3
    
    ! Get the info from the collection
    itest = rcollection%IquickAccess(1)
    isolution = rcollection%IquickAccess(2)
    dnu = rcollection%DquickAccess(1)
    dbeta1 = rcollection%DquickAccess(2)
    dbeta2 = rcollection%DquickAccess(3)
    dbeta3 = rcollection%DquickAccess(4)

    ! Okay, what type of system do we have here?
    select case(itest)
    case(301)
      ! L2-projection
      select case(isolution)
      case(0)
        Dcoefficients(1,:,:) = sin(SYS_PI*Dpoints(1,:,:)) &
            * sin(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:))
      case(1)
        Dcoefficients(1,:,:) = Dpoints(1,:,:)
      case(2)
        Dcoefficients(1,:,:) = Dpoints(2,:,:)
      case(3)
        Dcoefficients(1,:,:) = Dpoints(3,:,:)
      end select
                   
    case(302)
      ! Poisson-System
      select case(isolution)
      case(0)
        Dcoefficients(1,:,:) = 3.0_DP * SYS_PI**2 * sin(SYS_PI*Dpoints(1,:,:)) &
            * sin(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:))
      case(1,2,3)
        Dcoefficients(1,:,:) = 0.0_DP
      end select
    
    case(303)
      ! Convection-Diffusion-System
      select case(isolution)
      case(0)
        Dcoefficients(1,:,:) = 3.0_DP * dnu * SYS_PI**2 * sin(SYS_PI*Dpoints(1,:,:)) &
            * sin(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:)) &
          + dbeta1 * SYS_PI * cos(SYS_PI*Dpoints(1,:,:)) &
            * sin(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:)) &
          + dbeta2 * SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) &
            * cos(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:)) &
          + dbeta3 * SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) &
            * sin(SYS_PI*Dpoints(2,:,:)) * cos(SYS_PI*Dpoints(3,:,:))
      case(1)
        Dcoefficients(1,:,:) = dbeta1
      case(2)
        Dcoefficients(1,:,:) = dbeta2
      case(3)
        Dcoefficients(1,:,:) = dbeta3
      end select
          
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction3D (icomponent,cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                rdomainIntSubset,Dvalues,rcollection)
  
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
  ! Specifies which component of the vector field is to be evaluated
  integer, intent(IN) :: icomponent

  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN) :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>
  
!</subroutine>
    
    Dvalues = 0.0_DP

    select case(rcollection%IquickAccess(2))
    case(0)
      select case (cderivative)
      case (DER_FUNC3D)
        Dvalues(:,:) = sin(SYS_PI*Dpoints(1,:,:)) &
                     * sin(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:))
      case (DER_DERIV3D_X)
        Dvalues(:,:) = SYS_PI * cos(SYS_PI*Dpoints(1,:,:)) &
                     * sin(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:))
      case (DER_DERIV3D_Y)
        Dvalues(:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) &
                     * cos(SYS_PI*Dpoints(2,:,:)) * sin(SYS_PI*Dpoints(3,:,:))
      case (DER_DERIV3D_Z)
        Dvalues(:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) &
                     * sin(SYS_PI*Dpoints(2,:,:)) * cos(SYS_PI*Dpoints(3,:,:))
      end select
    
    case(1)
      select case (cderivative)
      case (DER_FUNC3D)
        Dvalues(:,:) = Dpoints(1,:,:)
      case (DER_DERIV3D_X)
        Dvalues(:,:) = 1.0_DP
      end select
    
    case(2)
      select case (cderivative)
      case (DER_FUNC3D)
        Dvalues(:,:) = Dpoints(2,:,:)
      case (DER_DERIV3D_Y)
        Dvalues(:,:) = 1.0_DP
      end select
    
    case(3)
      select case (cderivative)
      case (DER_FUNC3D)
        Dvalues(:,:) = Dpoints(3,:,:)
      case (DER_DERIV3D_Z)
        Dvalues(:,:) = 1.0_DP
      end select
    
    end select
  
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine getBoundaryValues3D (Icomponents,rdiscretisation,rmeshRegion,&
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(IN)                              :: rmeshRegion

  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(IN)                                         :: cinfoNeeded
  
  ! An array holding information about what type of DOF is currently processed.
  ! The information is build up as follows:
  ! Iwhere(1) = vertice number of the DOF, if the DOF is vertice-based, otherwise 0
  ! Iwhere(2) = edge number of the DOF, if the DOF is edge-based, otherwise 0
  ! Iwhere(3) = face number of the DOF, if the DOF is face-based, otherwise 0
  ! Iwhere(4) = currently processed element number.
  ! If Iwhere(1) = Iwhere(2) = Iwhere(3) = 0, then the DOF is element based.
  integer, dimension(4), intent(IN)                           :: Iwhere
  
  ! The coordinates of the point which is currently processed, given in
  ! reference coordinates of the currently processed cell type (edge,face,element).
  ! If the DOF is vertice-based, then Dwhere is undefined.
  ! If the DOF is edge-based or element-based in 1D, then Dwhere has dimension 1.
  ! If the DOF is face-based or element-based in 2D, then Dwhere has dimension 2.
  ! IF the DOF is element-based in 3D, then Dwhere has dimension 3.
  real(DP), dimension(:), intent(IN)                          :: Dwhere

  ! The coordinates of the point for which the boundary values are to be
  ! calculated.
  real(DP), dimension(:), intent(IN)                          :: Dcoords

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional                 :: rcollection

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
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

    select case(rcollection%IquickAccess(2))
    case(0)
      Dvalues(1) = 0.0_DP
    case(1)
      Dvalues(1) = Dcoords(1)
    case(2)
      Dvalues(1) = Dcoords(2)
    case(3)
      Dvalues(1) = Dcoords(3)
    end select

  end subroutine

end module

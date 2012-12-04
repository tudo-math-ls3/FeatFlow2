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
!# 1.) coeff_RHS_2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 2D case, Q2 bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 2.) coeff_RHS_Sin2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 2D case, sinus bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getBoundaryValues_2D
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 4.) getBoundaryValuesFBC_2D
!#     -> Returns analytic values in the inner of the domain on
!#        fictitious boundary objects
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcfassembly.inc'
!#
!# 5.) getBoundaryValuesMR_2D
!#     -> Returns discrete values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_discretebc.inc'
!#
!# 6.) getReferenceFunction_2D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_2D, Q2 bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 7.) getReferenceFunction_Sin2D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_Sin2D, sinus bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 8.) gethadaptMonitorFunction_2D
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

  subroutine coeff_RHS_2D (rdiscretisation,rform, &
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

    !    u(x,y) = 16*x*(1-x)*y*(1-y)
    ! => f(x,y) = 32 * (y*(1-y)+x*(1-x))
    Dcoefficients (1,:,:) = 0.0_DP
    !32.0_DP * &
    !                ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
    !                  Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_2D (cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = 16*x*(1-x)*y*(1-y)
    Dvalues (:,:) = Dpoints(1,:,:) * Dpoints(2,:,:)
  case (DER_DERIV_X)
    !    u(x,y)   = 16*x*(1-x)*y*(1-y)
    ! => u_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
    Dvalues (:,:) = Dpoints(2,:,:)
  case (DER_DERIV_Y)
    !    u(x,y)   = 16*x*(1-x)*y*(1-y)
    ! => u_y(x,y) = 16 * ( x*(1-x)*(1-y)-x*y*(1-x) )
    Dvalues (:,:) = Dpoints(1,:,:)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Sin2D (rdiscretisation,rform, &
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

    !    u(x,y) = SIN(PI * x) * SIN(PI * y)
    ! => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)
    Dcoefficients (1,:,:) = 2.0_DP * SYS_PI**2 &
                          * sin(SYS_PI * Dpoints(1,:,:)) &
                          * sin(SYS_PI * Dpoints(2,:,:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_Sin2D (cderivative,rdiscretisation, &
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

  select case (cderivative)
  case (DER_FUNC)
    ! u(x,y) = SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = SYS_PI * cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    REAL(DP) :: dx,dy
    
    CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP !dx*dy
    if (dx .eq. 0.0_dp) then
      dvalues(1) = -0.5_dp
    end if
    if (dx .eq. 1.0_dp) then
      dvalues(1) = 0.5_dp
    end if
    if (dy .eq. 0.0_dp) then
      dvalues(1) = dx-0.5_dp
    end if
    if (dy .eq. 1.0_DP) then
      Dvalues(1) = 0.0_DP
    end if
  
  end subroutine

  ! ***************************************************************************
  ! Only for poisson2d_method1_fbc: Values in a fictitious boundary component:

!<subroutine>

  subroutine getBoundaryValuesFBC_2D (Icomponents,rdiscretisation,&
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
      real(DP) :: ddistance, dxcenter, dycenter, dradius, dx, dy
      real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
      type(t_triangulation), pointer :: p_rtriangulation
      integer :: ipoint,idx

      
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
      dxcenter = 0.4
      dycenter = 0.4
      dradius  = 0.25
      
      ! Loop through the points where to evaluate:
      do idx = 1,Revaluation(1)%nvalues
      
        ! Get the number of the point to process
        ipoint = Revaluation(1)%p_Iwhere(idx)
        
        ! Get x- and y-coordinate
        dx = p_DvertexCoordinates(1,ipoint)
        dy = p_DvertexCoordinates(2,ipoint)
        
        ! Get the distance to the center
        ddistance = sqrt( (dx-dxcenter)**2 + (dy-dycenter)**2 )
        
        ! Point inside?
        if (ddistance .le. dradius) then
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here.
          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
        
        end if
        
      end do

      ! Definition of a 2nd circle
      dxcenter = 0.75
      dycenter = 0.75
      dradius  = 0.1
      
      ! Loop through the points where to evaluate:
      do idx = 1,Revaluation(1)%nvalues
      
        ! Get the number of the point to process
        ipoint = Revaluation(1)%p_Iwhere(idx)
        
        ! Get x- and y-coordinate
        dx = p_DvertexCoordinates(1,ipoint)
        dy = p_DvertexCoordinates(2,ipoint)
        
        ! Get the distance to the center
        ddistance = sqrt( (dx-dxcenter)**2 + (dy-dycenter)**2 )
        
        ! Point inside?
        if (ddistance .le. dradius) then
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here.
          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
        
        end if
        
      end do
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine getBoundaryValuesMR_2D (Icomponents,rdiscretisation,rmeshRegion,&
                                      cinfoNeeded,Dcoords,Dvalues,rcollection)
  
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
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in)                              :: rmeshRegion

  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
  ! The coordinates of the point for which the boundary values are to be
  ! calculated.
  real(DP), dimension(:), intent(in)                          :: Dcoords

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional                 :: rcollection

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

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP

  end subroutine

  ! ***************************************************************************
  ! Only for poisson2d_method1_hadapt: Monitor function for adaptive grid refinement
  
!<subroutine>

  subroutine gethadaptMonitorFunction_2D(rtriangulation,rsolution,ieltype,&
      ierrorestimator,rindicator)
  
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
    integer(I32), intent(in) :: ieltype
    
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
    type(t_scalarCubatureInfo) :: rcubatureInfo

    ! Initialise block discretisations
    call spdiscr_initBlockDiscr (rdiscrBlock,2,&
        rtriangulation, rsolution%p_rspatialDiscr%p_rboundary)
    call spdiscr_initBlockDiscr (rdiscrBlockRef,2,&
        rtriangulation, rsolution%p_rspatialDiscr%p_rboundary)

    ! What kind of element type is used for the FE solution
    select case(ieltype)
    case(1)
      ! Initialise spatial discretisations for gradient with P0-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E000, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E000, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_P1, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_P1, rdiscrBlockRef%RspatialDiscr(2))
      
    case(2)
      ! Initialise spatial discretisations for gradient with P1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_P1, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_P1, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P2-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E002, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E002, rdiscrBlockRef%RspatialDiscr(2))

    case(11)
      ! Initialise spatial discretisations for gradient with Q0-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E010, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E010, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with Q1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_Q1, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_Q1, rdiscrBlockRef%RspatialDiscr(2))

    case(13)
      ! Initialise spatial discretisations for gradient with Q1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_Q1, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_Q1, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with Q2-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E013, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E013, rdiscrBlockRef%RspatialDiscr(2))

    case(-1)
      ! Initialise spatial discretisations for gradient with P0/Q0-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_E000, EL_E010, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_E000, EL_E010, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P1/Q1-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_P1, EL_Q1, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_P1, EL_Q1, rdiscrBlockRef%RspatialDiscr(2))
      
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

    call spdiscr_createDefCubStructure(&  
        rsolution%p_rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G2)

    ! Recover smoothed gradient
    select case(ierrorestimator)
    case (1)
      call ppgrd_calcGradient (rsolution, rgradientRef)

    case (2)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_NODEPATCH,rcubatureInfo)

    case (3)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_ELEMPATCH,rcubatureInfo)

    case (4)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_FACEPATCH,rcubatureInfo)

    case DEFAULT
      call output_line("Invalid type of error estimator!",&
          OU_CLASS_ERROR,OU_MODE_STD,"getMonitorFunction")
      call sys_halt()
    end select

    ! Compute gradient error
    call pperr_blockErrorEstimate(rgradient,rgradientRef,PPERR_L2ERROR,&
        dgradientError,rcubatureInfo,relementError=rindicator)
    call output_line("!!gradient error!! = "//&
        trim(sys_sdEL(dgradientError,10)))

    ! Compute L2-norm of solution
    call pperr_scalar(PPERR_L2ERROR,dsolutionError,rsolution)

    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Prepare indicator for grid refinement/coarsening
    daux=sqrt((dsolutionError**2+dgradientError**2)/real(rindicator%NEQ,DP))
    call lsyssc_scaleVector(rindicator,1._DP/daux)
    
    ! Release temporal discretisation structure
    call spdiscr_releaseBlockDiscr(rdiscrBlock)
    call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
    
    ! Release vectors
    call lsysbl_releaseVector(rgradient)
    call lsysbl_releaseVector(rgradientRef)
    
  end subroutine

end module

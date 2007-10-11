!##############################################################################
!# ****************************************************************************
!# <name> poissoncallback </name>
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
!# </purpose>
!##############################################################################

MODULE anisotropicdiffusion_callback

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE mprimitives
  
  IMPLICIT NONE

CONTAINS

! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_Laplace (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, p_rcollection,&
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    ! analytic boundary boundary description etc.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The bilinear form which is currently being evaluated:
    TYPE(t_bilinearForm), INTENT(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 1.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_RHS (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    INTEGER :: isolution
    REAL(DP), DIMENSION(2,2) :: A

    ! Get information about the solution from the collection --
    ! as configured in the main program.
    isolution = p_rcollection%IquickAccess(1) 
    A(1,1) = p_rcollection%DquickAccess(1) 
    A(1,2) = p_rcollection%DquickAccess(2) 
    A(2,1) = p_rcollection%DquickAccess(3) 
    A(2,2) = p_rcollection%DquickAccess(4) 

    SELECT CASE (isolution)
    CASE (0,2)

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

    CASE (3)             
      Dcoefficients (1,:,:) = 32.0_DP * &
                      ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
                        Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )
                        
    CASE DEFAULT
    
      Dcoefficients (1,:,:) = 0.0_DP
                        
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE getReferenceFunction (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,p_rcollection, &
                Dvalues)
  
  USE basicgeometry
  USE triangulation
  USE collection
  USE scalarpde
  USE domainintegration
  
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
  INTEGER, INTENT(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  INTEGER, INTENT(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  INTEGER, INTENT(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                      :: p_rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  REAL(DP), DIMENSION(:,:), INTENT(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

    INTEGER :: isolution
    REAL(DP), DIMENSION(2,2) :: A

    ! Get information about the solution from the collection --
    ! as configured in the main program.
    isolution = p_rcollection%IquickAccess(1) 
    A(1,1) = p_rcollection%DquickAccess(1) 
    A(1,2) = p_rcollection%DquickAccess(2) 
    A(2,1) = p_rcollection%DquickAccess(3) 
    A(2,2) = p_rcollection%DquickAccess(4) 

  SELECT CASE (cderivative)
  CASE (DER_FUNC)
    ! u(x,y) = 16*x*(1-x)*y*(1-y)
    Dvalues (:,:) = 16.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) * &
                              Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
  CASE (DER_DERIV_X)
    !    u(x,y)   = 16*x*(1-x)*y*(1-y)
    ! => u_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
    Dvalues (:,:) = 16.0_DP * ( &
        Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
        Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(2,:,:)) )
  CASE (DER_DERIV_Y)
    !    u(x,y)   = 16*x*(1-x)*y*(1-y)
    ! => u_y(x,y) = 16 * ( x*(1-x)*(1-y)-x*y*(1-x) )
    Dvalues (:,:) = 16.0_DP * ( &
        Dpoints(1,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
        Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) )
  CASE DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  END SELECT
  

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE getBoundaryValues (Icomponents,rdiscretisation,rbcRegion,ielement, &
                                cinfoNeeded,iwhere,dwhere, p_rcollection, Dvalues)
  
  USE collection
  USE spatialdiscretisation
  USE discretebc
  
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
  INTEGER, DIMENSION(:), INTENT(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Boundary condition region that is currently being processed.
  ! (This e.g. defines the type of boundary conditions that are
  !  currently being calculated, as well as information about the current
  !  boundary segment 'where we are at the moment'.)
  TYPE(t_bcRegion), INTENT(IN)                                :: rbcRegion
  
  
  ! The element number on the boundary which is currently being processed
  INTEGER(I32), INTENT(IN)                                    :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  INTEGER, INTENT(IN)                                         :: cinfoNeeded
  
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
  INTEGER(I32), INTENT(IN)                                     :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  REAL(DP), INTENT(IN)                                        :: dwhere
    
  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                  :: p_rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
!</output>
  
!</subroutine>
    REAL(DP) :: dx,dy
    INTEGER :: isolution
    REAL(DP), DIMENSION(2,2) :: A

    ! Get information about the solution from the collection --
    ! as configured in the main program.
    isolution = p_rcollection%IquickAccess(1) 
    A(1,1) = p_rcollection%DquickAccess(1) 
    A(1,2) = p_rcollection%DquickAccess(2) 
    A(2,1) = p_rcollection%DquickAccess(3) 
    A(2,2) = p_rcollection%DquickAccess(4) 

    ! Get the X/Y-coordinates of the boundary point:
    !
    CALL boundary_getCoords(rdiscretisation%p_rboundary, &
        rbcRegion%rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

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
  
  END SUBROUTINE

  ! ***************************************************************************

!<function>

  REAL(DP) FUNCTION analyticalFunction (ifunction,dx,dy)
  
!<description>
  ! Calculates the value of an analytical function in a point (dx,dy).
  ! ifunction is an identifier for the function to be evaluated.
!</description>
  
!<input>
  ! Identifier for the function to evaluate.
  INTEGER, INTENT(IN) :: ifunction

  ! Point where to evaluate
  REAL(DP), INTENT(IN) :: dx,dy
!</input>

!<return>
  ! Value of the function.
!</return>
  
!</subroutine>

    REAL(DP) :: dvalue

    ! Return zero Dirichlet boundary values for all situations.
    SELECT CASE (ifunction)
    CASE (0)
      dvalue = 0.0_DP
      
    CASE (1)
      dvalue = 0.0_DP
      IF (dy .EQ. 0.0_DP) THEN
        IF ((dx .GE. 0.0_DP) .AND. (dx .LE. 0.2_DP)) dvalue = 1.0_DP
        IF ((dx .GE. 0.2_DP) .AND. (dx .LE. 0.3_DP)) &
          CALL mprim_linearRescale(dx,0.2_DP,0.3_DP,1.0_DP,0.5_DP,dvalue)
        IF ((dx .GE. 0.3_DP) .AND. (dx .LE. 1.0_DP)) dvalue = 0.5_DP
      ELSE IF (dx .EQ. 1.0_DP) THEN
        IF ((dy .GE. 0.0_DP) .AND. (dy .LE. 0.7_DP)) dvalue = 0.5_DP
        IF ((dy .GE. 0.7_DP) .AND. (dy .LE. 0.8_DP)) &
          CALL mprim_linearRescale(dy,0.7_DP,0.8_DP,0.5_DP,0.0_DP,dvalue)
        IF ((dy .GE. 0.8_DP) .AND. (dy .LE. 1.0_DP)) dvalue = 0.0_DP
      ELSE IF (dy .EQ. 1.0_DP) THEN
        IF ((dx .GE. 0.8_DP) .AND. (dx .LE. 1.0_DP)) dvalue = 0.0_DP
        IF ((dx .GE. 0.7_DP) .AND. (dx .LE. 0.8_DP)) &
          CALL mprim_linearRescale(dx,0.7_DP,0.8_DP,0.5_DP,0.0_DP,dvalue)
        IF ((dx .GE. 0.0_DP) .AND. (dx .LE. 0.7_DP)) dvalue = 0.5_DP
      ELSE IF (dx .EQ. 0.0_DP) THEN
        IF ((dy .GE. 0.3_DP) .AND. (dy .LE. 1.0_DP)) dvalue = 0.5_DP
        IF ((dy .GE. 0.2_DP) .AND. (dy .LE. 0.3_DP)) &
          CALL mprim_linearRescale(dy,0.2_DP,0.3_DP,1.0_DP,0.5_DP,dvalue)
        IF ((dy .GE. 0.0_DP) .AND. (dy .LE. 0.2_DP)) dvalue = 1.0_DP
      END IF
    CASE (2)
      dvalue = 16.0_DP*dx*(1-dx)*dy*(1-dy)
    CASE (4)
      IF ((dx .EQ. 0.0_DP) .OR. (dy .EQ. 0.0_DP) .OR. & 
          (dx .EQ. 1.0_DP) .OR. (dy .EQ. 1.0_DP)) THEN
        dvalue = 0.0_DP
      ELSE
        dvalue = 2.0_DP
      END IF
      
    CASE DEFAULT
      dvalue = 0.0_DP
      
    END SELECT
    
    analyticalFunction = dvalue
  
  END FUNCTION

END MODULE

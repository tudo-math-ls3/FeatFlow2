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
!#        equation for the 1D poisson-problem (i.e. poisson0). 3D case.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getReferenceFunction_3D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_3D for the 3D poisson-problem.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# </purpose>
!##############################################################################

MODULE poisson3d_callback

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
  
  IMPLICIT NONE

CONTAINS

! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_Laplace_3D (rdiscretisation,rform, &
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
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints
    
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

  SUBROUTINE coeff_RHS_3D (rdiscretisation,rform, &
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
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

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

    !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
    ! => f(x,y,z) = 128 * ( y*(1-y)*z*(1-z)
    !             + x*(1-x)*z*(1-z) + x*(1-x)*y*(1-y))
    Dcoefficients(1,:,:) = 128.0_DP * &
        ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))*&
            Dpoints(3,:,:)*(1.0_DP-Dpoints(3,:,:)) + &
          Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*&
            Dpoints(3,:,:)*(1.0_DP-Dpoints(3,:,:)) + &
          Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))*&
            Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)))
          

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE getReferenceFunction_3D (cderivative,rdiscretisation, &
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
  ! DIMENSION(dimension,npointsPerElement,nelements)
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

  SELECT CASE (cderivative)
  CASE (DER_FUNC3D)
    !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
    Dvalues(:,:) = 64.0_DP * Dpoints(1,:,:) * (1.0_DP - Dpoints(1,:,:))&
                           * Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))&
                           * Dpoints(3,:,:) * (1.0_DP - Dpoints(3,:,:))
  CASE (DER_DERIV3D_X)
    !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
    ! => u_x(x,y,z) = 64*(1-2*x)*y*(1-y)*z*(1-z)
    Dvalues(:,:) = 64.0_DP * (1.0_DP - 2.0_DP * Dpoints(1,:,:))&
                           * Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))&
                           * Dpoints(3,:,:) * (1.0_DP - Dpoints(3,:,:))
  CASE (DER_DERIV3D_Y)
    !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
    ! => u_y(x,y,z) = 64*(1-2*y)*x*(1-x)*z*(1-z)
    Dvalues(:,:) = 64.0_DP * (1.0_DP - 2.0_DP * Dpoints(2,:,:))&
                           * Dpoints(1,:,:) * (1.0_DP - Dpoints(1,:,:))&
                           * Dpoints(3,:,:) * (1.0_DP - Dpoints(3,:,:))
  CASE (DER_DERIV3D_Z)
    !    u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
    ! => u_y(x,y,z) = 64*(1-2*z)*x*(1-x)*y*(1-y)
    Dvalues(:,:) = 64.0_DP * (1.0_DP - 2.0_DP * Dpoints(3,:,:))&
                           * Dpoints(1,:,:) * (1.0_DP - Dpoints(1,:,:))&
                           * Dpoints(2,:,:) * (1.0_DP - Dpoints(2,:,:))
  CASE DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  END SELECT
  
  END SUBROUTINE

END MODULE

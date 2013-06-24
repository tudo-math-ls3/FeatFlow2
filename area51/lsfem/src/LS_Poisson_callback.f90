!##############################################################################
!# ****************************************************************************
!# <name> LS_Poisson_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the stokes problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) coeff_Stokes_2D
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 1.) coeff_Pressure_2D
!#     -> Returns the coefficients for the pressure matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS_X_2D, coeff_RHS_Y_2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation -- for the X-velocity as well as for the Y-velocity.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getBoundaryValues_2D
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

module LS_Poisson_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use spatialdiscretisation
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use genoutput
  use triangulation
  use derivatives
  use linearsystemscalar
  use linearsystemblock
  use element
  use paramlist
  
  implicit none

contains

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
  
  
! select case (cderivative)
! case (DER_FUNC)
!   ! phi(x,y) = SIN(PI * x) * SIN(PI * y)
!   Dvalues (:,:) = sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
! case (DER_DERIV_X)
!   !    phi(x,y)   = SIN(PI * x) * SIN(PI * y)
!   ! => phi_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
!   Dvalues (:,:) = SYS_PI * cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
! case (DER_DERIV_Y)
!  !    phi(x,y)   = SIN(PI * x) * SIN(PI * y)
!   ! => phi_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
!   Dvalues (:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
! case DEFAULT
!   ! Unknown. Set the result to 0.0.
!   Dvalues = 0.0_DP
! end select  

! ======================================================================
  ! Local variable
    real(DP) :: n
    ! parameter
    type(t_parlist) :: rparams  

    ! Ok, let us start. 
    ! reading data from the *.dat file 
    ! We want to solve our Poisson problem on level...
    call parlst_init(rparams)
    call parlst_readfromfile (rparams, "./data/ls.dat") 
    
    call parlst_getvalue_double (rparams, 'GENERAL', 'n', n, 1.0_DP) 
 
 select case (cderivative)
 case (DER_FUNC)
   ! phi(x,y) = SIN(n* PI * x)
   Dvalues (:,:) = sin(n*SYS_PI*Dpoints(1,:,:))*sin(n*SYS_PI*Dpoints(2,:,:))
 case (DER_DERIV_X)
   ! => phi_x(x,y) = n * PI * COS(n * PI * x)
   Dvalues (:,:) = n * SYS_PI * cos(n * SYS_PI*Dpoints(1,:,:))*sin(n*SYS_PI*Dpoints(2,:,:))
 case (DER_DERIV_Y)
   ! => phi_y(x,y) = 0
   Dvalues (:,:) = n * SYS_PI * cos(n * SYS_PI*Dpoints(2,:,:))*sin(n*SYS_PI*Dpoints(1,:,:))
 case DEFAULT
   ! Unknown. Set the result to 0.0.
   Dvalues = 0.0_DP
 end select   
  
! =========================================================================

  
!   select case (cderivative)
!   case (DER_FUNC)
!     ! phi(x,y) = SIN(x * (1-x)) * SIN(y * (1-y))
!     Dvalues (:,:) = sin(Dpoints(1,:,:)*(1-Dpoints(1,:,:))) * sin(Dpoints(2,:,:)*(1-Dpoints(2,:,:)))
!   case (DER_DERIV_X)
!     !    phi(x,y)   = SIN(x * (1-x)) * SIN(y * (1-y))
!     ! => phi_x(x,y) = [(1-2x)*COS(x * (1-x))] * SIN(y * (1-y))
!     Dvalues (:,:) = ((1-2*Dpoints(1,:,:)) * (cos(Dpoints(1,:,:)*(1-Dpoints(1,:,:))))) * sin(Dpoints(2,:,:)*(1-Dpoints(2,:,:)))
!   case (DER_DERIV_Y)
!    !    phi(x,y)   = SIN(x * (1-x)) * SIN(y * (1-y))
!     ! => phi_y(x,y) = [(1-2y)*COS(y * (1-y))] * SIN(x * (1-x))
!     Dvalues (:,:) = ((1-2*Dpoints(2,:,:)) * (cos(Dpoints(2,:,:)*(1-Dpoints(2,:,:))))) * sin(Dpoints(1,:,:)*(1-Dpoints(1,:,:)))
!   case DEFAULT
!     ! Unknown. Set the result to 0.0.
!     Dvalues = 0.0_DP
!   end select   
  
!   select case (cderivative)
!   case (DER_FUNC)
!     ! phi(x,y) = 8 * (x*x+y*y)
!     Dvalues (:,:) = 8.0_DP * (Dpoints(1,:,:)*Dpoints(1,:,:) + &
!                               Dpoints(2,:,:)*Dpoints(2,:,:))
!   case (DER_DERIV_X)
!     !    phi(x,y) = 8 * (x*x+y*y)
!     ! => phi_x(x,y) = 16 * x 
!     Dvalues (:,:) = 16.0_DP * Dpoints(1,:,:)
!   case (DER_DERIV_Y)
!     !    phi(x,y) = 8 * (x*x+y*y)
!     ! => phi_y(x,y) = 16 * y
!     Dvalues (:,:) = 16.0_DP * Dpoints(2,:,:)
!   case DEFAULT
!     ! Unknown. Set the result to 0.0.
!     Dvalues = 0.0_DP
!   end select  
  
!   select case (cderivative)
!   case (DER_FUNC)
!     ! u(x,y) = 16*x*(1-x)*y*(1-y)
!     Dvalues (:,:) = 4.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:))
!   case (DER_DERIV_X)
!     !    u(x,y)   = 16*x*(1-x)*y*(1-y)
!     ! => u_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
!     Dvalues (:,:) = 4.0_DP-8*Dpoints(1,:,:)
!   case (DER_DERIV_Y)
!     !    u(x,y)   = 16*x*(1-x)*y*(1-y)
!     ! => u_y(x,y) = 16 * ( x*(1-x)*(1-y)-x*y*(1-x) )
!     Dvalues (:,:) = 0.0_DP
!   case DEFAULT
!     ! Unknown. Set the result to 0.0.
!     Dvalues = 0.0_DP
!   end select
  
!   select case (cderivative)
!   case (DER_FUNC)
!     ! phi(x,y) = 16*x*(1-x)*y*(1-y)
!     Dvalues (:,:) = 16.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) * &
!                               Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))
!   case (DER_DERIV_X)
!     ! => phi_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
!     Dvalues (:,:) = 16.0_DP * ( &
!         Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
!         Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(2,:,:)) )
!   case (DER_DERIV_Y)
!     ! => phi_y(x,y) = 16 * ( x*(1-x)*(1-y)-x*y*(1-x) )
!     Dvalues (:,:) = 16.0_DP * ( &
!         Dpoints(1,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
!         Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) )
!   case DEFAULT
!     ! Unknown. Set the result to 0.0.
!     Dvalues = 0.0_DP
!   end select

  end subroutine

! ***************************************************************************
  
!<subroutine>

  subroutine getReferenceFunction_2Du (cderivative,rdiscretisation, &
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
!==========================================================  
  ! Local variable
    real(DP) :: n
    ! parameter
    type(t_parlist) :: rparams  

    ! Ok, let us start. 
    ! reading data from the *.dat file 
    ! We want to solve our Poisson problem on level...
    call parlst_init(rparams)
    call parlst_readfromfile (rparams, "./data/ls.dat") 
    
    call parlst_getvalue_double (rparams, 'GENERAL', 'n', n, 1.0_DP) 
!==========================================================  
  select case (cderivative)
  case (DER_FUNC)
    
    
    Dvalues (:,:) = -n*SYS_PI * cos(n*SYS_PI*Dpoints(1,:,:))*sin(n*SYS_PI*Dpoints(2,:,:))
    
    !    phi(x,y)   = SIN(x * (1-x)) * SIN(y * (1-y))
    ! => phi_x(x,y) = [(1-2x)*COS(x * (1-x))] * SIN(y * (1-y))
!     Dvalues (:,:) = -((1-2*Dpoints(1,:,:)) * (cos(Dpoints(1,:,:)*(1-Dpoints(1,:,:))))) * sin(Dpoints(2,:,:)*(1-Dpoints(2,:,:)))
  
    !    phi(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => phi_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
!    Dvalues (:,:) = -SYS_PI * cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))


    !    phi(x,y) = 8 * (x*x+y*y)
    ! => phi_x(x,y) = 16 * x 
!     Dvalues (:,:) = -16.0_DP * Dpoints(1,:,:)
    
    !    phi(x,y)   = 16*x*(1-x)*y*(1-y)
    ! => phi_x(x,y) = 16 * ( y*(1-x)*(1-y)-x*y*(1-y) )
!     Dvalues (:,:) = -16.0_DP * ( &
!        Dpoints(2,:,:) * (1.0_DP-Dpoints(1,:,:)) * (1.0_DP-Dpoints(2,:,:)) - &
!        Dpoints(1,:,:) * Dpoints(2,:,:) * (1.0_DP-Dpoints(2,:,:)) )

!     Dvalues (:,:) = -4.0_DP-8*Dpoints(1,:,:)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select

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
    ! DIMENSION(\#local DOF`s in test space,nelements)
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
  
  ! Local variable
    real(DP) :: kk, n
    ! parameter
    type(t_parlist) :: rparams  

    ! Ok, let us start. 
    ! reading data from the *.dat file 
    ! We want to solve our Poisson problem on level...
    call parlst_init(rparams)
    call parlst_readfromfile (rparams, "./data/ls.dat") 
    
    call parlst_getvalue_double (rparams, 'GENERAL', 'kk', kk, 1.0_DP)  
    call parlst_getvalue_double (rparams, 'GENERAL', 'n', n, 1.0_DP)    
    
!     phi(x,y) = 8*(x*x+y*y)
!     f(x,y) = -32
!     Dcoefficients (1,:,:) = -32.0_DP

!     Dcoefficients (1,:,:) = 8.0_DP

!     phi(x,y) = 16*x*y*((1-x)*(1-y))
!     f(x,y) = 32*(y*(1-y)+x*(1-x))
!   !</subroutine>
!      Dcoefficients (1,:,:) = 32.0_DP * &
!                           ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
!                            Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )+ &
!                            16.0_DP * kk* Dpoints(1,:,:)* &
!                            (1.0_DP-Dpoints(1,:,:)) * &
!                            Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))

!     Dcoefficients (1,:,:) = 0.0_DP
  
!        phi(x,y) = SIN(PI * x) * SIN(PI * y)
!     => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)
!    Dcoefficients (1,:,:) = 2.0_DP * SYS_PI**2 &
!                         * sin(SYS_PI * Dpoints(1,:,:)) &
!                         * sin(SYS_PI * Dpoints(2,:,:))
	 
!    Dcoefficients (1,:,:) = (2.0_DP * SYS_PI**2 + kk)&
!                         * sin(SYS_PI * Dpoints(1,:,:)) &
!                         * sin(SYS_PI * Dpoints(2,:,:))
 
    Dcoefficients (1,:,:) = (2 * n**2 * SYS_PI**2)&
                         * (sin(n*SYS_PI * Dpoints(1,:,:))&
                         *  sin(n*SYS_PI * Dpoints(2,:,:)))


!        phi(x,y)   = SIN(x * (1-x)) * SIN(y * (1-y))
!     Dcoefficients (1,:,:) = -( ( ((-2*cos(Dpoints(1,:,:)*(1-Dpoints(1,:,:)))) + &
!     (1-2*Dpoints(1,:,:))**2 * (-sin(Dpoints(1,:,:)*(1-Dpoints(1,:,:))))) * & 
!     sin(Dpoints(2,:,:)*(1-Dpoints(2,:,:))) ) + ( ((-2*cos(Dpoints(2,:,:)* &
!     (1-Dpoints(2,:,:)))) + (1-2*Dpoints(2,:,:))**2 * (-sin(Dpoints(2,:,:)*&
!     (1-Dpoints(2,:,:))))) * sin(Dpoints(1,:,:)*(1-Dpoints(1,:,:))) ) )
    
  end subroutine

  ! ***************************************************************************
  

!<subroutine>

  subroutine coeff_RHS_alfa (rdiscretisation,rform, &
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
    ! DIMENSION(\#local DOF`s in test space,nelements)
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
  
  ! Local variable
    real(DP) :: kk, alfa
    ! parameter
    type(t_parlist) :: rparams  

    ! Ok, let us start. 
    ! reading data from the *.dat file 
    ! We want to solve our Poisson problem on level...
    call parlst_init(rparams)
    call parlst_readfromfile (rparams, "./data/ls.dat") 
    
    call parlst_getvalue_double (rparams, 'GENERAL', 'kk', kk, 1.0_DP)  
    
    alfa = 1.0_DP/kk
   
!     phi(x,y) = 16*x*y*((1-x)*(1-y))
!     f(x,y) = 32*(y*(1-y)+x*(1-x))
!   !</subroutine>
!      Dcoefficients (1,:,:) = 32.0_DP * &
!                      ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
!                        Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )

!     Dcoefficients (1,:,:) = 0.0_DP
			 
    Dcoefficients (1,:,:) = alfa*(2.0_DP * SYS_PI**2 + kk)&
                         * sin(SYS_PI * Dpoints(1,:,:)) &
                         * sin(SYS_PI * Dpoints(2,:,:))
    
  end subroutine
 ! ***************************************************************************  

  

!<subroutine>

  subroutine coeff_RHS_alfatu (rdiscretisation,rform, &
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
    ! DIMENSION(\#local DOF`s in test space,nelements)
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
  
  ! Local variable
    real(DP) :: kk, alfa
    ! parameter
    type(t_parlist) :: rparams  

    ! Ok, let us start. 
    ! reading data from the *.dat file 
    ! We want to solve our Poisson problem on level...
    call parlst_init(rparams)
    call parlst_readfromfile (rparams, "./data/ls.dat") 
    
    call parlst_getvalue_double (rparams, 'GENERAL', 'kk', kk, 1.0_DP)  
    
    alfa = kk
   
!     phi(x,y) = 16*x*y*((1-x)*(1-y))
!     f(x,y) = 32*(y*(1-y)+x*(1-x))
!   !</subroutine>
!      Dcoefficients (1,:,:) = 32.0_DP * &
!                      ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
!                        Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )

!     Dcoefficients (1,:,:) = 0.0_DP
			 
    Dcoefficients (1,:,:) = alfa*(2.0_DP * SYS_PI**2 + kk)&
                         * sin(SYS_PI * Dpoints(1,:,:)) &
                         * sin(SYS_PI * Dpoints(2,:,:))
    
  end subroutine
 ! ***************************************************************************  


!<subroutine>

  subroutine coeff_RHSS (rdiscretisation,rform, &
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
    ! DIMENSION(\#local DOF`s in test space,nelements)
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

  ! Local variable
    real(DP) :: kk
    ! parameter
    type(t_parlist) :: rparams  

    ! Ok, let us start. 
    ! reading data from the *.dat file 
    ! We want to solve our Poisson problem on level...
    call parlst_init(rparams)
    call parlst_readfromfile (rparams, "./data/ls.dat") 
    
    call parlst_getvalue_double (rparams, 'GENERAL', 'kk', kk, 1.0_DP) 

  
!        phi(x,y) = SIN(PI * x) * SIN(PI * y)
!     => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)coeff_RHS_2D
    Dcoefficients (1,:,:) = kk * (2.0_DP * SYS_PI**2 + kk) &
                         * sin(SYS_PI * Dpoints(1,:,:)) &
                         * sin(SYS_PI * Dpoints(2,:,:))

    
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

 ! Local variable
    real(DP) :: n
    ! parameter
    type(t_parlist) :: rparams 
  
!</subroutine>

   ! To get the X/Y-coordinates of the boundary point, use:
   REAL(DP) :: dx,dy
   CALL boundary_getCoords(rdiscretisation%p_rboundary, &
         rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
          

  
!   !     This is the 1st Edge
!     if ((dwhere.ge.0.0_DP).and.(dwhere.le.1.0_DP)) then
!     Dvalues(1) = 8.0_DP*dx*dx
!     end if
! !     2nd Edge
!     if ((dwhere.ge.1.0_DP).and.(dwhere.le.2.0_DP)) then
!     Dvalues(1) = 8.0_DP*(1.0_DP+dy*dy)
!     end if
! !     3rd Edge
!     if ((dwhere.ge.2.0_DP).and.(dwhere.le.3.0_DP)) then
!     Dvalues(1) = 8.0_DP*(1.0_DP+dx*dx)
!     end if
! !     4th Edge
!     if ((dwhere.ge.3.0_DP).and.(dwhere.le.4.0_DP)) then
!     Dvalues(1) = 8.0_DP*(dy*dy)
!     end if

!!=======================================================
!
!    ! Ok, let us start. 
!    ! reading data from the *.dat file 
!    ! We want to solve our Poisson problem on level...
!    call parlst_init(rparams)
!    call parlst_readfromfile (rparams, "./data/ls.dat") 
!    
!    call parlst_getvalue_double (rparams, 'GENERAL', 'n', n, 1.0_DP)
! !     This is the 1st Edge
!     if ((dwhere.ge.0.0_DP).and.(dwhere.le.1.0_DP)) then
!     Dvalues(1) = sin(n*SYS_PI*dx)
!     end if
! !     2nd Edge
!     if ((dwhere.ge.1.0_DP).and.(dwhere.le.2.0_DP)) then
!     Dvalues(1) = sin(n*SYS_PI)
!     end if
! !     3rd Edge
!     if ((dwhere.ge.2.0_DP).and.(dwhere.le.3.0_DP)) then
!     Dvalues(1) = sin(n*SYS_PI*dx)
!     end if
! !     4th Edge
!     if ((dwhere.ge.3.0_DP).and.(dwhere.le.4.0_DP)) then
!     Dvalues(1) = 0.0_DP
!     end if
!!=======================================================

  Dvalues(1) = 0.0_DP

  end subroutine
  ! ***************************************************************************

!<subroutine>

  subroutine coeff_mix (rdiscretisation,rform, &
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
    ! DIMENSION(\#local DOF`s in test space,nelements)
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
 
!     phi(x,y) = 8*(x*x+y*y)
!     f(x,y) = -32
!     Dcoefficients (1,:,:) = -32.0_DP

!     phi(x,y) = 16*x*y*((1-x)*(1-y))
!     f(x,y) = 32*(y*(1-y)+x*(1-x))
!   !</subroutine>
!      Dcoefficients (1,:,:) = 32.0_DP * &
!                      ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
!                        Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )

  ! Dcoefficients = 0.0_DP
  
    !    phi(x,y) = SIN(PI * x) * SIN(PI * y)
    ! => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)
    !Dcoefficients (1,:,:) = -2.0_DP * SYS_PI**2 &
    !                      * sin(SYS_PI * Dpoints(1,:,:)) &
    !                      * sin(SYS_PI * Dpoints(2,:,:))
    
    
    !    phi(x,y)   = SIN(x * (1-x)) * SIN(y * (1-y))
    Dcoefficients (1,:,:) = ( ( ((-2*cos(Dpoints(1,:,:)*(1-Dpoints(1,:,:)))) + &
    (1-2*Dpoints(1,:,:))**2 * (-sin(Dpoints(1,:,:)*(1-Dpoints(1,:,:))))) * & 
    sin(Dpoints(2,:,:)*(1-Dpoints(2,:,:))) ) + ( ((-2*cos(Dpoints(2,:,:)* &
    (1-Dpoints(2,:,:)))) + (1-2*Dpoints(2,:,:))**2 * (-sin(Dpoints(2,:,:)*&
    (1-Dpoints(2,:,:))))) * sin(Dpoints(1,:,:)*(1-Dpoints(1,:,:))) ) )
    
  end subroutine



  ! ***************************************************************************

!<subroutine>

    subroutine bdr_coeff_mat (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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
  
  ! local variables
  real(DP) :: dnx, dny
  integer :: i,j
  !</subroutine>
  
   do i = 1,nelements
    do j = 1,npointsPerElement
        call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary, ibct, DpointPar(j,i), dnx, dny) 
        Dcoefficients (1,:,:) = dnx**2
        Dcoefficients (2,:,:) = dnx*dny
        Dcoefficients (3,:,:) = dnx*dny
        Dcoefficients (4,:,:) = dny**2
    end do    
   end do
    
  end subroutine
  
! ***************************************************************************

!<subroutine>

subroutine bdr_coeff_vec (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
    
    ! Local variables
    integer :: i,j

    Dcoefficients(1,:,:) = 10.0_DP

!    do i = 1,nelements
!        do j = 1,npointsPerElement
!            if ((DpointPar(j,i) .gt. 0.0_DP).and.(DpointPar(j,i) .le. 0.5_DP)) then
!                Dcoefficients(1,j,i) = 1.0_DP !8.0_DP*Dpoints(1,j,i)*Dpoints(1,j,i)
!            end if
!            if ((DpointPar(j,i) .gt. 0.5_DP).and.(DpointPar(j,i) .lt. 1.0_DP)) then
!                Dcoefficients(1,j,i) = 4.0_DP !8.0_DP*Dpoints(1,j,i)*Dpoints(1,j,i)
!            end if
!            if ((DpointPar(j,i) .gt. 1.0_DP).and.(DpointPar(j,i) .lt. 2.0_DP)) then
!                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*(1.0_DP+Dpoints(2,j,i)*Dpoints(2,j,i))
!            end if
!            if ((DpointPar(j,i) .gt. 2.0_DP).and.(DpointPar(j,i) .lt. 3.0_DP)) then
!                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*(1.0_DP+Dpoints(1,j,i)*Dpoints(1,j,i))
!            end if
!            if ((DpointPar(j,i) .gt. 3.0_DP).and.(DpointPar(j,i) .lt. 4.0_DP)) then
!                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*Dpoints(2,j,i)*Dpoints(2,j,i)
!            end if            
!        end do
!    end do  

    end subroutine

! ***************************************************************************

!<subroutine>

subroutine bdr_coeff_vecn (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
  !</subroutine>
    
    ! Local variables
    integer :: i,j
    real(DP) :: dnx, dny
    
!    Dcoefficients(1,:,:) = 0.0_DP

    do i = 1,nelements
        do j = 1,npointsPerElement
            if ((DpointPar(j,i) .gt. 0.0_DP).and.(DpointPar(j,i) .lt. 1.0_DP)) then
                call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, DpointPar(j,i), dnx, dny)
                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*Dpoints(1,j,i)*Dpoints(1,j,i)*dnx
                Dcoefficients(2,j,i) = 0.0_DP !8.0_DP*Dpoints(1,j,i)*Dpoints(1,j,i)*dny             
            end if
            if ((DpointPar(j,i) .gt. 1.0_DP).and.(DpointPar(j,i) .lt. 2.0_DP)) then
                call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, DpointPar(j,i), dnx, dny)
                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*(1.0_DP+Dpoints(2,j,i)*Dpoints(2,j,i))*dnx
                Dcoefficients(2,j,i) = 0.0_DP !8.0_DP*(1.0_DP+Dpoints(2,j,i)*Dpoints(2,j,i))*dny
            end if
            if ((DpointPar(j,i) .gt. 2.0_DP).and.(DpointPar(j,i) .lt. 3.0_DP)) then
                call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, DpointPar(j,i), dnx, dny)
                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*(1.0_DP+Dpoints(1,j,i)*Dpoints(1,j,i))*dnx
                Dcoefficients(2,j,i) = 0.0_DP !8.0_DP*(1.0_DP+Dpoints(1,j,i)*Dpoints(1,j,i))*dny
            end if
            if ((DpointPar(j,i) .gt. 3.0_DP).and.(DpointPar(j,i) .lt. 4.0_DP)) then
                call boundary_getNormalVec2D(rdiscretisation%p_rboundary, ibct, DpointPar(j,i), dnx, dny)
                Dcoefficients(1,j,i) = 0.0_DP !8.0_DP*Dpoints(2,j,i)*Dpoints(2,j,i)*dnx
                Dcoefficients(2,j,i) = 0.0_DP !8.0_DP*Dpoints(2,j,i)*Dpoints(2,j,i)*dny
            end if            
        end do
    end do  

    end subroutine

! ***************************************************************************

!<subroutine>

    subroutine bdr_coeff_matcurl_a (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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
  
  ! local variables
  real(DP) :: dnx, dny
  integer :: i,j
  !</subroutine>
  
   do i = 1,nelements
    do j = 1,npointsPerElement
        call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary, ibct, DpointPar(j,i), dnx, dny) 
        Dcoefficients (1,:,:) = -dnx*dny
    end do    
   end do
    
  end subroutine


! ***************************************************************************

!<subroutine>

    subroutine bdr_coeff_matcurl_b (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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
  
  ! local variables
  real(DP) :: dnx, dny
  integer :: i,j
  !</subroutine>
  
   do i = 1,nelements
    do j = 1,npointsPerElement
        call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary, ibct, DpointPar(j,i), dnx, dny) 
        Dcoefficients (1,:,:) = dny*dny
    end do    
   end do
    
  end subroutine
    
! ***************************************************************************

!<subroutine>

    subroutine bdr_coeff_matcurl_c (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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
  
  ! local variables
  real(DP) :: dnx, dny
  integer :: i,j
  !</subroutine>
  
   do i = 1,nelements
    do j = 1,npointsPerElement
        call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary, ibct, DpointPar(j,i), dnx, dny) 
        Dcoefficients (1,:,:) = dnx*dnx
    end do    
   end do
    
  end subroutine


! ***************************************************************************

!<subroutine>

    subroutine bdr_coeff_matcurl_d (rdiscretisationTrial,&
                  rdiscretisationTest, rform, nelements, npointsPerElement,&
                  Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
                  rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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
  
  ! local variables
  real(DP) :: dnx, dny
  integer :: i,j
  !</subroutine>
  
   do i = 1,nelements
    do j = 1,npointsPerElement
        call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary, ibct, DpointPar(j,i), dnx, dny) 
        Dcoefficients (1,:,:) = -dnx*dny
    end do    
   end do
    
  end subroutine


! ***************************************************************************

!<subroutine>

subroutine bdr_coeff_veccurl (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, ibct, DpointPar, &
                  IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
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

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>

    Dcoefficients(1,:,:) = 0.0_DP
    
    end subroutine

end module

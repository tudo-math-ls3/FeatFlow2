!##############################################################################
!# ****************************************************************************
!# <name> chemotaxis_callback </name>
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
!# 1.2) coeff_chemo
!#     -> Returns the coefficients for the chemo matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
	
!# 2.) coeff_RHS
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getBoundaryValues
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

module chemotaxis_callback

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
	
    ! Defining pi
    real(DP), parameter, public :: PI = 3.141592654_DP

contains
	
! *******************************************************************************
!# functions to implement the chemoattractant influence to the actual chemotaxis model
               


! The following fct is to specify a certain fct  CHI which maybe nonlinear
! Furtheron it may depend on u and c , e.g. CHI = CHI (u,c)
!       u and c should already are fct evaluated at a certain point
!       a  represents a constant
! For this test case we consider a constant CHI. E.g. CHI = a
! For the cherkur3 model this should be set to  = 1.0_DP / ( 1.0_DP + c )**2
! this simulates signal-dependent sensitivity
        function f_CHI(u,c,a) result (chi_result)
                implicit none
                 real(DP) :: u,c,a,chi_result, ALPHA
                ALPHA = 1.0_DP
                 chi_result = a!1.0_DP / ( 1.0_DP + c )**2 !a / ( (1.0_DP + ALPHA*c )**2)
         end function f_CHI

! the folowing two fcts. are used for the timedependent initial sol
! E.g. chemotaxiscoupldRHStimeconsttest.f90
        function Fu(t) result (f_resultu)
                implicit none
                 real(DP) :: t,f_resultu
                 f_resultu = 2*dexp(t)
         end function Fu
                      

        ! This functional is corresponding to te generic functional D ( . ) in the Hillen, Painter paper
        ! It determines the (non-) constant coefficient of the diffusion-operator
        function D ( u, D_1, N ) result (f_result)
            implicit none
            real(dP) :: u, N, D_1, f_result
            f_result = D_1 !* u **N
        end function D

        ! Here we specify the operators which influce the chemotaxis-sensitivity function
        ! This determines the density-dependent sensitivity.
        ! If = 1 there's no strong density-dependency
        ! If = ( 1 - u/GAMMA ) it simulates, that "...occupation ( by cells )of an area limits other cells from penetrating it ..." ( Hillen, Painter paper )
        function A( u, GAMMA ) result (f_result)
            implicit none
            real(DP) :: u, GAMMA, f_result
            f_result = 1.0_DP!( 1.0_DP - u/GAMMA )
        end function A

        ! This functional determines a signal-dependent sensitivity. (ref. Hillen, Painter paper)
        ! If = CHI there's no signal dependencies, e.g. the cell density continues to grow.
        ! If = CHI / ( 1 + ALPHA* c )**2 (or likewise) it simulates the competition-condition of
        ! chemoattractant molecules binding on corresponding receptors on the cell surface.
        ! (which happens to invoke the chemotaxis)
        function B( c, CHI, ALPHA ) result (f_result)
            implicit none
            real(DP) :: c, CHI, ALPHA, f_result
            f_result = CHI!1.0_DP!CHI / ( 1.0_DP + ALPHA* c )**2
        end function B

        ! This functional is corresponding to te generic functional f ( . ) in the Hillen, Painter paper
        ! It determines the kinetics of the cells. If = 0 then there's no kinetics
        function f ( u, R ) result (f_result)
            implicit none
            real(dP) :: u, R, f_result
            f_result = u * ( 1.0_DP - u )
        end function f

        ! This function is corresponding to the generic functional g ( . ) in the Hillen, Painter
        !paper
        ! This determines the saturating signal production rate. If =1 cells will continue
        !producing chemoattractants, even if the cell population blows-up.
        ! If = 1 / ( 1 + PHI * u ) "this would prevent excessive chemoattractant production as
        !the cell density increases" ( Hillen, Painter paper )
        function func_g ( u, PHI ) result (f_result)
            implicit none
            real(DP) :: u, PHI, f_result
            f_result = 1.0_DP ! / ( 1.0_DP + PHI * u )
        end function func_g

                      
        ! This should realize a small pertubation at the center of the domain
        function ic_pattern (x,y) result(f_result)
            implicit none
            intrinsic RANDOM_NUMBER
            real(DP) :: x, y, f_result, random_num
            if ( sqrt ( (x-8)**2 + (y-8)**2) <= 1.5_DP ) then
                CALL RANDOM_NUMBER (random_num)
                !f_result = random_num
		        f_result = 0.2_DP
                !f_result = rand(0) !1.1 * cos ( 4 * ( PI * sqrt ( (x-8)**2 + (y-8)**2 ) ) / 4 ) **2
            else
                f_result = 0.0_DP
            end if
        end function ic_pattern


        ! This function is used to describe a pattern forming IC for u
        function ic_cos( x, y ) result (f_result)
            implicit none
            real(DP) :: x, y, f_result
            if ( sqrt ( (x-5)**2 + (y-5)**2) <=4 ) then
                f_result = 5.0_DP * cos ( 2 * ( PI * sqrt ( (x-5)**2 + (y-5)**2 ) ) / 4 ) **2
            else
                f_result = 0.0_DP
            end if
        end function

        function Fc(t) result (f_resultc)
                implicit none
                real(DP) :: t,f_resultc
                f_resultc = dexp(t)
         end function Fc
               
        function f_sigma(x,y) result (f_result)
                implicit none
                intrinsic RANDOM_NUMBER
                real(DP) :: x, y, f_result, random_num
                CALL RANDOM_NUMBER (random_num)
                f_result = random_num
                !f_result = rand( 0 )
        end function f_sigma

! the following two fcts. are used for the nontrivial initial sol
! E.g. chemotaxiscoupldRHSnontrivialtest.f90
            function U_nontrivial(x,y,t) result (f_resultu)
                implicit none
                 real(DP) :: x,y,t,f_resultu
                 f_resultu = 2*C_nontrivial(x,y,t)
         end function U_nontrivial


        function C_nontrivial(x,y,t) result (f_resultc)
                implicit none
                real(DP) :: x,y,t,f_resultc
                f_resultc = (x-1)*x*(y-1)*y*dexp(t)
         end function C_nontrivial
               
! the following two fcts. are used for the nonconst initial sol
! E.g. chemotaxiscoupldRHStest.f90
        function U_nonconst(a,x,y) result (f_resultu)
                        implicit none
                        real(DP) :: a,x,y,f_resultu
                        f_resultu = a*(x-1)*x*y*(y-1)
                end function U_nonconst


        function C_nonconst(a,x,y) result (f_resultc)
                        implicit none
                        real(DP) :: a,x,y,f_resultc
                        f_resultc = a*(x-1)*x*y*(y-1)
                end function C_nonconst
                
! The following function is needed to asseble the RHS f in the
! u-part. See coeff_cherkur3_f for further details
! f_chemo  = ALPHA * grad * ( u / ( 1+c )^2 * grad(c) )
!               = ALPHA* { (u / ( 1+c )^2 * c_x )_x + (u / ( 1+c )^2 * c_y )_y }
        function f_chemo( x,y ) result (f_result_chemo)
                    implicit none
                    real(DP) :: x, y, f_result_chemo
                    f_result_chemo = 2* ( f_1stpart_chemo ( x,y ) + f_2ndpart_chemo ( x,y ) )
        end function f_chemo

! Deriving  ( u / ( 1+c )^2 * c_x )_x
        function f_1stpart_chemo( x,y ) result (f_result)
                    implicit none
                    real(DP) :: x, y, f_result
                    f_result = ( (2*(2*x-1)*(y**2-y) * (1+C_nonconst(1.0_DP,x,y))**2 - 2* U_nonconst(1.0_DP,x,y)&
                              *(1+C_nonconst(1.0_DP,x,y)) * ( 2*x-1 ) * ( y**2-y ) ) / ( 1+C_nonconst(1.0_DP,x,y) )**4 )&
                                *( 2*x-1 ) * ( y**2-y )+&
                                (U_nonconst(1.0_DP,x,y) / ( 1+C_nonconst(1.0_DP,x,y) )**2) * 2*(y**2-y)
        end function f_1stpart_chemo

! Deriving  ( u / ( 1+c )^2 * c_y )_y
        function f_2ndpart_chemo( x,y ) result (f_result)
                    implicit none
                    real(DP) :: x, y, f_result
                    f_result = ( (2*(2*y-1)*(x**2-x) * (1+C_nonconst(1.0_DP,x,y))**2 - 2* U_nonconst(1.0_DP,x,y)&
                              *(1+C_nonconst(1.0_DP,x,y)) * ( 2*y-1 ) * ( x**2-x ) ) / ( 1+C_nonconst(1.0_DP,x,y) )**4 )&
                                *( 2*y-1 ) * ( x**2-x )+&
                                (U_nonconst(1.0_DP,x,y) / ( 1+C_nonconst(1.0_DP,x,y) )**2) * 2*(x**2-x)
        end function f_2ndpart_chemo


	function c(x,y,a,b) result (c_result1)
		implicit none
		real(DP) :: c_result1, x, y, a, b

		c_result1 = a*dexp(-b*(x*x+y*y))
	end function c
	
	function c_x(x,y,a,b) result (c_result2)
		implicit none
		real(DP) :: c_result2, x, y, a, b

		c_result2 = -2*b*(x-0.5)*c(x,y, a, b)
	end function c_x

	function c_xx(x,y,a,b) result (c_result3)
		implicit none
		real(DP) :: c_result3, x, y, a, b

		c_result3 = -2*b*c(x,y,a,b)+4*b**2*(x-0.5)**2*c(x,y, a, b)
	end function c_xx


	function c_y(x,y,a,b) result (c_result4)
		implicit none
		real(DP) :: c_result4, x, y, a, b

		c_result4 = -2*b*(y-0.5)*c(x,y, a, b)
	end function c_y
	


	function c_yy(x,y,a,b) result (c_result5)
		implicit none
		real(DP) :: c_result5, x, y, a, b

		c_result5 = -2*b*c(x,y, a, b)+4*b**2*(y-0.5)**2*c(x,y, a, b)
	end function c_yy

  ! ***************************************************************************
  !<subroutine>

  subroutine coeff_chemotaxis (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)						!--->dtstep added
    
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
	
    real(DP) :: amp = 1.0_DP   ! Amplifies the chemotaxis
    integer :: i,j
    real(DP) :: a,b
    ! Determines whether or not we use an implicit scheme, e.g. =1 for impli
    !  =-1 for explicit
    integer :: impli = 1
    
    real(DP) :: dtimestep
    dtimestep = rcollection%Dquickaccess(1)

    a = 600
    b = 60
	DO i=lbound(Dpoints,2),ubound(Dpoints,2)						!--->return value modified
		DO j=lbound(Dpoints,3),ubound(Dpoints,3)
			! Berechne d_xx c + d_yy c an der entsprechenden Stelle (d.h. Pkt)
    			Dcoefficients(1,i,j) = impli * dtimestep * amp *( c_xx( Dpoints(1,i,j),Dpoints(2,i,j), a , b )&
				 + c_yy(Dpoints(1,i,j),Dpoints(2,i,j), a , b) )
			! Berechne d_x c an der entsprechenden Stelle (d.h. Pkt)
			Dcoefficients(2,i,j) = impli * dtimestep * amp *( c_x(Dpoints(1,i,j) , Dpoints(2,i,j), a , b) )
			! Berechne d_y c an der entsprechenden Stelle (d.h. Pkt)
			Dcoefficients(3,i,j) = impli * dtimestep * amp *c_y(Dpoints(1,i,j) , Dpoints(2,i,j), a , b)
		END DO
	END DO

  end subroutine

  ! ***************************************************************************



  ! ***************************************************************************
  !<subroutine>

 subroutine coeff_heatcond (rdiscretisationTrial,rdiscretisationTest,rform, &
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

  subroutine coeff_RHS (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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

  ! First of all get the vector saved in the collection to obtain the linform factor
  type(t_vectorScalar),pointer :: uvector
  logical :: bexists
  if(present(rcollection))then
  uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
  IF(.NOT.bexists)then
	call output_lbrk ()
        call output_line ("***************ERROR:***************")
        call output_line ("**********COLLECTION FAILED*********")
        call output_lbrk ()
  END IF
  ELSE
  Dcoefficients = 0.0_DP
      
END IF
! assign the corresponding vectorentry of uvector to the needed callback coefficient
   
  call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
! 	uvector => null()		# to release the pointer !?!?
	
  end subroutine
  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSu (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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

  ! First of all get the vector saved in the collection to obtain the linform factor
  type(t_vectorScalar),pointer :: uvector
  logical :: bexists
! assign the corresponding vectorentry of uvector to the needed callback coefficient
      if(present(rcollection))then
             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
            IF(.NOT.bexists)then
	           call output_lbrk ()
                   call output_line ("***************ERROR:***************")
                   call output_line ("**********COLLECTION FAILED*********")
                   call output_lbrk ()
            END IF
            call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
      else
            Dcoefficients = 0.002_DP ! = BETA*dtime*u
            print *,"**************SHOULD NOT APPEAR***************"
      end if


  end subroutine
  ! ***************************************************************************


  ! ***************************************************************************
!<subroutine>
! this subroutine is used for nonconstant test fcts
! E.g. u=2*c=2*x(x-1)y(y-1).
  subroutine coeff_RHSf (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
  type(t_vectorScalar),pointer :: uvector
  logical :: bexists

    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, SCALE_C, SCALE_U
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)
    D_1 = rcollection%DquickAccess(8)
    D_2 = rcollection%DquickAccess(9)
    CHI = rcollection%DquickAccess(7)
    SCALE_C = rcollection%DquickAccess(10)
    SCALE_U = rcollection%DquickAccess(11)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!       else
    ! this represents the fct. f for 2c=u=2*(x-1)x(y-1)y
    ! returning the analytical computed RHS for the  solution (e.g. cell density )PDE
    !-dtstep * { 4 ( y(y-1)+x(x-1) )  +  CHI [ (x²-x)² (6y²-6y+1)+(y²-y)² (6x²-6x+1 ) ] }
    Dcoefficients (1,:,:) = dtstep*( -2.0_DP * SCALE_U *D_1*( Dpoints(2,:,:)*( Dpoints(2,:,:)-1.0_DP)+Dpoints(1,:,:)*&
                                                ( Dpoints(1,:,:)-1.0_DP) )&
                                        + CHI * SCALE_C * SCALE_U * ( ( ( Dpoints(1,:,:)**2 - Dpoints(1,:,:) )**2) *&
                                                (6 * Dpoints(2,:,:)**2-6 * Dpoints(2,:,:)+1.0_DP) +&
                                                (( Dpoints(2,:,:)**2 - Dpoints(2,:,:))**2) *&
                                                (6 * Dpoints(1,:,:)**2-6 * Dpoints(1,:,:)+1.0_DP)))


!   saved
!     Dcoefficients (1,:,:) = -dtstep *  ( 4.0_DP * D_1 * ( Dpoints(2,:,:) * ( Dpoints(2,:,:)-1.0_DP) + Dpoints(1,:,:) * &
!                                     ( Dpoints(1,:,:)-1.0_DP) ) + CHI * ( ( ( Dpoints(1,:,:)**2 - Dpoints(1,:,:) )**2) *&
!  (6 * Dpoints(2,:,:)**2-6 * Dpoints(2,:,:)+1.0_DP) +  (( Dpoints(2,:,:)**2 - Dpoints(2,:,:))**2) *&
!  (6 * Dpoints(1,:,:)**2-6 * Dpoints(1,:,:)+1.0_DP)))
      

  end subroutine


  ! ***************************************************************************

!<subroutine>
! this subroutine is used for the heatconduction explicit euler method
  subroutine coeff_RHSee (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional : A collection structure to provide additional
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

    Dcoefficients (1,:,:) =  ( 2.0_DP  * ( Dpoints(2,:,:) * ( Dpoints(2,:,:)-1.0_DP) + Dpoints(1,:,:) * &
                                    ( Dpoints(1,:,:)-1.0_DP) ))


      

  end subroutine


  ! ***************************************************************************


!<subroutine>
! this subroutine is used for nonconstant test fcts
! E.g. u=2*c=2*x(x-1)y(y-1).
  subroutine coeff_RHSg (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    type(t_vectorScalar),pointer :: uvector
    logical :: bexists
    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, SCALE_C , SCALE_U
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)
    D_1 = rcollection%DquickAccess(8)
    D_2 = rcollection%DquickAccess(9)
    CHI = rcollection%DquickAccess(7)
    SCALE_C = rcollection%DquickAccess(10)
    SCALE_U = rcollection%DquickAccess(11)
!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! should rpresent the rhs fct g since 2*c=u=....
    ! returning the analytical computed RHS for the  chemoattractant PDE

           Dcoefficients(1,:,:) = dtstep*    (-2.0_DP * SCALE_C * D_2 * (  Dpoints(2,:,:)*( Dpoints(2,:,:)-&
                                            1.0_DP) + Dpoints(1,:,:) * ( Dpoints(1,:,:) - 1.0_DP ) ) &
                                            - SCALE_U * BETA * ( Dpoints(2,:,:) * ( Dpoints(2,:,:)- 1.0_DP) *&
                                            Dpoints(1,:,:) * ( Dpoints(1,:,:)- 1.0_DP  ) ) + &
                                            SCALE_C * SCALE_U * ALPHA * ((Dpoints(2,:,:)**2) * ( ( Dpoints(2,:,:)-1.0_DP)**2 ) *&
                                            (Dpoints(1,:,:)** 2)  * ( Dpoints(1,:,:)-1.0_DP) **2))





  end subroutine

  ! ***************************************************************************


 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_poly (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! local variables given by the collection
    real(DP) :: SCALE_X

    ! loop variables
    integer :: i,j

    SCALE_X = rcollection%DquickAccess(1)


!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned

!     DO i =1,nelements
!         DO j=1,npointsPerElement
!             Dvalues(i,j) =  U_nonconst ( SCALE_X , Dpoints(1,i,j) , Dpoints(2,i,j) )
!         END DO
!     END DO

    Dvalues(:,:) =  SCALE_X* ( Dpoints(1,:,:)**2 -  Dpoints(1,:,:) ) * ( Dpoints(2,:,:)**2&
 - Dpoints(2,:,:) )
  end subroutine

  ! ***************************************************************************


 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_const (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! local variables given by the collection
    real(DP) :: CONST

    CONST = rcollection%DquickAccess(1)

    Dvalues(:,:) =  CONST
  end subroutine

  ! ***************************************************************************


 ! ***************************************************************************


    !<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_pattern (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! loop-indices
    integer :: icub, iel

    DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            Dvalues( icub, iel ) = 1_DP
            !userPresc_cellsInitCond(Dpoints ( 1, icub, iel ), &
            !                                               Dpoints ( 2, icub, iel ), &
            !                                               Dpoints ( 1, icub, iel ))
        END DO
    END DO

  end subroutine

  ! ***************************************************************************


 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_cos (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! loop-indices
    integer :: icub, iel

    DO iel = 1, nelements
        DO icub = 1, npointsPerElement
            Dvalues( icub, iel ) =  ic_cos  ( Dpoints ( 1, icub, iel ), Dpoints ( 2, icub, iel ) )
        END DO
    END DO

  end subroutine

  ! ***************************************************************************


  ! ***************************************************************************


!<subroutine>
! this subroutine is used for constant test fcts.
! E.g. u=2*c= 2*1
  subroutine coeff_RHSconst (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    type(t_vectorScalar),pointer :: uvector
    logical :: bexists
    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, U_0, C_0
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)
    D_1 = rcollection%DquickAccess(4)
    D_2 = rcollection%DquickAccess(5)
    CHI = rcollection%DquickAccess(6)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! should rpresent the rhs fct g since 2*c=u=....
    ! returning the analytical computed RHS for the  chemoattractant PDE

           Dcoefficients(1,:,:) = dtstep*  ( ALPHA*2.0_DP -BETA*2.0_DP)





  end subroutine

  ! ***************************************************************************

  ! ***************************************************************************

!<subroutine>
! this subroutine is used for the timedependent solution
! E.g. u=2*c=2*exp(-t)
  subroutine coeff_RHS_timedepend_g (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    type(t_vectorScalar),pointer :: uvector
    logical :: bexists
    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, STEPS,STARTTIME,TEND,DTIME
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    DTIME = rcollection%DquickAccess(6)
    TEND = STEPS*dtstep + STARTTIME

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! should rpresent the rhs fct g. E.g. g= c_t -D_2*Laplacian c -BETA*u + ALPHA*u*c
    ! returning the analytical computed RHS for the  chemoattractant PDE

           Dcoefficients(1,:,:) = &
                            dtstep*  ( Fc(DTIME)-BETA*Fu(DTIME)+ALPHA*Fu(DTIME)*Fc(DTIME))





  end subroutine

  ! ***************************************************************************

  ! ***************************************************************************

!<subroutine>
! this subroutine is used for the timedependent solution
! E.g. u=2*c=2*exp(-t)
  subroutine coeff_RHS_nontrivial_g (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    type(t_vectorScalar),pointer :: uvector
    logical :: bexists
    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, STEPS,STARTTIME,TEND,DTIME
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    DTIME = rcollection%DquickAccess(6)
    D_2 = rcollection%DquickAccess(9)
    TEND = STEPS*dtstep + STARTTIME

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! should rpresent the rhs fct g. E.g. g= c_t -D_2*Laplacian c -BETA*u + ALPHA*u*c
    ! returning the analytical computed RHS for the  chemoattractant PDE

           Dcoefficients(1,:,:) = dtstep* (Dpoints(2,:,:)*( Dpoints(2,:,:)-1.0_DP) * Dpoints(1,:,:) *&
                                                     ( Dpoints(1,:,:) - 1.0_DP )*dexp(DTIME) &
                                        - 2.0_DP *D_2 *dexp(DTIME)* (  Dpoints(2,:,:)*( Dpoints(2,:,:)-&
                                            1.0_DP) + Dpoints(1,:,:) * ( Dpoints(1,:,:) - 1.0_DP ) ) &
                                        - BETA *2.0_DP* (dexp(DTIME)* Dpoints(2,:,:) * ( Dpoints(2,:,:)- 1.0_DP) *&
                                            Dpoints(1,:,:) * ( Dpoints(1,:,:)- 1.0_DP  ) ) &
                                        + ALPHA *2.0_DP*dexp(2*DTIME)* ((Dpoints(2,:,:)**2) * ( ( Dpoints(2,:,:)-1.0_DP)**2 ) *&
                                            (Dpoints(1,:,:)** 2)  * ( Dpoints(1,:,:)-1.0_DP) **2))


  end subroutine

 ! ***************************************************************************

  ! ***************************************************************************

!<subroutine>
! this subroutine is used for the nontrivial reference  solution
! E.g. u=2*c=2*(x-1)*x*(y-1)*y*exp(t)
  subroutine coeff_RHS_nontrivial_f (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    type(t_vectorScalar),pointer :: uvector
    logical :: bexists
    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, STEPS,STARTTIME,TEND,DTIME
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    DTIME = rcollection%DquickAccess(6)
    CHI = rcollection%DquickAccess(7)
    D_1 = rcollection%DquickAccess(8)
    TEND = STEPS*dtstep + STARTTIME

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! this represents the fct. f for u=2c=2*(x-1)x(y-1)yexp(t)
    ! returning the analytical computed RHS for the  solution (e.g. cell density )PDE
    !-dtstep * { 4 ( y(y-1)+x(x-1) )  +  CHI [ (x²-x)² (6y²-6y+1)+(y²-y)² (6x²-6x+1 ) ] }
    Dcoefficients (1,:,:) = dtstep *  (2*Dpoints(2,:,:)*( Dpoints(2,:,:)-1.0_DP) *&
                                                 Dpoints(1,:,:) * ( Dpoints(1,:,:) - 1.0_DP )*dexp(DTIME) &
                                    -4.0_DP * D_1 * dexp(DTIME)*&
                                        ( Dpoints(2,:,:) * ( Dpoints(2,:,:)-1.0_DP) + Dpoints(1,:,:) * &
                                        ( Dpoints(1,:,:)-1.0_DP) )&
                                    + dexp(2*DTIME)*CHI * 2.0_DP * &
                                        ( ( ( Dpoints(1,:,:)**2 - Dpoints(1,:,:) )**2) * (6 * Dpoints(2,:,:)**2-&
                                        6 * Dpoints(2,:,:)+1.0_DP) +  (( Dpoints(2,:,:)**2 - Dpoints(2,:,:))**2) *&
                                        (6 * Dpoints(1,:,:)**2-6 * Dpoints(1,:,:)+1.0_DP)))




  end subroutine
  ! ***************************************************************************

 ! ***************************************************************************


!<subroutine>
! this routine is used for timdependent sol.
! E.g. u=2*c=2*exp(-t)
  subroutine coeff_RHS_timedepend_f (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    type(t_vectorScalar),pointer :: uvector
    logical :: bexists
    real(DP) :: dtstep
    real(DP) :: CHI, D_1, D_2, ALPHA, BETA, STEPS,STARTTIME, TEND, DTIME
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    DTIME = rcollection%DquickAccess(6)
    TEND = STEPS*dtstep + STARTTIME

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! should rpresent the rhs fct f. E.g. f = u_t - D_1 Laplacian u + CHI* grad dot ( u* grad c )
    ! returning the analytical computed RHS for the  chemoattractant PDE

           Dcoefficients(1,:,:) =  dtstep* Fu(DTIME)





  end subroutine

  ! ***************************************************************************

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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

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
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return zero Dirichlet boundary values for all situations.
  Dvalues(1) = 0.0_DP
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine

  ! ***************************************************************************

  ! ***************************************************************************

!<subroutine>

! is used to set the BC for the constant sol.
! To achieve U=2*c=2*1 we have certain constant BCs
  subroutine getBoundaryValuesU (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  

!</subroutine>


!   ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 2 Dirichlet boundary values for all situations.
  Dvalues(1) = 2.0_DP
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine

  ! ***************************************************************************


  ! ***************************************************************************

!<subroutine>
! is used to set the BC for the constant sol.
! To achieve U=2*c=2*1 we have certain constant BCs
  subroutine getBoundaryValuesC (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

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
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 2 Dirichlet boundary values for all situations.
  Dvalues(1) = 1.0_DP
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine

  ! ***************************************************************************

!   ! ***************************************************************************
! !<subroutine>
! Just messing with the BC. Since we need 'em, but it's not clear
! wether we use time-dependent ( ! ) DBC or NBC instead.
  subroutine getBoundaryValuestimeU (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
    real(DP) :: dtstep
    real(DP) :: STEPS,STARTTIME, TEND, dtime


!</subroutine>
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    TEND = STEPS*dtstep + STARTTIME
    dtime = rcollection%DquickAccess(6)
  ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 1 Dirichlet boundary values for all situations.
  Dvalues(1) = Fu(DTIME)
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine



  ! ***************************************************************************
!<subroutine>
! Just messing with the BC. Since we need 'em, but it's not clear
! wether we use time-dependent ( ! ) DBC or NBC instead.
  subroutine getBoundaryValuestimeC (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
    real(DP) :: dtstep
    real(DP) :: STEPS,STARTTIME, TEND, dtime

!</subroutine>

! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    TEND = STEPS*dtstep + STARTTIME
    dtime = rcollection%DquickAccess(6)



  ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 1 Dirichlet boundary values for all situations.
  Dvalues(1) = Fc(DTIME)
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine
 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
! Just messing with the BC. Since we need 'em, but it's not clear
! wether we use time-dependent ( ! ) DBC or NBC instead.
  subroutine getBoundaryValues_nontrivial (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
    real(DP) :: dtstep
    real(DP) :: STEPS,STARTTIME, TEND, dtime

!</subroutine>

! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    STEPS = rcollection%DquickAccess(4)
    STARTTIME = rcollection%DquickAccess(5)
    TEND = STEPS*dtstep + STARTTIME
    dtime = rcollection%DquickAccess(6)



  ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 1 Dirichlet boundary values for all situations.
  Dvalues(1) = 0.0_DP
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine
 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
! Just messing with the BC. Since we need 'em, but it's not clear
! wether we use time-dependent ( ! ) DBC or NBC instead.
  subroutine getBoundaryValues_constC (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
    real(DP) :: C_0, U_0

!</subroutine>

! assign the corresponding vectorentry of uvector to the needed callback coefficient

    C_0 = rcollection%DquickAccess(1)
    U_0 = rcollection%DquickAccess(2)



  ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 1 Dirichlet boundary values for all situations.
  Dvalues(1) = C_0
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine
 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
! Just messing with the BC. Since we need 'em, but it's not clear
! wether we use time-dependent ( ! ) DBC or NBC instead.
  subroutine getBoundaryValues_constU (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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
  integer, dimension(:), intent(IN)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer(I32), intent(IN)                                    :: ielement
  
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
  integer(I32), intent(IN)                                     :: iwhere

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
  type(t_collection), intent(IN), optional      :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(OUT)                         :: Dvalues
!</output>
  
    real(DP) :: C_0, U_0

!</subroutine>

! assign the corresponding vectorentry of uvector to the needed callback coefficient

    C_0 = rcollection%DquickAccess(1)
    U_0 = rcollection%DquickAccess(2)



  ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! REAL(DP) :: dx,dy
  !
  ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return 1 Dirichlet boundary values for all situations.
  Dvalues(1) = U_0
  
  !if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) &
  !  Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine
 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
    ! This subroutine is used by chemotaxis_cherkur_3
    ! Here we calculate  the RHS fct. of the c-part.
    ! Since the RHS depends on the cell density u
    ! we need a collection structure and the fevl routine
    ! to calculate u at the cubaturepts.
  subroutine coeff_hillenM6_RHS_c(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: PHI

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector

    ! Setting the params contained in the collection
    PHI = rcollection%DquickAccess(1)

    ! Fetching the vector
!    rvector = rcollection%p_rvectorQuickAccess1
    rvector = collct_getvalue_vecsca (rcollection, "cbvector",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term w*u_n^2/(SIGMA+u_n^2) for the RHS of c of the third
    ! chertock kurganov example
    Dcoefficients(1,:,:) = DvaluesFevl(1,:,:) / ( 1 + PHI * DvaluesFevl(1,:,:) )

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************



 ! ***************************************************************************
!<subroutine>
    ! This routine is used for a generic chemoattractant equation, which is presented by Hillen and Painter
    ! It specifies the dependency of the RHS to the cell distribution
    ! The chemoattractant equation reads:
    ! c_t = \Delta c + u*g(u) - c
    ! So this callback function returns u*g(u) used by a linearform.
  subroutine coeff_hillenX_RHS_c(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: PHI
    integer :: icub, iel

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector

    ! Setting the params contained in the collection
    PHI = rcollection%DquickAccess(1)

    ! Fetching the vector
!    rvector = rcollection%p_rvectorQuickAccess1
    rvector = collct_getvalue_vecsca (rcollection, "cbvector",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    DO icub = 1, npointsPerElement
        DO iel = 1, nelements
            ! function func_g returns 1
            Dcoefficients(1,icub,iel) = DvaluesFevl(1,icub,iel) *  func_g ( DvaluesFevl(1,icub,iel), PHI )
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************

 ! ***************************************************************************

! This callback routine is used for test-purpose. It simulates the source term f in the pattern model
! as a linear-term. So it designs the non-linearity into a linearity, which normally comes up with
!  accuracy-pbs !!!!!!!!!!!!
  subroutine coeff_pattern_RHS_u(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: PHI
    integer :: icub, iel

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rcell


    ! Fetching the vector
!    rvector = rcollection%p_rvectorQuickAccess1
    rcell = collct_getvalue_vecsca (rcollection, "cbvector",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rcell, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    DO icub = 1, npointsPerElement
        DO iel = 1, nelements
            Dcoefficients(1,icub,iel) = DvaluesFevl(1,icub,iel)**2 * ( 1.0_DP - DvaluesFevl(1,icub,iel)  )
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************




 ! ***************************************************************************
!<subroutine>
    ! This subroutine is used by chemotaxis_cherkur_3
    ! Here we calculate  the RHS fct. of the c-part.
    ! Since the RHS depends on the cell density u
    ! we need a collection structure and the fevl routine
    ! to calculate u at the cubaturepts.
  subroutine coeff_hillen_RHS_c(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: PHI

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector

    ! Setting the params contained in the collection
    PHI = rcollection%DquickAccess(1)

    ! Fetching the vector
!    rvector = rcollection%p_rvectorQuickAccess1
    rvector = collct_getvalue_vecsca (rcollection, "cbvector",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term w*u_n^2/(SIGMA+u_n^2) for the RHS of c of the third
    ! chertock kurganov example
    Dcoefficients(1,:,:) = DvaluesFevl(1,:,:)

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************




 ! ***************************************************************************
!<subroutine>
    ! This subroutine is used by chemotaxis_cherkur_3
    ! Here we calculate  the RHS fct. of the c-part.
    ! Since the RHS depends on the cell density u
    ! we need a collection structure and the fevl routine
    ! to calculate u at the cubaturepts.
  subroutine coeff_cherkur3_RHS_c(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: W, SIGMA, BETA

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector

    ! Setting the params contained in the collection
    W = rcollection%DquickAccess(1)
    SIGMA = rcollection%DquickAccess(2)
    BETA = rcollection%DquickAccess(3)

    ! Fetching the vector
!    rvector = rcollection%p_rvectorQuickAccess1
    rvector = collct_getvalue_vecsca (rcollection, "cbvector",0,'')

    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term w*u_n^2/(SIGMA+u_n^2) for the RHS of c of the third
    ! chertock kurganov example
    Dcoefficients(1,:,:) = BETA*( W*DvaluesFevl(1,:,:)*DvaluesFevl(1,:,:) ) / ( SIGMA+DvaluesFevl(1,:,:)*DvaluesFevl(1,:,:) )

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************




 ! ***************************************************************************
!<subroutine>
    ! This subroutine is used by chemotaxis_cherkur_3
    ! Here we calculate  the RHS fct. of the c-part.
    ! Since the RHS depends on the cell density u
    ! we need a collection structure and the fevl routine
    ! to calculate u at the cubaturepts.
  subroutine coeff_marrocco_RHS_c(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: ALPHA

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:) , allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector, rvector_star

    ! Setting the params contained in the collection
    ALPHA = rcollection%DquickAccess(1)

    ! Fetching the vector
!    rvector = rcollection%p_rvectorQuickAccess1
    rvector = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_star = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    ! Alllocate some memory for the array, since it'll be written by
    ! the fevl calls
    allocate(DvaluesFevl(1,npointsPerElement , nelements))

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_star, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 2)

    ! calculate the term w*u_n^2/(SIGMA+u_n^2) for the RHS of c of the third
    ! chertock kurganov example
    Dcoefficients(1,:,:) = ALPHA*DvaluesFevl(1,:,:)*dexp( - DvaluesFevl(1,:,:) / DvaluesFevl(2,:,:) )

    deallocate(DvaluesFevl)
  end subroutine

 ! ***************************************************************************



 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_cherkur_3.f90.
    ! calculates the chemosensitivity function, since in this model  it's no more
    ! constant, but nonlinear.
    ! We'll iterate once, so there's no loop and no residual-check.
  subroutine coeff_cherkur3_u(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c , rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, ALPHA

    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector
!     rvector_c => rcollection%p_rvectorQuickAccess1
!     rvector_u => rcollection%p_rvectorQuickAccess2

    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')

    dtstep = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)

    ! Fetching the values of rvector in the cubature pts.
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_FUNC2D, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC2D, DvaluesFevl, 2)
!     call fevl_evaluate_sim4(rvector_c, &
!                                  rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 3)
!     call fevl_evaluate_sim4(rvector_c, &
!                                  rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 4)

    ! calculate the term u_n / (1+c_n+1)^2   * c_n+1_x  (resp. c_n+1_y ) for the
    ! LHS of u_n+1 of the third chertock kurganov example

    ! first term
    Dcoefficients(1,:,:) =  -dtstep * ALPHA * ( DvaluesFevl(2,:,:) / (1+DvaluesFevl(1,:,:)* DvaluesFevl(1,:,:) ) )
    !second term
    Dcoefficients(2,:,:) =  -dtstep * ALPHA * ( DvaluesFevl(2,:,:) / (1+DvaluesFevl(1,:,:)* DvaluesFevl(1,:,:) ) )
    deallocate(DvaluesFevl)
  end subroutine

 ! ***************************************************************************


 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of an analytic given fct
  subroutine coeff_anprj_ic_nonconst (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables given by the collection
    real(DP) :: A_COEFF , B_COEFF


    A_COEFF = rcollection%DquickAccess(1)
    B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned
    Dcoefficients(1,:,:) =  A_COEFF * exp ( -B_COEFF* ( &
                ( Dpoints(1,:,:) - 0.5_DP )**2 + ( Dpoints(2,:,:) - 0.5_DP)**2)  )

  end subroutine

  ! ***************************************************************************



 ! ***************************************************************************
!<subroutine>
    ! This cb fct is used for the analytic projection of an analytic given fct
  subroutine coeff_anprj_ic_reference (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    ! local variables given by the collection
    real(DP) :: COEFF


    COEFF = rcollection%DquickAccess(1)
!     B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned
    Dcoefficients(1,:,:) =  COEFF * ( Dpoints(1,:,:) **2 - Dpoints(1,:,:) ) * ( Dpoints(2,:,:) **2 - Dpoints(2,:,:) )

  end subroutine

  ! ***************************************************************************


 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_exp (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! local variables given by the collection
    real(DP) :: COEFF_A, COEFF_B


    COEFF_A = rcollection%DquickAccess(1)
    COEFF_B = rcollection%DquickAccess(2)
!     B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned
    Dvalues(:,:) =  COEFF_A * dexp ( - COEFF_B *  ( ( Dpoints(1,:,:)  - 0.5_DP )** 2 +  ( Dpoints(2,:,:)  - 0.5_DP ) **2 ) )

  end subroutine

  ! ***************************************************************************




 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_expdistinct (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! local variables given by the collection
    real(DP) :: MASS


    MASS = rcollection%DquickAccess(1)

    ! The values of an bell-shaped exp fct. are returned
    Dvalues(:,:) =  ( MASS / ( PI * 0.01_DP ) )  *  dexp ( - ( ( Dpoints(1,:,:)  - 0.5_DP )** 2 +&
                          ( Dpoints(2,:,:)  - 0.5_DP ) **2 ) / 0.01_DP )

  end subroutine

  ! ***************************************************************************



 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the shifted exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_FCT.f90 )
    subroutine coeff_anprj_ic_exp_shifted (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! local variables given by the collection
    real(DP) :: COEFF_A, COEFF_B


    COEFF_A = rcollection%DquickAccess(1)
    COEFF_B = rcollection%DquickAccess(2)
!     B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned
    Dvalues(:,:) =  COEFF_A * dexp ( - COEFF_B *  ( ( Dpoints(1,:,:)  - 0.6_DP )** 2 +  ( Dpoints(2,:,:)  - 0.75_DP ) **2 ) )

  end subroutine

  ! ***************************************************************************


 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_2exp (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! local variables given by the collection
    real(DP) :: COEFF_A, COEFF_B


    COEFF_A = rcollection%DquickAccess(1)
    COEFF_B = rcollection%DquickAccess(2)
!     B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned
    Dvalues(:,:) =  COEFF_A * dexp ( -COEFF_B *  ( ( Dpoints(1,:,:)  - 0.35_DP )** 2 +&
                          ( Dpoints(2,:,:)  - 0.35_DP ) **2 ) ) +&
                         COEFF_A * dexp ( - COEFF_B *  ( ( Dpoints(1,:,:)  - 0.65_DP )** 2 +&
                          ( Dpoints(2,:,:)  - 0.65_DP ) **2 ) )

  end subroutine

  ! ***************************************************************************



 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_rand (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    
  intrinsic RANDOM_NUMBER
  
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    real(DP) :: random_num
     
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

    ! local variables given by the collection
    real(DP) :: COEFF_A, COEFF_B


    COEFF_A = rcollection%DquickAccess(1)
    COEFF_B = rcollection%DquickAccess(2)
!     B_COEFF = rcollection%DquickAccess(2)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! The values of an bell-shaped exp fct. are returned
    CALL RANDOM_NUMBER (random_num)
    !Dvalues(:,:) =  0.9_DP + 0.2_DP * rand(0)
    Dvalues(:,:) =  0.9_DP + 0.2_DP * random_num
    
  end subroutine

  ! ***************************************************************************




 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_hillen (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    
  intrinsic RANDOM_NUMBER
  
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    real(DP) :: random_num
    
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
    CALL RANDOM_NUMBER (random_num)
    ! Some spiky pattern
    !Dvalues(:,:) =  1.0_DP + 0.1_DP * rand(0)
    Dvalues(:,:) =  1.0_DP + 0.1_DP * random_num
  end subroutine

  ! ***************************************************************************




 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_ic_hillen2a (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! The values of an bell-shaped exp fct. at the unit square corner are returned
    Dvalues(:,:) =  1.0_DP + 0.1_DP * dexp ( -10.0_DP * ( &
                ( Dpoints(1,:,:) - 1.0_DP )**2 + ( Dpoints(2,:,:) - 1.0_DP)**2)  )

  end subroutine

  ! ***************************************************************************




 ! ***************************************************************************

!<subroutine>

    subroutine coeff_anprj_ic_one (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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

    ! returning a zero IC
    Dvalues(:,:) =  1.0_DP

  end subroutine

  ! ***************************************************************************



 ! ***************************************************************************


!<subroutine>
    ! This cb fct is used for the analytic projection of the exponential test fct
    ! E.g. this is used for every test file, which does not use a analytic given fct as
    ! a reference sol. ( like chemotaxis_cherkur_TVD_test.f90 )
    subroutine coeff_anprj_u_star (cderivative,rdiscretisation, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dvalues,rcollection)
    
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

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


    ! The values of an bell-shaped exp fct. are returned
    Dvalues(:,:) =  0!to be implemented

  end subroutine

  ! ***************************************************************************



 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_cherkur_3.f90.
    ! calculates the chemosensitivity function, since in this model  it's no more
    ! constant, but nonlinear.
    ! We'll iterate once, so there's no loop and no residual-check.
  subroutine coeff_cherkurbilf(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c, rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, CHI

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector
!     rvector_c => rcollection%p_rvectorQuickAccess1
!     rvector_u => rcollection%p_rvectorQuickAccess2
    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI


    ! Fetching the values of rvector_c in the cubature pts.
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 2)

    ! These calls are neccessary to fit the signature of f_CHI
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 4)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) =   f_CHI ( DvaluesFevl(3,icub,iel) ,&
                                                 DvaluesFevl(4,icub,iel) , CHI ) *  DvaluesFevl(1,icub,iel)
            ! second term
            Dcoefficients(2,icub,iel) =   f_CHI ( DvaluesFevl(3,icub,iel) ,&
                                                  DvaluesFevl(4,icub,iel) , CHI )* DvaluesFevl(2,icub,iel)
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************





 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_cherkur_3.f90.
    ! calculates the chemosensitivity function, since in this model  it's no more
    ! constant, but nonlinear.
    ! We'll iterate once, so there's no loop and no residual-check.
  subroutine coeff_hillenbilf(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c, rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, CHI, GAMMA, ALPHA

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector
!     rvector_c => rcollection%p_rvectorQuickAccess1
!     rvector_u => rcollection%p_rvectorQuickAccess2
    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)
    GAMMA = rcollection%DquickAccess(3)
    ALPHA = rcollection%DquickAccess(4)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI


    ! Fetching the values of rvector_c in the cubature pts.
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 2)

    ! These calls are neccessary to fit the signature of f_CHI
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 4)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) =   CHI * ( ( DvaluesFevl(3,icub,iel) * ( 1 - DvaluesFevl(3,icub,iel) / GAMMA ) ) /&
                                                    ( 1 + ALPHA * DvaluesFevl(4,icub,iel) )**2 ) *  DvaluesFevl(1,icub,iel)
            ! second term
            Dcoefficients(2,icub,iel) =   CHI * ( ( DvaluesFevl(3,icub,iel) * ( 1 - DvaluesFevl(3,icub,iel) / GAMMA ) ) /&
                                                     ( 1 + ALPHA * DvaluesFevl(4,icub,iel) )**2 ) *  DvaluesFevl(2,icub,iel)
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************





 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_hillen*
    ! this is induced by hillen's set of models
    ! The chemo sensitivity function reads
    ! \nabla ( A(u)*B(c)* \nabla c))
    ! where A and B are operators which influence the cell density distribution
    ! ( e.g. evolution of peaks ).
    ! A and B are specified as functions ( on top of this module )
  subroutine coeff_hillenX(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c, rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, CHI, GAMMA, ALPHA

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector

    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)
    GAMMA = rcollection%DquickAccess(3)
    ALPHA = rcollection%DquickAccess(4)



    ! Fetching the values of rvector_c in the cubature pts.
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 2)

    ! These calls are neccessary to fit the signature of f_CHI
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 4)


    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) =   A ( DvaluesFevl(3,icub,iel), GAMMA ) * &
                                                    B ( DvaluesFevl(4,icub,iel), CHI, ALPHA ) * &
                                                    DvaluesFevl(1,icub,iel)
            ! second term
            Dcoefficients(2,icub,iel) =   A ( DvaluesFevl(3,icub,iel), GAMMA ) * &
                                                    B ( DvaluesFevl(4,icub,iel), CHI, ALPHA )  * &
                                                    DvaluesFevl(2,icub,iel)
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************






 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_hillen2a*
    ! this is induced by hillen's model (M2a)
    ! The chemo sensitivity function reads
    ! \nabla ( u * CHI / ( 1 + ALPHA * c )**2  \cdot \nabla c)
  subroutine coeff_hillen2a(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_c, rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, CHI, GAMMA, ALPHA

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector

    rvector_c = collct_getvalue_vecsca (rcollection, "cbvector1",0,'')
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)
    GAMMA = rcollection%DquickAccess(3)
    ALPHA = rcollection%DquickAccess(4)



    ! Fetching the values of rvector_c in the cubature pts.
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_X, DvaluesFevl, 1)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_DERIV_Y, DvaluesFevl, 2)

    ! These calls are neccessary to fit the signature of f_CHI
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    call fevl_evaluate_sim4(rvector_c, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 4)


    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) =   CHI * ( DvaluesFevl(3,icub,iel)  /&
                                                    ( 1 + ALPHA * DvaluesFevl(4,icub,iel) )**2 ) *  DvaluesFevl(1,icub,iel)
            ! second term
            Dcoefficients(2,icub,iel) =   CHI * (  DvaluesFevl(3,icub,iel) /&
                                                     ( 1 + ALPHA * DvaluesFevl(4,icub,iel) )**2 ) *  DvaluesFevl(2,icub,iel)
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************




 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_cherkur_3.f90.
    ! calculates the chemosensitivity function, since in this model  it's no more
    ! constant, but nonlinear.
    ! We'll iterate once, so there's no loop and no residual-check.
  subroutine coeff_hillen_laplace(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, D_1, N

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(2,npointsPerElement,nelements))

    ! Fetching the vector
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    D_1 = rcollection%DquickAccess(2)
    N = rcollection%DquickAccess(3)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI

    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) = dtstep * D ( DvaluesFevl(1,icub,iel), D_1, N )
            ! second term
            Dcoefficients(2,icub,iel) = dtstep * D ( DvaluesFevl(1,icub,iel), D_1, N )
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************







 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_cherkur_3.f90.
    ! calculates the chemosensitivity function, since in this model  it's no more
    ! constant, but nonlinear.
    ! We'll iterate once, so there's no loop and no residual-check.
  subroutine coeff_pattern_growthterm(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep
    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(2,npointsPerElement,nelements))

    ! Fetching the vector
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI

    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 1)

    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            Dcoefficients(1,icub,iel) = dtstep *  DvaluesFevl(1,icub,iel) * ( 1.0_DP - DvaluesFevl(1,icub,iel) )
        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************




 ! ***************************************************************************
!<subroutine>
    ! Used for chemotaxis_cherkur_3.f90.
    ! calculates the chemosensitivity function, since in this model  it's no more
    ! constant, but nonlinear.
    ! We'll iterate once, so there's no loop and no residual-check.
  subroutine coeff_m8bilf(rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial, IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest
    ! The linear form which is currently to be evaluated:
    type(t_bilinearForm), intent(IN)                              :: rform
    
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    integer, dimension(:,:), intent(IN) :: IdofsTrial

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! This array contains the output of the FE evaluations, e.g. the values which
    ! we ' ll deriving the Dvalues with
    real(DP), dimension(:,:,:), allocatable :: DvaluesFevl

    ! This is the vector which is of interest
    type(t_vectorScalar) :: rvector_u

    ! Some params passed by the collection structure
    real(DP) :: dtstep, R

    integer :: icub, iel


    ! allocate some memory for the calls of Fevl
    allocate (DvaluesFevl(4,npointsPerElement,nelements))
    ! Fetching the vector
!     rvector_c => rcollection%p_rvectorQuickAccess1
!     rvector_u => rcollection%p_rvectorQuickAccess2
    rvector_u = collct_getvalue_vecsca (rcollection, "cbvector2",0,'')
    dtstep = rcollection%DquickAccess(1)
    R = rcollection%DquickAccess(2)

    ! If CHI =CHI (u,c) then we should call fevl also for certain functional descriptors
    ! to get the needed evaluations of u and c for computing CHI
    ! For the sake of simplicity, we' re now only considering a const CHI


    ! These calls are neccessary to fit the signature of f_CHI
    call fevl_evaluate_sim4(rvector_u, &
                                 rdomainIntSubset, DER_FUNC, DvaluesFevl, 3)
    ! calculate the term u_n / (1+c_{n+1})^2   * c_{n+1}_x  (resp. c_{n+1}_y ) for the
    ! LHS of u_n+1 of the third chertock kurganov example

    DO iel = 1,nelements
        DO icub = 1,npointsPerElement
            ! first term
            Dcoefficients(1,icub,iel) = f ( DvaluesFevl(3,icub,iel), R )

        END DO
    END DO

    deallocate(DvaluesFevl)

  end subroutine

 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
    ! This subroutine is used by chemotaxis_cherkur_bilf_nonlin_mdl3
    ! Here we calculate  the RHS for the c-part.
    ! This RHS is captured for fit a feasible sol of an analytically given fct.
  subroutine coeff_cherkur3_g(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) :: W, SIGMA, D_2 , BETA, dtstep

    ! for the loop
    integer :: i, j

    ! Setting the params contained in the collection
    W = rcollection%DquickAccess(1)
    SIGMA = rcollection%DquickAccess(2)
    D_2 = rcollection%DquickAccess(3)
    BETA = rcollection%DquickAccess(4)
    dtstep = rcollection%DquickAccess(5)

    ! calculate the term w*u_n^2/(SIGMA+u_n^2) for the RHS of c of the third
    ! chertock kurganov example

    DO i = 1,npointsPerElement
        DO j = 1,nelements

            Dcoefficients(1,i,j) = dtstep*( -2.0_DP * D_2 * (  Dpoints(2,i,j)*( Dpoints(2,i,j)-&
                                                    1.0_DP) + Dpoints(1,i,j) * ( Dpoints(1,i,j) - 1.0_DP ) )-&
                                                        BETA * ( W*( U_nonconst(2.0_DP,Dpoints(1,i,j),Dpoints(2,i,j))**2 ) ) / &
                                                        ( SIGMA + ( U_nonconst(2.0_DP,Dpoints(1,i,j),Dpoints(2,i,j))**2 ) ) &
                                                    )
        END DO
    END DO

  end subroutine

 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
    ! This subroutine is used by chemotaxis_cherkur_bilf_nonlin_mdl3
    ! Here we calculate  the RHS for the u-part.
    ! This RHS is captured for fit a feasible sol of an analytically given fct.
  subroutine coeff_cherkur3_f(rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
!</subroutine>

    ! local variables

    ! constants of the PDE
    real(DP) ::  D_1 , ALPHA, dtstep

    ! for the loop
    integer :: i, j
    ! Setting the params contained in the collection
    D_1 = rcollection%DquickAccess(1)
    ALPHA = rcollection%DquickAccess(2)
    dtstep = rcollection%DquickAccess(3)

    ! calculate the term w*u_n^2/(SIGMA+u_n^2) for the RHS of c of the third
    ! chertock kurganov example

   DO i = 1,npointsPerElement
        DO j = 1,nelements
            Dcoefficients(1,i,j) = dtstep*( -2.0_DP * D_1 * (  Dpoints(2,i,j)*( Dpoints(2,i,j)-&
                                                    1.0_DP) + Dpoints(1,i,j) * ( Dpoints(1,i,j) - 1.0_DP ) )+&
                                                        ALPHA * f_chemo( Dpoints(1,i,j) , Dpoints(2,i,j))&
                                                    )
        END DO
    END DO

  end subroutine

 ! ***************************************************************************


 ! ***************************************************************************
!<subroutine>
! this subroutine is used for nonconstant test fcts in the cherkur model
! E.g. u=2*c=2*x(x-1)y(y-1).
  subroutine coeff_cherkurRHSg (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    real(DP) :: dtstep
    real(DP) :: SCALE_C , SCALE_U
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    SCALE_C = rcollection%DquickAccess(2)
    SCALE_U = rcollection%DquickAccess(3)
!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!
!       else

    ! should rpresent the rhs fct g since 2*c=u=....
    ! returning the analytical computed RHS for the  chemoattractant PDE

           Dcoefficients(1,:,:) = dtstep*    (-2.0_DP * SCALE_C  * (  Dpoints(2,:,:)*( Dpoints(2,:,:)-&
                                            1.0_DP) + Dpoints(1,:,:) * ( Dpoints(1,:,:) - 1.0_DP ) ) &
                                            - SCALE_U  * ( Dpoints(2,:,:) * ( Dpoints(2,:,:)- 1.0_DP) *&
                                            Dpoints(1,:,:) * ( Dpoints(1,:,:)- 1.0_DP  ) ) + &
                                            SCALE_C  * ( Dpoints(2,:,:) * ( Dpoints(2,:,:)- 1.0_DP) *&
                                            Dpoints(1,:,:) * ( Dpoints(1,:,:)- 1.0_DP  ) ))

  end subroutine

  ! ***************************************************************************


  ! ***************************************************************************
!<subroutine>
! this subroutine is used for nonconstant test fcts in the cherkur model
! E.g. u=2*c=2*x(x-1)y(y-1).
  subroutine coeff_cherkurRHSf (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    use feevaluation
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! NOTE: The current simulation time is available via the collection
    ! rcollection in two ways:
    !  a) Parameter "TIME"
    !  b) Via the quick-access element Dquickaccess (1)
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
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT)      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  ! First of all get the vector saved in the collection to obtain the linform factor
    real(DP) :: dtstep
    real(DP) :: CHI, SCALE_C, SCALE_U
! assign the corresponding vectorentry of uvector to the needed callback coefficient

    dtstep = rcollection%DquickAccess(1)
    CHI = rcollection%DquickAccess(2)
    SCALE_C = rcollection%DquickAccess(3)
    SCALE_U = rcollection%DquickAccess(4)

!       if(present(rcollection))then
!             uvector => collct_getvalue_vecsca (rcollection, "cbvector",0,'',bexists)
!
!
!             IF(.NOT.bexists)then
! 	          call output_lbrk ()
!                   call output_line ("***************ERROR:***************")
!                   call output_line ("**********COLLECTION FAILED*********")
!                   call output_lbrk ()
!             END IF
!             call fevl_evaluate_sim4 (uvector, rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!       else
    ! this represents the fct. f for 2c=u=2*(x-1)x(y-1)y
    ! returning the analytical computed RHS for the  solution (e.g. cell density )PDE
    !-dtstep * { 4 ( y(y-1)+x(x-1) )  +  CHI [ (x²-x)² (6y²-6y+1)+(y²-y)² (6x²-6x+1 ) ] }
    Dcoefficients (1,:,:) = dtstep*( -2.0_DP * SCALE_U *( Dpoints(2,:,:)*( Dpoints(2,:,:)-1.0_DP)+Dpoints(1,:,:)*&
                                                ( Dpoints(1,:,:)-1.0_DP) )&
                                        + CHI * SCALE_C * SCALE_U * ( ( ( Dpoints(1,:,:)**2 - Dpoints(1,:,:) )**2) *&
                                                (6 * Dpoints(2,:,:)**2-6 * Dpoints(2,:,:)+1.0_DP) +&
                                                (( Dpoints(2,:,:)**2 - Dpoints(2,:,:))**2) *&
                                                (6 * Dpoints(1,:,:)**2-6 * Dpoints(1,:,:)+1.0_DP)))


!   saved
!     Dcoefficients (1,:,:) = -dtstep *  ( 4.0_DP * D_1 * ( Dpoints(2,:,:) * ( Dpoints(2,:,:)-1.0_DP) + Dpoints(1,:,:) * &
!                                     ( Dpoints(1,:,:)-1.0_DP) ) + CHI * ( ( ( Dpoints(1,:,:)**2 - Dpoints(1,:,:) )**2) *&
!  (6 * Dpoints(2,:,:)**2-6 * Dpoints(2,:,:)+1.0_DP) +  (( Dpoints(2,:,:)**2 - Dpoints(2,:,:))**2) *&
!  (6 * Dpoints(1,:,:)**2-6 * Dpoints(1,:,:)+1.0_DP)))
      

  end subroutine

  ! ***************************************************************************

end module


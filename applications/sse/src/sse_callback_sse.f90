!##############################################################################
!# ****************************************************************************
!# <name> sse_callback_sse </name>
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
!# 1.) coeff_MatrixA_SSEre   /
!#     coeff_MatrixA_SSEim   /
!#     coeff_MatrixC1_SSE2re /
!#     coeff_MatrixC1_SSE2im /
!#     coeff_MatrixC2_SSE2re /
!#     coeff_MatrixC2_SSE2im
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 2.) coeff_MatrixA_Bdr_SSEre /
!#     coeff_MatrixA_Bdr_SSEim /
!#     coeff_MatrixD1_Bdr_SSE2 /
!#     coeff_MatrixD2_Bdr_SSE2
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 3.) coeff_MatrixD_SSEre /
!#     coeff_MatrixD_SSEim
!#     -> Returns analytic valyes for the system matrix D.
!#
!# 4.) coeff_MatrixD_Bdr_SSEre /
!#     coeff_MatrixD_Bdr_SSEim
!#     -> Returns analytic valyes for the system matrix D.
!#
!# 5.) coeff_RHS_SSEre /
!#     coeff_RHS_SSEim /
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 6.) coeff_RHS_Bdr_SSEre   /
!#     coeff_RHS_Bdr_SSEim   /
!#     coeff_RHSb1_Bdr_SSEre /
!#     coeff_RHSb1_Bdr_SSEim /
!#     coeff_RHSb2_Bdr_SSEre /
!#     coeff_RHSb2_Bdr_SSEim
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 7.) getBoundaryValues_SSEre   /
!#     getBoundaryValues_SSEim   /
!#     getBoundaryValues3_SSE2   /
!#     getBoundaryValues4_SSE2   /
!#     getBoundaryValues5_SSE2   /
!#     getBoundaryValues6_SSE2
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#
!# 8.) getReferenceFunction_SSEre /
!#     getReferenceFunction_SSEim
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 9.) getReferenceDerivX_SSEre /
!#     getReferenceDerivX_SSEim
!#     -> Returns the values of the analytic derivatives in x-direction.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the recovered FE gradient in comparison to the
!#        analytic derivative
!#
!# 10.) getReferenceDerivY_SSEre /
!#      getReferenceDerivY_SSEim
!#      -> Returns the values of the analytic derivatives in x-direction.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#         $H_1$-error of the recovered FE gradient in comparison to the
!#         analytic derivative
!#
!# 11.) getReferenceDerivXX_SSEre /
!#      getReferenceDerivXX_SSEim
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#         $H_1$-error of the FE function in comparison to the analytic
!#         function
!#
!# 12.) getReferenceDerivXY_SSEre /
!#      getReferenceDerivXY_SSEim
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 13.) getReferenceDerivYX_SSEre /
!#      getReferenceDerivYX_SSEim /
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 14.) getReferenceDerivYY_SSEre /
!#      getReferenceDerivYY_SSEim
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 15.) getAnalyticValues_SSE
!#      -> Returns the values of the analytic function and its derivatives.
!#
!# 16.) getAnalyticVelocities_SSE
!#      -> Returns the values of the analytic velocities.
!#
!# </purpose>
!##############################################################################

module sse_callback_sse

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
  use sse_base_sse

  implicit none

  private

  public :: coeff_MatrixA_SSEre
  public :: coeff_MatrixA_SSEim
  public :: coeff_MatrixC1_SSE2re
  public :: coeff_MatrixC1_SSE2im
  public :: coeff_MatrixC2_SSE2re
  public :: coeff_MatrixC2_SSE2im
  public :: coeff_MatrixA_Bdr_SSEre
  public :: coeff_MatrixA_Bdr_SSEim
  public :: coeff_MatrixD1_Bdr_SSE2
  public :: coeff_MatrixD2_Bdr_SSE2
  public :: coeff_MatrixD_SSEre
  public :: coeff_MatrixD_SSEim
  public :: coeff_MatrixD_Bdr_SSEre
  public :: coeff_MatrixD_Bdr_SSEim
  public :: coeff_RHS_SSEre
  public :: coeff_RHS_SSEim
  public :: coeff_RHS_Bdr_SSEre
  public :: coeff_RHS_Bdr_SSEim
  public :: coeff_RHSb1_Bdr_SSEre
  public :: coeff_RHSb1_Bdr_SSEim
  public :: coeff_RHSb2_Bdr_SSEre
  public :: coeff_RHSb2_Bdr_SSEim
  public :: getBoundaryValues_SSEre
  public :: getBoundaryValues_SSEim
  public :: getBoundaryValues3_SSE2
  public :: getBoundaryValues4_SSE2
  public :: getBoundaryValues5_SSE2
  public :: getBoundaryValues6_SSE2
  public :: getReferenceFunction_SSEre
  public :: getReferenceFunction_SSEim
  public :: getReferenceDerivX_SSEre
  public :: getReferenceDerivX_SSEim
  public :: getReferenceDerivY_SSEre
  public :: getReferenceDerivY_SSEim
  public :: getReferenceDerivXX_SSEre
  public :: getReferenceDerivXX_SSEim
  public :: getReferenceDerivXY_SSEre
  public :: getReferenceDerivXY_SSEim
  public :: getReferenceDerivYX_SSEre
  public :: getReferenceDerivYX_SSEim
  public :: getReferenceDerivYY_SSEre
  public :: getReferenceDerivYY_SSEim
  public :: getAnalyticValues_SSE
  public :: getAnalyticVelocities_SSE

contains

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_SSEre(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) = -0.5_DP * real(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = -0.5_DP * real( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(3,ipoint,iel) = -0.5_DP * real(-cimg*(cCalpha1-cCalpha2)) ! C3
        Dcoefficients(4,ipoint,iel) = -0.5_DP * real(       cCalpha1+cCalpha2 ) ! C4
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_SSEim(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(3,ipoint,iel) = 0.5_DP * aimag(-cimg*(cCalpha1-cCalpha2)) ! C3
        Dcoefficients(4,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C4
        Dcoefficients(5,ipoint,iel) =        - dtidalfreq
      end do
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC1_SSE2re(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * real(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * real(-cimg*(cCalpha1-cCalpha2)) ! C3
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC1_SSE2im(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag(-cimg*(cCalpha1-cCalpha2)) ! C3
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC2_SSE2re(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * real( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(2,ipoint,iel) = 0.5_DP * real(       cCalpha1+cCalpha2 ) ! C4
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixC2_SSE2im(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag( cimg*(cCalpha1-cCalpha2)) ! C2
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag(       cCalpha1+cCalpha2 ) ! C4
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Bdr_SSEre(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) =&
            0.5_DP * real(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,1) ! C1
        Dcoefficients(2,ipoint,iel) =&
            0.5_DP * real( cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,2) ! C2
        Dcoefficients(3,ipoint,iel) =&
            0.5_DP * real(-cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,1) ! C3
        Dcoefficients(4,ipoint,iel) =&
            0.5_DP * real(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,2) ! C4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Bdr_SSEim(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cCalpha1
        cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                               ds*sinh(-calpha1*dh)-&
                    calpha1*dh*ds*cosh( calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cCalpha2
        cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                               ds*sinh(-calpha2*dh)-&
                    calpha2*dh*ds*cosh( calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) =&
            -0.5_DP * aimag(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,1) ! C1
        Dcoefficients(2,ipoint,iel) =&
            -0.5_DP * aimag( cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,2) ! C2
        Dcoefficients(3,ipoint,iel) =&
            -0.5_DP * aimag(-cimg*(cCalpha1-cCalpha2)) * Dnormal(ipoint,iel,1) ! C3
        Dcoefficients(4,ipoint,iel) =&
            -0.5_DP * aimag(       cCalpha1+cCalpha2 ) * Dnormal(ipoint,iel,2) ! C4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD1_Bdr_SSE2(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))
    
    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)
      
      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Store first normal component
        Dcoefficients(1,ipoint,iel) = Dnormal(ipoint,iel,1)
      end do
    end do
    
    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD2_Bdr_SSE2(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))
    
    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)
      
      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Store second normal component
        Dcoefficients(1,ipoint,iel) = Dnormal(ipoint,iel,2)
      end do
    end do
    
    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_SSEre(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
                              ds*sinh(calpha1*dz)+&
                   calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
                              ds*sinh(calpha2*dz)+&
                   calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
                
        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) = -0.5_DP * real(       cDalpha1+cDalpha2 ) ! D1
        Dcoefficients(2,ipoint,iel) = -0.5_DP * real( cimg*(cDalpha1-cDalpha2)) ! D2
        Dcoefficients(3,ipoint,iel) = -0.5_DP * real(-cimg*(cDalpha1-cDalpha2)) ! D3
        Dcoefficients(4,ipoint,iel) = -0.5_DP * real(       cDalpha1+cDalpha2 ) ! D4
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_SSEim(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
                              ds*sinh(calpha1*dz)+&
                   calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
                              ds*sinh(calpha2*dz)+&
                   calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))

        ! Compute imaginary parts of the coefficients
        Dcoefficients(1,ipoint,iel) = 0.5_DP * aimag(       cDalpha1+cDalpha2 ) ! D1
        Dcoefficients(2,ipoint,iel) = 0.5_DP * aimag( cimg*(cDalpha1-cDalpha2)) ! D2
        Dcoefficients(3,ipoint,iel) = 0.5_DP * aimag(-cimg*(cDalpha1-cDalpha2)) ! D3
        Dcoefficients(4,ipoint,iel) = 0.5_DP * aimag(       cDalpha1+cDalpha2 ) ! D4
        Dcoefficients(5,ipoint,iel) =        - dtidalfreq
      end do
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_Bdr_SSEre(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
                              ds*sinh(calpha1*dz)+&
                   calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
                              ds*sinh(calpha2*dz)+&
                   calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients
        Dcoefficients(1,ipoint,iel) =&
            0.5_DP * real(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,1) ! D1
        Dcoefficients(2,ipoint,iel) =&
            0.5_DP * real( cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,2) ! D2
        Dcoefficients(3,ipoint,iel) =&
            0.5_DP * real(-cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,1) ! D3
        Dcoefficients(4,ipoint,iel) =&
            0.5_DP * real(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,2) ! D4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixD_Bdr_SSEim(rdiscretisationTrial,rdiscretisationTest,&
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

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    complex(DP) :: cDalpha1,cDalpha2,calpha1,calpha2
    real(DP) :: dAv,dh,ds,dz
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))

    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Get global z-value from collection
    dz = rcollection%DquickAccess(1)

    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)

      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute coefficients calpha1 and calpha2
        calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
        calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)

        ! Compute coefficient cDalpha1
        cDalpha1 = dgravaccel/(dAv*(calpha1**3))*&
            ((calpha1**2)*dAv*dz*sinh(calpha1*dh)-&
                              ds*sinh(calpha1*dz)+&
                   calpha1*dz*ds*cosh(calpha1*dh))/&
            (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))

        ! Compute coefficient cDalpha2
        cDalpha2 = dgravaccel/(dAv*(calpha2**3))*&
            ((calpha2**2)*dAv*dz*sinh(calpha2*dh)-&
                              ds*sinh(calpha2*dz)+&
                   calpha2*dz*ds*cosh(calpha2*dh))/&
            (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
        
        ! Compute real parts of the coefficients multiplied by -1
        Dcoefficients(1,ipoint,iel) =&
            -0.5_DP * aimag(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,1) ! D1
        Dcoefficients(2,ipoint,iel) =&
            -0.5_DP * aimag( cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,2) ! D2
        Dcoefficients(3,ipoint,iel) =&
            -0.5_DP * aimag(-cimg*(cDalpha1-cDalpha2)) * Dnormal(ipoint,iel,1) ! D3
        Dcoefficients(4,ipoint,iel) =&
            -0.5_DP * aimag(       cDalpha1+cDalpha2 ) * Dnormal(ipoint,iel,2) ! D4
      end do
    end do

    ! Deallocate temporal memory
    deallocate(Dnormal)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_SSEre(rdiscretisation,rform, &
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_SSEim(rdiscretisation,rform, &
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_SSEre(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = dforcing

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Bdr_SSEim(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb1_Bdr_SSEre(rdiscretisation, rform,&
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

     ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))
    
    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
  
    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)
      
      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)
        
        Dcoefficients (1,ipoint,iel) = dforcing * Dnormal(ipoint,iel,1)
      end do
    end do
    
    ! Deallocate temporal memory
    deallocate(Dnormal)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb1_Bdr_SSEim(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = 0.0_DP
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb2_Bdr_SSEre(rdiscretisation, rform,&
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

     ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dnormal
    integer :: iel,ipoint

    ! Allocate temporal memory
    allocate(Dnormal(npointsPerElement,nelements,NDIM2D))
    
    ! Compute the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints,&
        Dnormal(:,:,1), Dnormal(:,:,2), 1)
    
    ! Loop over all elements
    do iel=1,size(Dcoefficients,3)
      
      ! Loop over all points per element
      do ipoint=1,size(Dcoefficients,2)

        Dcoefficients (1,ipoint,iel) = dforcing * Dnormal(ipoint,iel,2)
      end do
    end do
    
    ! Deallocate temporal memory
    deallocate(Dnormal)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHSb2_Bdr_SSEim(rdiscretisation, rform,&
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

    Dcoefficients (1,:,:) = 0.0_DP
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_SSEre(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Return Dirichlet boundary values for all situations.
    Dvalues(1) = dforcing
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_SSEim(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues3_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)
    
    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)
    
    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)
    
    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                           ds*sinh(-calpha1*dh)-&
                calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))
    
    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                           ds*sinh(-calpha2*dh)-&
                calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
    
    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4
    
    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!real(cc2*dnx+cc4*dny)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues4_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)
    
    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)
    
    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)
    
    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                           ds*sinh(-calpha1*dh)-&
                calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))
    
    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                           ds*sinh(-calpha2*dh)-&
                calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
    
    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4
    
    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!aimag(cc2*dnx+cc4*dny)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues5_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)
    
    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)
    
    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)
    
    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                           ds*sinh(-calpha1*dh)-&
                calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))
    
    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                           ds*sinh(-calpha2*dh)-&
                calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
    
    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4
    
    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!real(-cc1*dnx-cc3*dny)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues6_SSE2(Icomponents,rdiscretisation,rboundaryRegion,&
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

    ! Local variables
    complex(DP) :: cCalpha1,cCalpha2,calpha1,calpha2,cc1,cc2,cc3,cc4
    real(DP) :: dAv,dh,ds,dx,dy,dnx,dny

    ! Compute the coordinates of the point on the boundary
    call boundary_getCoords(rdiscretisation%p_rboundary, 1, dwhere, dx, dy)

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)
    
    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)
    
    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)

    ! Compute coefficients calpha1 and calpha2
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel)/dAv)
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel)/dAv)
    
    ! Compute coefficient cCalpha1
    cCalpha1 = dgravaccel/(dAv*(calpha1**3))*&
        (-(calpha1**2)*dAv*dh*sinh( calpha1*dh)-&
                           ds*sinh(-calpha1*dh)-&
                calpha1*dh*ds*cosh( calpha1*dh))/&
        (calpha1*dAv*sinh(calpha1*dh)+ds*cosh(calpha1*dh))
    
    ! Compute coefficient cCalpha2
    cCalpha2 = dgravaccel/(dAv*(calpha2**3))*&
        (-(calpha2**2)*dAv*dh*sinh( calpha2*dh)-&
                           ds*sinh(-calpha2*dh)-&
                calpha2*dh*ds*cosh( calpha2*dh))/&
        (calpha2*dAv*sinh(calpha2*dh)+ds*cosh(calpha2*dh))
    
    ! Compute the coefficients
    cc1 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C1
    cc2 = 0.5_DP * ( cimg*(cCalpha1-cCalpha2)) ! C2
    cc3 = 0.5_DP * (-cimg*(cCalpha1-cCalpha2)) ! C3
    cc4 = 0.5_DP * (       cCalpha1+cCalpha2 ) ! C4
    
    ! Compute the normal vector in the point on the boundary
    call boundary_getNormalVec2D(rdiscretisation%p_rboundary, 1, dwhere, dnx, dny)

    ! Return real part of boundary condition for x-component
    Dvalues(1) = 0.0_DP!aimag(-cc1*dnx-cc3*dny)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  ! local variables
  complex(DP) :: cC,calpha,cr1,cr2
  real(DP) :: dAv,dh,ds
  integer :: iel,ipoint
  
  select case (cderivative)
  case (DER_FUNC)
    
    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = real((cr1*exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                    cr2*exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                   (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_X)
    
    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = real(cr1*cr2*(exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                            exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                       (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_Y)
    Dvalues = 0.0_DP

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  ! local variables
  complex(DP) :: cC,calpha,cr1,cr2
  real(DP) :: dAv,dh,ds
  integer :: iel,ipoint
  
  select case (cderivative)
  case (DER_FUNC)
    
    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = aimag((cr1*exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                     cr2*exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                    (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_X)
    
    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = aimag(cr1*cr2*(exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                             exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                        (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_Y)
    Dvalues = 0.0_DP

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivX_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  ! local variables
  complex(DP) :: cC,calpha,cr1,cr2
  real(DP) :: dAv,dh,ds
  integer :: iel,ipoint

  select case (cderivative)
  case (DER_FUNC)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = real(cr1*cr2*(exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                            exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                       (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_X)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = real(cr1*cr2*(cr2*exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                            cr1*exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                           (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_Y)
    Dvalues = 0.0_DP

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivX_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  ! local variables
  complex(DP) :: cC,calpha,cr1,cr2
  real(DP) :: dAv,dh,ds
  integer :: iel,ipoint

  select case (cderivative)
  case (DER_FUNC)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = aimag(cr1*cr2*(exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                             exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                        (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_X)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = aimag(cr1*cr2*(cr2*exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                             cr1*exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                            (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_Y)
    Dvalues = 0.0_DP

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivY_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivY_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXX_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  ! local variables
  complex(DP) :: cC
  real(DP) :: dAv,calpha,dh,cr1,cr2,ds
  integer :: iel,ipoint

  select case (cderivative)
  case (DER_FUNC)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = real(cr1*cr2*(cr2*exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                            cr1*exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                           (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_X)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = real((cr1**2)*(cr2**2)*(exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                                      exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                                 (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_Y)
    Dvalues = 0.0_DP

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXX_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  ! local variables
  complex(DP) :: cC,calpha,cr1,cr2
  real(DP) :: dAv,dh,ds
  integer :: iel,ipoint

  select case (cderivative)
  case (DER_FUNC)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = aimag(cr1*cr2*(cr2*exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                             cr1*exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                            (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_X)

    do iel=1,size(Dpoints,3)
      do ipoint=1,npointsPerElement
        
        ! Compute bottom profile
        dh = sse_bottomProfile(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))

        ! Compute bottom stress
        ds = sse_bottomStress(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute vertical eddy viscosity
        dAv = sse_eddyViscosity(Dpoints(1,ipoint,iel),Dpoints(2,ipoint,iel))
        
        ! Compute coefficients
        calpha = sqrt(-cimg*dtidalfreq/dAv)
        cC     = dgravaccel/(dAv*(calpha**3))*&
            ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
        cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
        cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

        Dvalues(ipoint,iel) = aimag((cr1**2)*(cr2**2)*(exp(cr1*dlength+cr2*Dpoints(1,ipoint,iel))-&
                                                       exp(cr1*Dpoints(1,ipoint,iel)+cr2*dlength))/&
                                                  (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength)))
      end do
    end do

  case (DER_DERIV_Y)
    Dvalues = 0.0_DP

  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXY_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivXY_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYX_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYX_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYY_SSEre(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceDerivYY_SSEim(cderivative,rdiscretisation, &
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

#if defined(CASE_SSE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  elemental subroutine getAnalyticValues_SSE(dx,dy,&
      cN,cN_x,cN_y,cN_xx,cN_xy,cN_yy)

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
    complex(DP), intent(out) :: cN

    ! Values of the first derivatives
    complex(DP), intent(out) :: cN_x,cN_y

    ! Values of the second derivatives
    complex(DP), intent(out) :: cN_xx,cN_xy,cN_yy
!</output>
!</subroutine>

#if defined(CASE_SSE_ALEX)
    ! local variables
    complex(DP) :: cC,calpha,cr1,cr2
    real(DP) :: dAv,dh,ds

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)
    
    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)
    
    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)
    
    ! Compute coefficients
    calpha = sqrt(-cimg*dtidalfreq/dAv)
    cC     = dgravaccel/(dAv*(calpha**3))*&
        ((ds*sin(calpha*dh))/(calpha*dAv*sin(calpha*dh)-ds*cos(calpha*dh))+dh*calpha)
    cr1    = 0.5_DP * (1.0_DP/dlengthB + sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))
    cr2    = 0.5_DP * (1.0_DP/dlengthB - sqrt( 1.0_DP/(dlengthB**2) - 4.0_DP*cimg*dtidalfreq/cC))

    ! Solution values    
    cN = (cr1*exp(cr1*dlength+cr2*dx)-&
          cr2*exp(cr1*dx+cr2*dlength))/&
         (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength))

    ! Values of first derivatives
    cN_x = cr1*cr2*(cr2*exp(cr1*dlength+cr2*dx)-&
                    cr1*exp(cr1*dx+cr2*dlength))/&
                   (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength))
    cN_y = cmplx(0.0_DP,0.0_DP)

    ! Values of second derivatives
    cN_xx = (cr1**2)*(cr2**2)*(exp(cr1*dlength+cr2*dx)-&
                               exp(cr1*dx+cr2*dlength))/&
                          (cr1*exp(cr1*dlength)-cr2*exp(cr2*dlength))
    cN_xy = cmplx(0.0_DP,0.0_DP)
    cN_yy = cmplx(0.0_DP,0.0_DP)
#else
#error 'Test case is undefined.'
#endif

  end subroutine getAnalyticValues_SSE

  ! ***************************************************************************

!<subroutine>

  elemental subroutine getAnalyticVelocities_SSE(dx,dy,dz,cU,cV,cW)

!<description>
    ! This function computes the values of the analytic function and
    ! its derivatives.
!</description>

!<input>
    ! Coordinates
    real(DP), intent(in) :: dx,dy,dz
!</input>

!<output>
    ! Velocity values
    complex(DP), intent(out) :: cU,cV,cW
!</output>
!</subroutine>

#if defined(CASE_SSE_ALEX)
    ! local variables
    complex(DP) :: calpha,cbeta,calpha1,calpha2,cr1,cr2
    complex(DP) :: cN,cN_x,cN_y,cN_xx,cN_xy,cN_yy
    real(DP) :: dAv,dh,ds

    ! Compute bottom profile
    dh = sse_bottomProfile(dx,dy)
    
    ! Compute bottom stress
    ds = sse_bottomStress(dx,dy)
    
    ! Compute vertical eddy viscosity
    dAv = sse_eddyViscosity(dx,dy)
    
    ! Compute sea-surface elevation and its derivatives
    call getAnalyticValues_SSE(dx,dy,cN,cN_x,cN_y,cN_xx,cN_xy,cN_yy)

    ! Compute coefficients
    calpha1 = sqrt(cimg*(dtidalfreq+dcoraccel/dAv))
    calpha2 = sqrt(cimg*(dtidalfreq-dcoraccel/dAv))

    cr1 = dgravaccel/(dAv*calpha1**2)*((ds*cosh(calpha1*dz))/&
        (calpha1*dAv*sinh(calpha1*dh) + ds*cosh(calpha1*dh))-1.0_DP)*(cN_x+cimg*cN_y)
    cr2 = dgravaccel/(dAv*calpha2**2)*((ds*cosh(calpha2*dz))/&
        (calpha2*dAv*sinh(calpha2*dh) + ds*cosh(calpha2*dh))-1.0_DP)*(cN_x-cimg*cN_y)

    ! Horizontal velocities
    cU = 0.5_DP*(cr1+cr2)
    cV = 0.5_DP*(cr1-cr2)/cimg

    ! Vertical velocity
    cbeta  = sqrt(cimg*dtidalfreq/dAv);
    calpha = ds/(dAv*cbeta*sinh(cbeta*dheight)+ds*cosh(cbeta*dheight))
    cW     = dgravaccel/(cimg*dtidalfreq)*(cN_x/dlengthB - cN_xx)*&
        (dz-calpha/cbeta*sinh(cbeta*dz)) + cimg*dtidalfreq*cN

#else
#error 'Test case is undefined.'
#endif

  end subroutine getAnalyticVelocities_SSE
  
end module sse_callback_sse

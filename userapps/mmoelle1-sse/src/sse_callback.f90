!##############################################################################
!# ****************************************************************************
!# <name> sse_callback </name>
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
!# 1.) coeff_MatrixA_Poisson /
!#     coeff_MatrixA_Real /
!#     coeff_MatrixA_Aimag
!#     -> Returns analytic valyes for the system matrix A.
!#
!# 2.) coeff_MatrixD_Real /
!#     coeff_MatrixD_Aimag
!#     -> Returns analytic valyes for the system matrix D.
!#
!# 3.) coeff_RHS_Poisson /
!#     coeff_RHS_Real /
!#     coeff_RHS_Aimag /
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 4.) coeff_RHSBdr_Real /
!#     coeff_RHSBdr_Aimag
!#     -> Returns analytical values for the right hand side of the equation.
!#
!# 5.) getBoundaryValues_Poisson /
!#     getBoundaryValues_Real /
!#     getBoundaryValues_Aimag
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#
!# 6.) getReferenceFunction_Poisson /
!#     getReferenceFunction_Real /
!#     getReferenceFunction_Aimag
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 7.) getReferenceDerivX_Poisson /
!#     getReferenceDerivX_Real /
!#     getReferenceDerivX_Aimag
!#     -> Returns the values of the analytic derivatives in x-direction.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the recovered FE gradient in comparison to the
!#        analytic derivative
!#
!# 8.) getReferenceDerivY_Poisson /
!#     getReferenceDerivY_Real /
!#     getReferenceDerivY_Aimag
!#     -> Returns the values of the analytic derivatives in x-direction.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the recovered FE gradient in comparison to the
!#        analytic derivative
!#
!# 9.) getReferenceDerivXX_Poisson /
!#     getReferenceDerivXX_Real /
!#     getReferenceDerivXX_Aimag
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 10.) getReferenceDerivXY_Poisson /
!#      getReferenceDerivXY_Real /
!#      getReferenceDerivXY_Aimag
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 11.) getReferenceDerivYX_Poisson /
!#      getReferenceDerivYX_Real /
!#      getReferenceDerivYX_Aimag /
!#      -> Returns the values of the analytic function and its derivatives.
!#      -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 12.) getReferenceDerivYY_Poisson /
!#      getReferenceDerivYY_Real /
!#      getReferenceDerivYY_Aimag
!#     -> Returns the values of the analytic function and its derivatives.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# </purpose>
!##############################################################################

module sse_callback

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
  public :: coeff_MatrixA_Real
  public :: coeff_MatrixA_Aimag
  public :: coeff_MatrixD_Real
  public :: coeff_MatrixD_Aimag
  public :: coeff_RHS_Poisson
  public :: coeff_RHS_Real
  public :: coeff_RHS_Aimag
  public :: coeff_RHSBdr_Real
  public :: coeff_RHSBdr_Aimag
  public :: getBoundaryValues_Poisson
  public :: getBoundaryValues_Real
  public :: getBoundaryValues_Aimag
  public :: getReferenceFunction_Poisson
  public :: getReferenceFunction_Real
  public :: getReferenceFunction_Aimag
  public :: getReferenceDerivX_Poisson
  public :: getReferenceDerivX_Real
  public :: getReferenceDerivX_Aimag
  public :: getReferenceDerivY_Poisson
  public :: getReferenceDerivY_Real
  public :: getReferenceDerivY_Aimag
  public :: getReferenceDerivXX_Poisson
  public :: getReferenceDerivXX_Real
  public :: getReferenceDerivXX_Aimag
  public :: getReferenceDerivXY_Poisson
  public :: getReferenceDerivXY_Real
  public :: getReferenceDerivXY_Aimag
  public :: getReferenceDerivYX_Poisson
  public :: getReferenceDerivYX_Real
  public :: getReferenceDerivYX_Aimag
  public :: getReferenceDerivYY_Poisson
  public :: getReferenceDerivYY_Real
  public :: getReferenceDerivYY_Aimag

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

    Dcoefficients(:,:,:) = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_MatrixA_Real(rdiscretisationTrial,rdiscretisationTest,&
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

  subroutine coeff_MatrixA_Aimag(rdiscretisationTrial,rdiscretisationTest,&
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

  subroutine coeff_MatrixD_Real(rdiscretisationTrial,rdiscretisationTest,&
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

  subroutine coeff_MatrixD_Aimag(rdiscretisationTrial,rdiscretisationTest,&
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

  subroutine coeff_RHS_Poisson(rdiscretisation,rform, &
      nelements,npointsPerElement,Dpoints, &
      IdofsTest,rdomainIntSubset,&
      Dcoefficients,rcollection)
    
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

  subroutine coeff_RHS_Real(rdiscretisation,rform, &
      nelements,npointsPerElement,Dpoints, &
      IdofsTest,rdomainIntSubset,&
      Dcoefficients,rcollection)
    
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

  subroutine coeff_RHS_Aimag(rdiscretisation,rform, &
      nelements,npointsPerElement,Dpoints, &
      IdofsTest,rdomainIntSubset,&
      Dcoefficients,rcollection)
    
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

  subroutine coeff_RHSBdr_Real(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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

  subroutine coeff_RHSBdr_Aimag(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_Real(Icomponents,rdiscretisation,rboundaryRegion,&
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

  subroutine getBoundaryValues_Aimag(Icomponents,rdiscretisation,rboundaryRegion,&
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

  subroutine getReferenceFunction_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceFunction_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivX_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivX_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivY_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivY_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivXX_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivXX_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivXY_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivXY_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivYX_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivYX_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivYY_Real(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
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

  subroutine getReferenceDerivYY_Aimag(cderivative,rdiscretisation, &
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

#if defined(CASE_ALEX)
  select case (cderivative)
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
#else
#error 'Test case is undefined.'
#endif
  
  end subroutine

end module sse_callback

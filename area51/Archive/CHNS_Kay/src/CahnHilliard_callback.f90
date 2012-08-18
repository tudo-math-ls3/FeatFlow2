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
!# 2.) coeff_RHS_phi
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 2') coeff_RHS_w
!#
!# 3.) getBoundaryValues
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 4.) CH_iniconPhi
!# 5.) CH_iniconChemP
!#
!# 6.) coeff_NonlinearMass
!# 7.) coeff_VarMass
!# 8.) coeff_VarLaplace
!#
!#
!# </purpose>
!##############################################################################

module CahnHilliard_callback

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use mprimitives
  use derivatives

  use feevaluation

  use CahnHilliard_basic
  use ccbasic

  IMPLICIT NONE

CONTAINS

! ***************************************************************************
  !<subroutine>
  subroutine CH_iniconPhi(x, y, f_val)
  
!<description>
  ! Give you the initial soluiont of Cahn-Hilliard equation, for c
!</description>

!<inputoutput>
  real(DP), intent(in) :: x, y
  real(DP), intent(out) :: f_val
!</inputoutput>

!</subroutine>

 ! local variables
  real(DP) :: delta

!    delta=(x-0.5_DP)**2+(y-0.5_DP)**2+1.0_DP
!    f_val=delta
 
    delta=(x-0.5_DP)**2/0.01_DP+(y-0.5_DP)**2/0.0225_DP-1.0_DP
	f_val=tanh(delta)

!     f_val=1.0_DP

!    if ((x .ge. 0.35) .and. (x .le. 0.65) .and. (y .ge. 0.35) .and. (y .le. 0.65)) then
!       f_val=1.0_DP
!    else
!       f_val=-1.0_DP
!    end if

!    delta=(x-0.5_DP)**2/0.01_DP+(y-0.5_DP)**2/0.0225_DP-1.0_DP
!    f_val=tanh(delta)

! if ((((x - 0.4)**2  + (y - 0.5)**2)  .le.  0.01_DP) .or.  &
!           (((x - 0.6)**2  + (y - 0.5)**2) .le.  0.01_DP)) then
!        f_val = 1.0_DP
!      else
!        f_val = -1.0_DP
!       end if
	
  end subroutine

! ***************************************************************************
  !<subroutine>
  subroutine CH_iniconChemP(x, y, f_val)
  
!<description>
  ! Give you the initial soluiont of Cahn-Hilliard equation, for w
!</description>

!<inputoutput>
  real(DP), intent(in) :: x, y
  real(DP), intent(out) :: f_val
!</inputoutput>

!</subroutine>

 ! local variables
  real(DP) :: delta, phi
  real(DP) :: eps=0.02_DP

!    delta=(x-0.5_DP)**2+(y-0.5_DP)**2+1.0_DP
!    f_val=delta
 
 ! The chemical potential is
  ! w = (phi^3 - phi)/eps - eps * Laplace \phi

!  f_val=(phi**3-phi)/ep-ep*(-80000*phi*(1-phi**2)*(x-0.5_DP)**2+&
!          2600/9-2600/9*phi**2-1280000/81*phi*(1-phi**2)*(y-0.5_DP)**2)
  ! The following one is based on matlab, the above one is based on Chain rule.
!Feng's setting
!  phi= tanh(100*(x-0.5_DP)**2+400.0_DP/9.0_DP*(y-0.5_DP)**2-1.0_DP)
!  f_val=2.0_DP*phi*(1.0_DP-phi**2)*(200.0_DP*x-100.0_DP)**2-2600.0_DP/9.0_DP+&
!    2600.0_DP/9.0_DP*phi**2+2.0_DP*phi*(1.0_DP-phi**2)*(800.0_DP/9.0_DP*y-400.0_DP/9.0_DP)**2+&
!	(phi**3-phi)/eps**2

 ! Kay's setting
   phi=tanh((x-0.5_DP)**2/0.01_DP+(y-0.5_DP)**2/0.0225_DP-1.0_DP)
   f_val=(phi**3-phi)/eps-eps*(-2.0_DP*phi*(1.0_DP-phi**2)*(200.0_DP*x-100.0_DP)**2+&
         2600.0_DP/9.0_DP-2600.0_DP/9.0_DP*phi**2-&
         2.0_DP*phi*(1-phi**2)*(800.0_DP/9.0_DP*y-400.0_DP/9.0_DP)**2)
 
 

!     f_val=1.0_DP

!    if ((x .ge. 0.35) .and. (x .le. 0.65) .and. (y .ge. 0.35) .and. (y .le. 0.65)) then
!       f_val=1.0_DP
!    else
!       f_val=-1.0_DP
!    end if

!    delta=(x-0.5_DP)**2/0.01_DP+(y-0.5_DP)**2/0.025_DP-1.0_DP
!    f_val=tanh(delta)

! if ((((x - 0.4)**2  + (y - 0.5)**2)  .le.  0.01_DP) .or.  &
!           (((x - 0.6)**2  + (y - 0.5)**2) .le.  0.01_DP)) then
!        f_val = 1.0_DP
!      else
!        f_val = -1.0_DP
!       end if
	
  end subroutine


!***************************************************************************
! MCai
!<subroutine>
  subroutine coeff_nonlinear (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration


  !<description>
    ! This subroutine is to calculate nonlinear term in assembling defect:
	! Specifically, it corresponds to (C^3-C, \psi), \psi is test function.
	! We use linearform to do this L^2 inner product. See also ACNS applications.
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
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), DIMENSION(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), OPTIONAL      :: rcollection

  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.

    real(DP), DIMENSION(:,:,:), intent(OUT)          :: Dcoefficients
  
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints,nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues

    integer :: i, j

	! rInnervector is to save rACvector from previous step
    type(t_vectorblock),pointer :: rInnervector
    real(DP) :: f_val

    !  rOutervector is to save rNSvector
!    type(t_vectorblock) :: rOutvector
  !</output>

!</subroutine>


! we treat the source term explicitly (f(\phi(t_n)), \psi)
  if (present(rcollection)) then

!Mcai, we need to modify something here   ! Refer to nonlinear_calcDefect
      rInnervector=>rcollection%p_rvectorQuickAccess2

      allocate(Dvalues(1,npointsPerElement,nelements))

! We need rInnervector
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues, 1)

! Dcoefficients only have 1 itemcount
      do i = 1,nelements
         do j = 1, npointsPerElement
! Debug, check Dvalues are correct
!              print *, Dvalues(1, j, i)
	       call cubicfun_val(Dvalues(1,j,i), f_val)
! Dcoefficients only have 1 item
           Dcoefficients(1,j,i)=f_val
         end do
      end do

      deallocate(Dvalues)
  end if

  end subroutine

  subroutine cubicfun_val(phi, f_val)
  real(DP), intent(IN) :: phi
  real(DP), intent(INOUT) :: f_val

     f_val=(phi**2-1.0_DP)*phi
  end subroutine

! ***************************************************************************

  !<subroutine>
  subroutine coeff_NonlinearPrec (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly: CH_assemblyMatrix.
	! It is used to construct preconditioner for nonlinear system, see
	! CH_generateNonlinearMat
    !
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
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
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
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rInnervector
    real(DP) :: f_val
	real(DP) :: eps=0.02_DP
  !</output>

  !</subroutine>

!  Dcoefficients = 0.0_DP

! we treat the source term explicitly (f(\phi(t_n)), \psi)
    if (present(rcollection)) then

!Mcai, we need to modify something here
      rInnervector=>rcollection%p_rvectorQuickAccess1
      allocate(Dvalues(1,npointsPerElement,nelements))

! We need function values at quad points
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues, 1)

! Dcoefficients only have 1 itemcount
      do i = 1,nelements
        do j = 1, npointsPerElement
! Debug, check Dvalues are correct
!           print *, Dvalues(1, j, i)
	      call Cubicderivative_fun(Dvalues(1,j,i), f_val)
! Dcoefficients only have 1 item: our formulation is based on Kay's paper
          Dcoefficients(1,j,i)=f_val/eps
        end do
      end do

      deallocate(Dvalues)
    end if

  end subroutine

  subroutine Cubicderivative_fun(c, f_val)
  real(DP), intent(IN) :: c
  ! f_val is the mobility
  real(DP), intent(INOUT) :: f_val
  ! local variable
  ! We should use derivative, Newton like nonlinear solver: Jacobian.
     f_val=3.0_DP*c**2-1.0_DP
    
  end subroutine

!**********************************************************************************
 subroutine coeff_Conv (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    ! In this example, we compute the poisson example with a nonconstant
    ! coefficient depending on a finite element function. The FE function is
    ! passed to this routine via the collection structure rcollection.
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
    integer, intent(IN)                        :: nelements
    
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
    real(DP), dimension(:,:,:), intent(OUT)          :: Dcoefficients
  !</output>
  !</subroutine>
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues_X, Dvalues_Y

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rOutervector

  !</output>

  !</subroutine>

!  Dcoefficients = 0.0_DP

    if (present(rcollection)) then
!Mcai, p_rvectorQuickAccess2 is for rNSvector
      rOutervector=>rcollection%p_rvectorQuickAccess1

      allocate(Dvalues_X(1,npointsPerElement,nelements))
      allocate(Dvalues_Y(1,npointsPerElement,nelements))

! We need function values at quad points: RvectorBlock(1) for 1st component of vel
      call fevl_evaluate_sim4 (rOutervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues_X, 1)
      call fevl_evaluate_sim4 (rOutervector%RvectorBlock(2), &
	                         rdomainIntSubset, DER_FUNC, Dvalues_Y, 1)


! Dcoefficients only have 1 itemcount
      do i = 1,nelements
        do j = 1, npointsPerElement
          ! MCai: {\bf u} \cdot \grad \phi
          ! Dcoefficients have 2 items?
          Dcoefficients(1,j,i)=Dvalues_X(1,j,i)
    	  Dcoefficients(2,j,i)=Dvalues_Y(1,j,i)
        end do
      end do

      deallocate(Dvalues_X)
      deallocate(Dvalues_Y)
    end if
  
  end subroutine

! ***************************************************************
   subroutine Coeff_conv1 (rdiscretisationTest,rdiscretisationTrial,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset,&
                  Dcoefficients, rcollection)
    
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
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The trilinear form which is currently being evaluated:
    type(t_trilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer(PREC_ELEMENTIDX), intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional ,target     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>

  !</subroutine>
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues_X

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rOutervector

  !</output>

  !</subroutine>

!  Dcoefficients = 0.0_DP

    if (present(rcollection)) then
!Mcai, p_rvectorQuickAccess2 is for rNSvector
      rOutervector=>rcollection%p_rvectorQuickAccess2

      allocate(Dvalues_X(1,npointsPerElement,nelements))

! We need function values at quad points: RvectorBlock(1) for 1st component of vel
      call fevl_evaluate_sim4 (rOutervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues_X, 1)

! Dcoefficients only have 1 itemcount
      do i = 1,nelements
        do j = 1, npointsPerElement
          ! MCai: {\bf u} \cdot \grad \phi
          ! Dcoefficients have 2 items?
          Dcoefficients(1,j,i)=Dvalues_X(1,j,i)
        end do
      end do

      deallocate(Dvalues_X)
    end if

  end subroutine

! ***************************************************************
   subroutine Coeff_conv2 (rdiscretisationTest,rdiscretisationTrial,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset,&
                  Dcoefficients, rcollection)
    
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
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The trilinear form which is currently being evaluated:
    type(t_trilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
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
    type(t_collection), intent(INOUT), optional ,target     :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>

  !</subroutine>
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues_Y

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rOutervector

  !</output>

  !</subroutine>

!  Dcoefficients = 0.0_DP

    if (present(rcollection)) then

!MCai, p_rvectorQuickAccess2 is for rNSvector
      rOutervector=>rcollection%p_rvectorQuickAccess2
      allocate(Dvalues_Y(1,npointsPerElement,nelements))

! We need function values at quad points:RvectorBlock(2) for 2nd component of vel
      call fevl_evaluate_sim4 (rOutervector%RvectorBlock(2), &
	                         rdomainIntSubset, DER_FUNC, Dvalues_Y, 1)

! Dcoefficients only have 1 itemcount
      do i = 1,nelements
        do j = 1, npointsPerElement
          ! MCai: {\bf u} \cdot \grad \phi
          ! Dcoefficients have 2 items?
          Dcoefficients(1,j,i)=Dvalues_Y(1,j,i)
        end do
      end do
		
      deallocate(Dvalues_Y)
    end if

  end subroutine


! ***************************************************************************

  !<subroutine>

  subroutine coeff_VarMass (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It gives the coeff
	! for variable dependent mass matrix: the coefficient depends on c.
    !
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
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
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
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rInnervector
    real(DP) :: f_val
  !</output>

  !</subroutine>


! we treat the source term explicitly (f(\phi(t_n)), \psi)
    if (present(rcollection)) then

!Mcai, we need to modify something here
      rInnervector=>rcollection%p_rvectorQuickAccess1
      allocate(Dvalues(1,npointsPerElement,nelements))

! We need function values at quad points
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues, 1)

! Dcoefficients only have 1 itemcount
      do i = 1,nelements
        do j = 1, npointsPerElement
! Debug, check Dvalues are correct
!           print *, Dvalues(1, j, i)
	      call Density_fun(Dvalues(1,j,i), f_val)
! Dcoefficients only have 1 item
          Dcoefficients(1,j,i)=f_val
        end do
      end do

      deallocate(Dvalues)
        
    end if

  end subroutine

  subroutine Density_fun(c, f_val)
  real(DP), intent(IN) :: c
  ! f_val is the mobility
  real(DP), intent(INOUT) :: f_val
  ! local variable
	! Two constant densities
	real(DP) :: rho_1=1.0_DP
	real(DP) :: rho_2=1.0_DP

    f_val=rho_1*(1.0_DP+c)/2.0_DP+rho_2*(1.0_DP-c)/2.0_DP
    
  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_VarLaplace (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It is used to assemble
	! variable coefficient Laplace matrix, the coefficient depends on c.
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
    INTEGER, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest
    
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
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rInnervector
    real(DP) :: f_val
  !</output>

  !</subroutine>

!  Dcoefficients = 0.0_DP

! we treat the source term explicitly (f(\phi(t_n)), \psi)
    if (present(rcollection)) then

!Mcai, we need to modify something here
      rInnervector=>rcollection%p_rvectorQuickAccess1
      allocate(Dvalues(1,npointsPerElement,nelements))

! We need function values at quad points
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues, 1)

! Dcoefficients only have 2 itemcount
      do i = 1,nelements
        do j = 1, npointsPerElement
! Debug, check Dvalues are correct
!           print *, Dvalues(1, j, i)
	      call Mobility_fun(Dvalues(1,j,i), f_val)
! Dcoefficients have 2 items: because it corresponds to Laplace
          Dcoefficients(1,j,i)=f_val
		  Dcoefficients(2,j,i)=f_val
        end do
      end do
      deallocate(Dvalues)
    end if

  end subroutine

! MCai, we need to rewrite Laplace term for CH model, because it is \grad (w/\rho(c))
! we need to use chain role to calculate it.

  subroutine Mobility_fun(c, f_val)
  real(DP), intent(IN) :: c
  ! f_val is the mobility
  real(DP), intent(INOUT) :: f_val
  ! local variable
    real(DP) :: b, rho
	! Two constant densities
	real(DP) :: rho_1=1.0_DP
	real(DP) :: rho_2=1.0_DP

	b=max(1.0_DP-c**2, 0.0_DP)
	b=b**2
 
    rho=1.0_DP + (c-1.0_DP)*(rho_1-rho_2)/(2.0_DP*rho_1)
    f_val=b/(rho**2)

  end subroutine

!************************************************************************************
  subroutine CH_initCollectForAssembly (rCHproblem,rcollection)
  
!<description>
  ! This subroutine is an auxiliary subroutine called by the CC2D framework
  ! and has usually not to be changed by the user.
  !
  ! The subroutine prepares the collection rcollection to be passed to callback
  ! routines for assembling boundary conditions or RHS vectors. It is
  ! called directly prior to the assembly to store problem-specific information
  ! in the quick-access arrays of the collection.
  ! Basically speaking, the routine stores information about whether thw problem
  ! is stationary, nonstationary or about the current simulation time.
!</description>

!<input>
  ! Problem structure with all problem relevant data.
  type(t_CHproblem), intent(IN) :: rCHproblem
!</input>

!<inputoutput>
  ! Collection structure to be initialised.
  type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
  
!</subroutine>

    ! In a nonstationary simulation, save the simulation time as well as the
    ! minimum and maximum time to the quick-access array of the collection,
    ! so it can be accessed in the callback routines!
    rcollection%Iquickaccess(1) = rCHproblem%itimedependence
    select case (rCHproblem%itimedependence)
    case (0)
      ! Stationary simulation
      rcollection%Dquickaccess(1) = 0.0_DP
      rcollection%Dquickaccess(2) = 0.0_DP
      rcollection%Dquickaccess(3) = 0.0_DP
    case (1)
      rcollection%Dquickaccess(1) = rCHproblem%rtimedependence%dtime
      rcollection%Dquickaccess(2) = rCHproblem%rtimedependence%dtimemin
      rcollection%Dquickaccess(3) = rCHproblem%rtimedependence%dtimemax
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_phi (rdiscretisation,rform, &
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
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), DIMENSION(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), DIMENSION(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients(:,:,:) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_w (rdiscretisation,rform, &
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
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), DIMENSION(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), DIMENSION(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 0.0_DP

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
  integer, DIMENSION(:), intent(IN)                           :: Icomponents

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
!  type(t_collection), intent(IN), OPTIONAL      :: rcollection
!=======
  type(t_collection), intent(INOUT), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), DIMENSION(:), intent(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

  ! To get the X/Y-coordinates of the boundary point, use:
  !
  ! real(DP) :: dx,dy
  !
  ! call boundary_getCoords(rdiscretisation%p_rboundary, &
  !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

  ! Return zero Dirichlet boundary values for all situations.
  Dvalues(1) = 0.0_DP
  
  IF ((dwhere .GE. 0.0_DP) .AND. (dwhere .LE. 1.0_DP)) &
    Dvalues(1) = mprim_getParabolicProfile (dwhere,1.0_DP,1.0_DP)

  end subroutine

 ! ***************************************************************************
  !<subroutine>
  ! MCai, we do not need the following subroutine

  subroutine coeff_CahnHilliard (rdiscretisationTrial,rdiscretisationTest,rform, &
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
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), DIMENSION(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, DIMENSION(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), DIMENSION(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 1.0_DP

  end subroutine


end module

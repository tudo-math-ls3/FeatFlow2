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
!# 0.) AllenCahn_inicon
!#     -> Returns the solution of initial solution
!#
!# 1.) coeff_Laplace
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2^*) coeff_RHS0
!# 2^#) coeff_RHS1
!#  2.) coeff_RHS
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
!# 4.) coeff_Conv1
!# 5.) coeff_Conv2
!# 6.) coeff_Poly
!#
!# </purpose>
!##############################################################################

MODULE AllenCahn_callback

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

  use feevaluation
! find the following module in this subdirectory.
  use fe_cub_evaluation

  IMPLICIT NONE

!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  real(DP), parameter :: gamma     =  0.01_DP

  ! Preconditioning with inverse mass matrix (not yet implemented)
  real(DP), parameter :: eps   = 0.01_DP

CONTAINS

  
  ! ***************************************************************************

!<subroutine>

  subroutine AllenCahn_inicon(x, y, f_val)
  
!<description>
  ! Give you the initial soluiont of Allen-Cahn equation
!</description>

!<inputoutput>
  REAL(DP), intent(in) :: x, y
  REAL(DP), intent(out) :: f_val
!</inputoutput>

!</subroutine>

 ! local variables
  REAL(DP) :: delta

! A circle, for debugging.
!    delta=(x-0.5_DP)**2+(y-0.5_DP)**2+1.0_DP
!    f_val=delta


! FengHeLiu's MathComp paper, elliptic bubble will become to a circular one
!    delta=(x-0.5_DP)**2/0.01_DP+(y-0.5_DP)**2/0.0225_DP-1.0_DP
!    f_val=tanh(delta)

   delta=(x-0.5_DP)**2/0.09_DP+(y-0.5_DP)**2/0.09_DP-1.0_DP
    f_val=tanh(delta/eps)

! Feng's example
!     delta=x**2/0.01_DP+y**2/0.0225_DP-1.0_DP
!     f_val=tanh(delta)
!

   ! Feng Shen
!    delta=sqrt((x-0.5_DP)**2+(y-0.5_DP)**2)-50.0_DP/128.0_DP
!    f_val=-tanh(delta/eps)



!   f_val=1.0_DP

! A circle, which will vanish, similar to Feng Shen's JCP paper
!    delta=(x-0.5_DP)**2/0.01_DP+(y-0.5_DP)**2/0.0225_DP-1.0_DP
!    f_val=tanh(delta)

! delta=sqrt(x**2+y**2) -100.0_DP/128.0_DP
! f_val=-tanh(delta/eps)

! constant, for debuggin
!     f_val=1.0_DP

! Square bubble
!    if ((x .ge. 0.35) .and. (x .le. 0.65) .and. (y .ge. 0.35) .and. (y .le. 0.65)) then
!       f_val=1.0_DP
!    else
!       f_val=-1.0_DP
!    end if


! Kisssing bubble
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
	! Specifically, it corresponds to (C**3-C, \psi), \psi is test function.
	! We use linearform to do this L**2 inner product. See also ACNS applications.
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
      rInnervector=>rcollection%p_rvectorQuickAccess1

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
           Dcoefficients(1,j,i)=gamma/(eps**2)*f_val
         end do
      end do

      deallocate(Dvalues)
  end if

  end subroutine

  subroutine cubicfun_val(phi, f_val)
  real(DP), intent(IN) :: phi
  real(DP), intent(INOUT) :: f_val

     f_val=(phi**2- 1.0_DP)*phi
  end subroutine
! ***************************************************************************
  !<subroutine>
  subroutine coeff_Jacobian(rdiscretisationTrial,rdiscretisationTest,rform, &
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
! Note that we do not need gamma and eps here, because rnonlinearIteration%dgamma include
! gamma/eps**2

  !</subroutine>
! we treat the nonlinear term as:
!    (\gamma/eps^2)((\phi(t_n)^2-1)\phi(t_n+1), \psi)

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
! Dcoefficients only have 1 item
          Dcoefficients(1,j,i)=gamma/(eps**2)*(3.0_DP*Dvalues(1,j,i)**2-1.0_DP)
        end do
        
      end do
      deallocate(Dvalues)
    end if
  end subroutine


! ***************************************************************************
  !<subroutine>
  subroutine coeff_Poly(rdiscretisationTrial,rdiscretisationTest,rform, &
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
  !</output>

  !</subroutine>
! we treat the nonlinear term as:
!    (\gamma/eps^2)((\phi(t_n)^2-1)\phi(t_n+1), \psi)

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
! Dcoefficients only have 1 item
          Dcoefficients(1,j,i)=gamma/(eps**2)*(Dvalues(1,j,i)**2-1.0_DP)
        end do
        
      end do
      deallocate(Dvalues)
    end if
  end subroutine

!**********************************************************************************
 subroutine coeff_Conv1 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
  !</subroutine>
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues_X

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rOutervector

  !</output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)          :: Dcoefficients

  !</subroutine>

!  Dcoefficients = 0.0_DP

    if (present(rcollection)) then
!Mcai, p_rvectorQuickAccess2 is for rNSvector
      rOutervector=>rcollection%p_rvectorQuickAccess1

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

!**********************************************************************************
 subroutine coeff_Conv2 (rdiscretisationTrial,rdiscretisationTest,rform, &
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
  !</subroutine>
    ! local variables
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints, nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues_Y

    integer :: i, j

	! rInnervector is to save rvector from previous step
    type(t_vectorblock),pointer :: rOutervector

  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)          :: Dcoefficients

  !</subroutine>

!  Dcoefficients = 0.0_DP

    if (present(rcollection)) then

!MCai, p_rvectorQuickAccess2 is for rNSvector
      rOutervector=>rcollection%p_rvectorQuickAccess1
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

  subroutine coeff_AnalyticSolution_AC (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the AC problem \phi in the
    ! initial solution vector.
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
    ! DIMENSION(\#local DOF's in test space,nelements)
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
  
    ! REAL(DP) :: dtime
    ! REAL(DP) :: dx,dy
    !
    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    !
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    !
    ! dtime = 0.0_DP
    ! IF (PRESENT(rcollection)) dtime = rcollection%Dquickaccess(1)
    !
    ! -----
    ! In the basic implementation, we just call ffunction_TargetX to get the
    ! values -- so the target function (which is normally used for
    ! calculating the error to a reference function) is also the definition
    ! for the initial solution.
    ! If necessary, this behaviour can be changed here.
    
    call AC5ffunction_TargetX (DER_FUNC,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dcoefficients(1,:,:),rcollection)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine AC5ffunction_TargetX (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine doesn't make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
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

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
  integer :: i, j
!</output>

!</subroutine>

    do i = 1,nelements
       do j = 1, npointsPerElement
         call AllenCahn_inicon(Dpoints(1,j, i), Dpoints(2, j, i), Dvalues(j, i))
       end do
    end do

     ! Example:
    ! IF (cderivative .EQ. DER_FUNC) THEN
    !   Dvalues(:,:) = (-dtime**2/100.+dtime/5.)*(Dpoints(1,:,:))
    ! END IF
  end subroutine

  ! ***************************************************************************

  ! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_AllenCahn (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
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
    ! analytic boundary boundary description etc.; trial space.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisationTest
    
    ! The bilinear form which is currently being evaluated:
    TYPE(t_bilinearForm), INTENT(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                  :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 1.0_DP

  END SUBROUTINE

  ! ***************************************************************************

! MCai
!***************************************************************************
! calculate the coeff_RHS from AC problem

!<subroutine>
  SUBROUTINE coeff_RHS0 (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration

    use ccbasic

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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.

    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)        :: Dcoefficients
  
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints,nelements).
    real(DP), dimension(:,:,:), allocatable :: Dvalues

    integer :: i, j

	! rInnervector is to save rACvector from previous step
    type(t_vectorblock),pointer :: rInnervector


    ! rOutervector is to save rNSvector
!    type(t_vectorblock) :: rOutvector
  !</output>

  !</subroutine>

  Dcoefficients = 0.0_DP


  END SUBROUTINE

! MCai
!***************************************************************************
! calculate the coeff_RHS from AC problem

!<subroutine>
  SUBROUTINE coeff_RHS1 (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration

    use ccbasic

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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.

    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  
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

!  Dcoefficients = 0.0_DP

! we treat the source term explicitly (f(\phi(t_n)), \psi)
  if (present(rcollection)) then

!Mcai, we need to modify something here
      rInnervector=>rcollection%p_rvectorQuickAccess1
      allocate(Dvalues(1,npointsPerElement,nelements))

!      call AC5_initSolution(rACproblem, rACvector)
      
!
! We need rACvector
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues, 1)

! Dcoefficients only have 1 itemcount
      do i = 1,nelements
         do j = 1, npointsPerElement
! Debug, check Dvalues are correct
!        	  print *, Dvalues(1, j, i)
    	      call fun_val(Dvalues(1,j,i), f_val)
! Dcoefficients only have 1 item
              Dcoefficients(1,j,i)=-gamma/(eps**2)*f_val
         end do
      end do

      deallocate(Dvalues)
!MCai, we shall also release the pointers: rInnervector
        
  end if
 !     call lsysbl_releasevector(rinnervector)
  end subroutine

  subroutine fun_val(phi, f_val)
  real(DP), intent(IN) :: phi
  real(DP), intent(INOUT) :: f_val

     f_val=(phi**2-1.0_DP)*phi

  end subroutine

! MCai
!***************************************************************************
! calculate the coeff_RHS from NS problem
  SUBROUTINE coeff_RHS (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration

    use ccbasic

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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER, DIMENSION(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection

  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.

    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints,nelements).
	! for \nabla phi(t_n)
    real(DP), dimension(:,:,:), allocatable :: Dvalues_derX
    real(DP), dimension(:,:,:), allocatable :: Dvalues_derY
	! for two components of velocity
    real(DP), dimension(:,:,:), allocatable :: Dvalues_X
    real(DP), dimension(:,:,:), allocatable :: Dvalues_Y
    integer :: i, j

	! rInnervector is to save rACvector from previous step
    ! type(t_vectorblock) :: rInnervector
    ! rOutvector is to save rNSvector
    type(t_vectorblock),pointer :: rInnervector
    type(t_vectorblock),pointer :: rOutvector
!    real(DP) :: dtime

!***************rdomainIntSubset_NS****************************
!    TYPE(t_domainIntSubset)  :: rdomainIntSubset_NS
!**************************************************************
  
!</output>

  !</subroutine>

  if (present(rcollection)) then

!MCai, we need to modify something here
      rInnervector=>rcollection%p_rvectorQuickAccess1
      rOutvector => rcollection%p_rvectorQuickAccess2

      allocate(Dvalues_derX(1,npointsPerElement,nelements))
      allocate(Dvalues_derY(1,npointsPerElement,nelements))

      allocate(Dvalues_X(1,npointsPerElement,nelements))
      allocate(Dvalues_Y(1,npointsPerElement,nelements))

! We need dphi/dx and dphi/dy
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_DERIV2D_X, Dvalues_derX, 1)
! MCai, debug.
! Check Dvalues_derX is a small number....
      call fevl_evaluate_sim4 (rInnervector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_DERIV2D_Y, Dvalues_derY, 1)

!****Here, we need use conforming discretisation to recover "rdomainIntSubset"
! for nonconforming discretisation********************************************
! We may need x-velocity and y-velocity
      call fevl_evaluate_sim4 (rOutvector%RvectorBlock(1), &
	                         rdomainIntSubset, DER_FUNC, Dvalues_X, 1)
      call fevl_evaluate_sim4 (rOutvector%RvectorBlock(2), &
	                         rdomainIntSubset, DER_FUNC, Dvalues_Y, 1)
!****We should use the following code, rather than the above*****************
! How to get rdomainIntSubset_NS????????? An example can be found in
!  Postprocessing/pprocgradients/
!      call fevl_evaluate_sim4 (rOutvector%RvectorBlock(1), &
!	                         rdomainIntSubset_NS, DER_FUNC, Dvalues_X, 1)
!      call fevl_evaluate_sim4 (rOutvector%RvectorBlock(2), &
!	                         rdomainIntSubset_NS, DER_FUNC, Dvalues_Y, 1)
!****************************************************************************

! Dcoefficients have 1 itermcount: given by {\bf u} \cdot \nabla \phi(t_n)
      do i = 1,nelements
          do j = 1, npointsPerElement
! MCai
!          print *, Dvalues_derX(1,j,i)
!          print *, Dvalues_derY(1,j,i)
!          print *, Dvalues_X(1,j,i)
!          print *, Dvalues_Y(1,j,i)

          Dcoefficients (1,j,i)=-Dvalues_derX(1,j,i)*Dvalues_X(1, j, i)&
                                       -Dvalues_derY(1,j,i)*Dvalues_Y(1, j, i)
          end do
       end do
       
       deallocate(Dvalues_derX)
       deallocate(Dvalues_derY)
       deallocate(Dvalues_X)
       deallocate(Dvalues_Y)
!MCai, we shall also release the pointers: rInnervector and rOutvector

  end if
! call lsysbl_releasevector(routvector)

  END SUBROUTINE

! MCai
!***************************************************************************
! calculate the coeff_RHS3, this code is for mass conservation
!
  subroutine coeff_RHS3 (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration

    use ccbasic

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
    INTEGER, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, intent(IN)                                         :: npointsPerElement
    
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
    type(t_collection), intent(INOUT), OPTIONAL, target  :: rcollection

  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.

    real(DP), DIMENSION(:,:,:), intent(OUT)          :: Dcoefficients
  
    ! Values of the FE function at the points specified by Dpoints.
    ! DIMENSION(itermCount, npoints,nelements).

	! We need to get the integration, which is a scalar(independent of x, y)
    real(DP) :: Dvalues_Integ
    integer :: i, j

	! rInnervector is to save rACvector from previous step
    ! type(t_vectorblock) :: rInnervector
    ! rOutervector is to save rNSvector
    type(t_vectorblock),pointer :: rInnervector
    type(t_vectorScalar) :: rtmpvectorScalar
    real(DP) :: dtime
  
!</output>

  !</subroutine>

  if (present(rcollection)) then

!MCai, we need to modify something here?
! Through rcollection, we can get rcollection%p_rvectorQuickAccess1
    rInnervector=>rcollection%p_rvectorQuickAccess1

    call lsyssc_createVecByDiscr (rInnervector%RvectorBlock(1)%p_rspatialDiscr, &
                          rtmpvectorScalar,.TRUE.,ST_DOUBLE)
    ! Make rtmpvectorScalar to be 0
    call lsyssc_clearvector(rtmpvectorScalar)
  
! We use a technique similar to cal L^1 error to get the integration of f(\phi(t_n))
! fun2d_Target will provided the function values at the cubature points.
    call fun_integral2d_conf(rtmpvectorScalar, PPERR_L1ERROR,Dvalues_Integ,&
                                  rdiscretisation,fun2d_Target,rcollection)

    ! MCai,
    ! Dcoefficients is a scalar
	! We also need to divide this by the volume of the doamin if the domain is not 1.0
!    print *, Dvalues_Integ
    Dcoefficients (:,:,:) = (gamma/eps**2)*Dvalues_Integ


    call lsyssc_releaseVector(rtmpvectorScalar)
  end if

  end subroutine
!
! We need use rcollection to get rACvector, and use rACvector to evaluate
! the function value of \phi^3-\phi at any given (x,y), but the problem here is
! that Dpoints may only corresponds to quadrature points.
  
! ***************************************************************************
  
!<subroutine>
!MCai, the output of this subroutine is Dvalues(npointsPerElement, nelements)
! we through rcollection to get rACvector, we need to cal
! \Phi(x,y)^3 -\Phi(x,y) evaluated at quadrature points.

  subroutine  fun2d_Target(cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTrial,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the postprocessing.
  ! It should return values of the analytical solution (if it is known).
  ! These are compared with the calculated solution to calculate the
  ! error in the X-velocity.
  !
  ! If the analytical solution is unknown, this routine doesn't make sense.
  ! In this case, error analysis should be deactivated in the .DAT files!
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
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsreal.
  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTrial

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
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

! MCai
    type(t_vectorBlock), pointer :: rInnervector

    real(DP), dimension(:,:,:), allocatable  :: DvaluesTmp

   
    ! In a nonstationary simulation, one can get the simulation time
    ! with the quick-access array of the collection.
    if (present(rcollection)) then
      ! MCai,
       rInnervector=>rcollection%p_rvectorQuickAccess1
      ! npoints<===> npointsPerElement.....
	  allocate(DvaluesTmp(1, npointsPerElement, nelements))
 
      Dvalues(:,:) = 0.0_DP
   

      if (cderivative .EQ. DER_FUNC) then
      ! write(*,*)'we need to specify Dvalues(npointsPerElement,nelements)'
      ! fevl_evaluate_cubicpoly is in fe_cub_evaluation.f90
	    call fevl_evaluate_cubicpoly (rInnervector%RvectorBlock(1), &
                 rdomainIntSubset, cderivative, DvaluesTmp, 1)

        Dvalues(:,:)=DvaluesTmp(1, :, :)
	  else
	    write(*,*) 'This code is only for evaluating fun value, not derivative'
      end if

      deallocate(DvaluesTmp)
    end if

  end subroutine

!<subroutine>

 
! ***************************************************************************

END MODULE

!##############################################################################
!# ****************************************************************************
!# <name> newtonderivative </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module implements the Newton derivative of certain functions.
!#
!# Routines in this module:
!#
!# 1.) nder_minMaxProjByCubature
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using a cubature approach.
!#
!# 2.) nder_minMaxProjByMass
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using a modification of the mass matrix.
!#
!# 3.) nder_minMaxProjByApproxDer
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using a finite difference approach in each entry.
!#
!# </purpose>
!##############################################################################

module newtonderivative

  use fsystem
  use genoutput
  use storage
  
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use element
  use derivatives
  use feevaluation
  use scalarpde
  use bilinearformevaluation
  use matrixmodification
  
  implicit none
  
  private
  
  public :: nder_minMaxProjByCubature
  public :: nder_minMaxProjByMass
  public :: nder_minMaxProjByApproxDer
  
contains

  ! ***************************************************************************
  
  !<subroutine>

  subroutine coeff_MinMaxProj (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the Newton derivative of the MinMax projection.
    !
    ! Calculates the values in the cubature points.
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
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    type(t_vectorScalar), pointer :: p_rvector, p_rvectorMin, p_rvectorMax
    real(dp), dimension(:,:), allocatable :: Dfunc, DfuncMin, DfuncMax
    integer(I32) :: celement
    real(DP) :: dwFct, dweight, dwMin, dwMax
    integer :: iel, ipt
    
    ! Get the bounds and the multiplier from the collection
    dwFct = rcollection%DquickAccess(1)
    dweight = rcollection%DquickAccess(2)
    dwMin = rcollection%DquickAccess(3)
    dwMax = rcollection%DquickAccess(4)
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! Lower/upper bound specified?
    nullify(p_rvectorMin)
    nullify(p_rvectorMax)
    
    if (associated(rcollection%p_rvectorQuickAccess2)) then
      p_rvectorMin => rcollection%p_rvectorQuickAccess2%RvectorBlock(1)
    end if

    if (associated(rcollection%p_rvectorQuickAccess3)) then
      p_rvectorMax => rcollection%p_rvectorQuickAccess3%RvectorBlock(1)
    end if
  
    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rvector, &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    ! Now check the function values lambda.
    ! If dwMin*fctmin < dweight*u < dwMax*fctMax, return dweight.
    ! Otherwise, return 0.
    
    ! We have to take care of 4 cases, depending on which parameters
    ! are specified.
    
    if (associated(p_rvectorMin) .and. associated(p_rvectorMax)) then
    
      ! Both functions nonconstant.
      ! Allocate memory for the function values in the cubature points.
      allocate(DfuncMin(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      allocate(DfuncMax(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      
      ! Calculate the function value of the solution vector in all
      ! our cubature points.
      !
      ! WARNING: We assume the element to be the same as rfunction!
      call fevl_evaluate_sim (p_rvectorMin, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMin)

      call fevl_evaluate_sim (p_rvectorMax, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMax)
          
      do iel = 1,ubound(Dcoefficients,3)
        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax*DfuncMax(ipt,iel))) then
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          else
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end if
        end do
      end do

      ! Release memory
      deallocate(DfuncMin)
      deallocate(DfuncMax)
          
    else if (associated(p_rvectorMin)) then
    
      ! Minimum function nonconstant.
      ! Allocate memory for the function values in the cubature points.
      allocate(DfuncMin(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      
      ! Calculate the function value of the solution vector in all
      ! our cubature points.
      !
      ! WARNING: We assume the element to be the same as rfunction!
      call fevl_evaluate_sim (p_rvectorMin, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMin)

      do iel = 1,ubound(Dcoefficients,3)
        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax)) then
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          else
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end if
        end do
      end do

      ! Release memory
      deallocate(DfuncMin)

    else if (associated(p_rvectorMax)) then
    
      ! Max functions nonconstant.
      ! Allocate memory for the function values in the cubature points.
      allocate(DfuncMax(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      
      ! Calculate the function value of the solution vector in all
      ! our cubature points.
      !
      ! WARNING: We assume the element to be the same as rfunction!

      call fevl_evaluate_sim (p_rvectorMax, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMax)
          
      do iel = 1,ubound(Dcoefficients,3)
        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax*DfuncMax(ipt,iel))) then
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          else
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end if
        end do
      end do

      ! Release memory
      deallocate(DfuncMax)

    else
    
      ! Constant Min/Max functions
      do iel = 1,ubound(Dcoefficients,3)
        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax)) then
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          else
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end if
        end do
      end do
    
    end if
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine

  ! ***************************************************************************
  
  !<subroutine>

  subroutine coeff_MinMaxProjColl (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the Newton derivative of the MinMax projection.
    !
    ! Additionally to the computation of the operator, this routine
    ! collects the elements of the active set for later re-assembly.
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
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    type(t_vectorScalar), pointer :: p_rvector, p_rvectorMin, p_rvectorMax
    real(dp), dimension(:,:), allocatable :: Dfunc, DfuncMin, DfuncMax
    integer(I32) :: celement
    real(DP) :: dweight, dwMin, dwMax, dwFct
    integer :: iel, ipt, nptsInactive
    integer, dimension(:), pointer :: p_IelementList
    
    ! Get the bounds and the multiplier from the collection
    dwFct = rcollection%DquickAccess(1)
    dweight = rcollection%DquickAccess(2)
    dwMin = rcollection%DquickAccess(3)
    dwMax = rcollection%DquickAccess(4)
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)
    
    ! Get the element list
    call storage_getbase_int (rcollection%IquickAccess(2),p_IelementList)

    ! Lower/upper bound specified?
    nullify(p_rvectorMin)
    nullify(p_rvectorMax)
    
    if (associated(rcollection%p_rvectorQuickAccess2)) then
      p_rvectorMin => rcollection%p_rvectorQuickAccess2%RvectorBlock(1)
    end if

    if (associated(rcollection%p_rvectorQuickAccess3)) then
      p_rvectorMax => rcollection%p_rvectorQuickAccess3%RvectorBlock(1)
    end if
  
    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rvector, &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    ! Now check the function values lambda.
    ! If dwMin*fctmin < dwFct*u < dwMax*fctMax, return dwFct*dweight.
    ! Otherwise, return 0.
    
    ! We have to take care of 4 cases, depending on which parameters
    ! are specified.
    
    if (associated(p_rvectorMin) .and. associated(p_rvectorMax)) then
    
      ! Both functions nonconstant.
      ! Allocate memory for the function values in the cubature points.
      allocate(DfuncMin(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      allocate(DfuncMax(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      
      ! Calculate the function value of the solution vector in all
      ! our cubature points.
      !
      ! WARNING: We assume the element to be the same as rfunction!
      call fevl_evaluate_sim (p_rvectorMin, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMin)

      call fevl_evaluate_sim (p_rvectorMax, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMax)
          
      do iel = 1,ubound(Dcoefficients,3)

        ! Count the inactive points.
        nptsInactive = 0

        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax*DfuncMax(ipt,iel))) then
            nptsInactive = nptsInactive + 1
          end if
        end do
        
        ! All points inactive? Ok.
        ! Partially active? Remember the element.
        ! Completely active? Return 0 everywhere.
        if (nptsInactive .eq.  ubound(Dcoefficients,2)) then
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          end do
        else
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end do
          if (nptsInactive .gt. 0) then
            rcollection%IquickAccess(3) = rcollection%IquickAccess(3) + 1
            p_IelementList(rcollection%IquickAccess(3)) = &
                rdomainIntSubset%p_Ielements(iel)
          end if
        end if
        
      end do

      ! Release memory
      deallocate(DfuncMin)
      deallocate(DfuncMax)
          
    else if (associated(p_rvectorMin)) then
    
      ! Minimum function nonconstant.
      ! Allocate memory for the function values in the cubature points.
      allocate(DfuncMin(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      
      ! Calculate the function value of the solution vector in all
      ! our cubature points.
      !
      ! WARNING: We assume the element to be the same as rfunction!
      call fevl_evaluate_sim (p_rvectorMin, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMin)

      do iel = 1,ubound(Dcoefficients,3)
        ! Count the inactive points.
        nptsInactive = 0

        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax)) then
            nptsInactive = nptsInactive + 1
          end if
        end do
        
        ! All points inactive? Ok.
        ! Partially active? Remember the element.
        ! Completely active? Return 0 everywhere.
        if (nptsInactive .eq.  ubound(Dcoefficients,2)) then
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          end do
        else
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end do
          if (nptsInactive .gt. 0) then
            rcollection%IquickAccess(3) = rcollection%IquickAccess(3) + 1
            p_IelementList(rcollection%IquickAccess(3)) = &
                rdomainIntSubset%p_Ielements(iel)
          end if
        end if
      end do

      ! Release memory
      deallocate(DfuncMin)

    else if (associated(p_rvectorMax)) then
    
      ! Max functions nonconstant.
      ! Allocate memory for the function values in the cubature points.
      allocate(DfuncMax(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      
      ! Calculate the function value of the solution vector in all
      ! our cubature points.
      !
      ! WARNING: We assume the element to be the same as rfunction!

      call fevl_evaluate_sim (p_rvectorMax, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, DfuncMax)
          
      do iel = 1,ubound(Dcoefficients,3)
        ! Count the inactive points.
        nptsInactive = 0

        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax*DfuncMax(ipt,iel))) then
            nptsInactive = nptsInactive + 1
          end if
        end do
        
        ! All points inactive? Ok.
        ! Partially active? Remember the element.
        ! Completely active? Return 0 everywhere.
        if (nptsInactive .eq.  ubound(Dcoefficients,2)) then
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          end do
        else
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end do
          if (nptsInactive .gt. 0) then
            rcollection%IquickAccess(3) = rcollection%IquickAccess(3) + 1
            p_IelementList(rcollection%IquickAccess(3)) = &
                rdomainIntSubset%p_Ielements(iel)
          end if
        end if
      end do

      ! Release memory
      deallocate(DfuncMax)

    else
    
      ! Constant Min/Max functions
      do iel = 1,ubound(Dcoefficients,3)
        ! Count the inactive points.
        nptsInactive = 0

        do ipt = 1,ubound(Dcoefficients,2)
          ! Check if the dual variable is in the bounds for the control.
          if ((dwFct*Dfunc(ipt,iel) .gt. dwMin) .and. &
              (dwFct*Dfunc(ipt,iel) .lt. dwMax)) then
            nptsInactive = nptsInactive + 1
          end if
        end do
        
        ! All points inactive? Ok.
        ! Partially active? Remember the element.
        ! Completely active? Return 0 everywhere.
        if (nptsInactive .eq.  ubound(Dcoefficients,2)) then
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = dwFct*dweight
          end do
        else
          do ipt = 1,ubound(Dcoefficients,2)
            Dcoefficients(1,ipt,iel) = 0.0_DP
          end do
          if (nptsInactive .gt. 0) then
            rcollection%IquickAccess(3) = rcollection%IquickAccess(3) + 1
            p_IelementList(rcollection%IquickAccess(3)) = &
                rdomainIntSubset%p_Ielements(iel)
          end if
        end if
      end do
    
    end if
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine nder_minMaxProjByCubature (dweight,rmatrix,rcubatureInfo,&
      dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax)
  
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> min( dwmin*rfunctionMin  max(dwmax*rfunctionMax, dweight*rfunction))
  ! If rfunctionMin/rfunctionMax are not specified, they are assumed
  ! to be =1.
  ! The operator is added to rmatrix.
  !
  ! NOTE: THis assumes rfunctionMin and rfunctionMax to be
  ! discretised in the same space as rfunction!
!</description>

!<input>
  ! Weight in front of the operator when being added to rmatrix.
  real(DP), intent(in) :: dweight
  
  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFct
  
  ! An FE function
  type(t_vectorScalar), intent(in) :: rfunction
  
  ! Weight for the lower bound. If rfunctionMin is not specified,
  ! this is the lower bound.
  real(DP), intent(in) :: dwMin

  ! Weight for the upper bound. If rfunctionMax is not specified,
  ! this is the upper bound.
  real(DP), intent(in) :: dwMax

  ! Cubature info structure that defines how to apply cubature.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
  
  ! OPTIONAL: Function specifying the lower bound.
  type(t_vectorScalar), intent(in), optional :: rfunctionMin

  ! OPTIONAL: Function specifying the upper bound.
  type(t_vectorScalar), intent(in), optional :: rfunctionMax
!</input>

!<inputoutput>
  ! The matrix which receives the operator.
  type(t_matrixScalar), intent(inout) :: rmatrix
!<inputoutput>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_bilinearForm) :: rform
    type(t_vectorBlock), target :: rvecFct, rvecMin, rvecMax
    
    ! Set up a bilinear form for the assembly of the
    ! modified mass matrices.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC

    ! In this case, we have nonconstant coefficients.
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff(:) = .false.
    
    ! Prepare a collection structure to be passed to the callback
    ! routine. We attach the function(s) in the quick-access variables
    ! so the callback routine can access it.
    !
    ! All scalar vectors must be converted to block vectors.
    call lsysbl_createVecFromScalar (rfunction,rvecFct)
    rcollection%p_rvectorQuickAccess1 => rvecFct
    
    if (present(rfunctionMin)) then
      call lsysbl_createVecFromScalar (rfunctionMin,rvecMin)
      rcollection%p_rvectorQuickAccess2 => rvecMin
    else
      nullify(rcollection%p_rvectorQuickAccess2)
    end if

    if (present(rfunctionMax)) then
      call lsysbl_createVecFromScalar (rfunctionMax,rvecMax)
      rcollection%p_rvectorQuickAccess3 => rvecMax
    else
      nullify(rcollection%p_rvectorQuickAccess3)
    end if
    
    rcollection%DquickAccess(1) = dwFct
    rcollection%DquickAccess(2) = dweight
    rcollection%DquickAccess(3) = dwMin
    rcollection%DquickAccess(4) = dwMax

    ! Call the assembly routine to calculate the operator.
    call bilf_buildMatrixScalar (rform,.false.,rmatrix,&
        rcubatureInfo,coeff_MinMaxProj,rcollection)    

    ! Release memory
    if (present(rfunctionMin)) then
      call lsysbl_releaseVector (rvecMin)
    end if

    if (present(rfunctionMax)) then
      call lsysbl_releaseVector (rvecMax)
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine nder_minMaxProjByMass (rmassMatrix,dweight,rmatrix,&
      dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax)
      
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> min( dwmin*rfunctionMin  max(dwmax*rfunctionMax, dweight*rfunction))
  ! by using filtering of a mass matrix.
  ! If rfunctionMin/rfunctionMax are not specified, they are assumed
  ! to be =1.
  ! All rows where the corresponding DOFs violate the bounds are
  ! treated as zero.
  !
  ! This routine makes only sense for Lagrangian type finite elements!
!</description>
  
!<input>
  ! Weight in front of the mass matrix when being added to rmatrix.
  real(DP), intent(in) :: dweight
  
  ! Mass matrix.
  type(t_matrixScalar), intent(in) :: rmassMatrix
  
  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFct
  
  ! An FE function
  type(t_vectorScalar), intent(in) :: rfunction
  
  ! Weight for the lower bound. If rfunctionMin is not specified,
  ! this is the lower bound.
  real(DP), intent(in) :: dwMin

  ! Weight for the upper bound. If rfunctionMax is not specified,
  ! this is the upper bound.
  real(DP), intent(in) :: dwMax
  
  ! OPTIONAL: Function specifying the lower bound.
  type(t_vectorScalar), intent(in), optional :: rfunctionMin

  ! OPTIONAL: Function specifying the upper bound.
  type(t_vectorScalar), intent(in), optional :: rfunctionMax
!</input>
  
!<inputoutput>
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>  
  
!</subroutine>
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata,p_DdataMin,p_DdataMax
    integer, dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du
    type(t_matrixScalar) :: rmassCopy
    
    ! Duplicate the mass matrix
    call lsyssc_duplicateMatrix (rmassMatrix,rmassCopy,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    
    ! Get the vector data
    call lsyssc_getbase_double (rfunction,p_Ddata)
    
    ! Lower/upper bound specified?
    nullify(p_DdataMin)
    nullify(p_DdataMax)
    
    if (present(rfunctionMin)) then
      call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
    end if

    if (present(rfunctionMax)) then
      call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
    end if
    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rfunction%NEQ))
    
    ! We have to take care of 4 cases, depending on which parameters
    ! are specified.
    if (associated(p_DdataMin) .and. associated(p_DdataMax)) then
      
      ! Both variable bounds specified
      nviolate = 0
      do i=1,rfunction%NEQ
        du = dwFct*p_Ddata(i)
        if ((du .le. dwMin*p_DdataMin(i)) .or. (du .ge. dwMax*p_DdataMax(i))) then
          nviolate = nviolate + 1
          p_Idofs(nviolate) = i
        end if
      end do

    else if (associated(p_DdataMin)) then
      
      ! Minimum specified
      nviolate = 0
      do i=1,rfunction%NEQ
        du = dwFct*p_Ddata(i)
        if ((du .le. dwMin*p_DdataMin(i)) .or. (du .ge. dwMax)) then
          nviolate = nviolate + 1
          p_Idofs(nviolate) = i
        end if
      end do
    
    else if (associated(p_DdataMax)) then
    
      ! Maximum specified
      nviolate = 0
      do i=1,rfunction%NEQ
        du = dwFct*p_Ddata(i)
        if ((du .le. dwMin) .or. (du .ge. dwMax*p_DdataMax(i))) then
          nviolate = nviolate + 1
          p_Idofs(nviolate) = i
        end if
      end do
    
    else
    
      ! None specified. Constant bounds.
      nviolate = 0
      do i=1,rfunction%NEQ
        du = dwFct*p_Ddata(i)
        if ((du .le. dwMin) .or. (du .ge. dwMax)) then
          nviolate = nviolate + 1
          p_Idofs(nviolate) = i
        end if
      end do
     
    end if
    
    if (nviolate .gt. 0) then
      ! Filter the matrix
      call mmod_replaceLinesByZero (rmassCopy,p_Idofs(1:nviolate))
    end if
    
    ! Sum up
    call lsyssc_matrixLinearComb (rmassCopy,rmatrix,dwFct*dweight,1.0_DP,&
      .false.,.false.,.true.,.true.)

    ! Release memory
    call lsyssc_releaseMatrix (rmassCopy)
    deallocate(p_Idofs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine nder_minMaxProjByApproxDer (dweight,rmatrix,dh,&
      dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax)
      
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> min( dwmin*rfunctionMin  max(dwmax*rfunctionMax, dweight*rfunction))
  ! by using an approximative derivative based on a finite difference
  ! approach in each entry.
  ! If rfunctionMin/rfunctionMax are not specified, they are assumed
  ! to be =1.
  !
  ! This routine makes only sense for Lagrangian type finite elements!
  ! This routine only works for Format-9 matrices!
!</description>
  
  
!<input>
  ! Weight in front of the mass matrix when being added to rmatrix.
  real(DP), intent(in) :: dweight
  
  ! Finite-difference step length for the approximative derivative.
  real(dp), intent(in) :: dh
  
  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFct
  
  ! An FE function
  type(t_vectorScalar), intent(in) :: rfunction
  
  ! Weight for the lower bound. If rfunctionMin is not specified,
  ! this is the lower bound.
  real(DP), intent(in) :: dwMin

  ! Weight for the upper bound. If rfunctionMax is not specified,
  ! this is the upper bound.
  real(DP), intent(in) :: dwMax
  
  ! OPTIONAL: Function specifying the lower bound.
  type(t_vectorScalar), intent(in), optional :: rfunctionMin

  ! OPTIONAL: Function specifying the upper bound.
  type(t_vectorScalar), intent(in), optional :: rfunctionMax
!</input>
  
!<inputoutput>
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>  
  
!</subroutine>
  
    ! The projection operator is given by:
    !
    !          a, if u <= a
    !  P(u) =  u, if a <= u <= b
    !          b, if u >= b
    !
    ! The Frechet derivative of this operator can be calculated by
    !
    !   (P(u+h)-P(u-h))/2 = DP(u)h
    !
    ! with h being an arbitrary function <> 0.
    ! Rewriting this in a discrete sense yields
    !
    !   ( P(u + h e_i) - P(u - h e_i) ) / (2h)  =  DP(u) h e_i
    !
    ! with u being the vector of the corresponding FE function.
    ! Let us denote the matrix B:=DP(u), then we obtain by this formula:
    !
    !   B_ij  =  [ ( P(u + h e_j) - P(u - h e_j) ) / (2h) ]_i
    !
    ! If we now treat the P()-operator coponent-wise instead of vector
    ! wise, this means:
    !
    !   B_ij  =  ( P(u_j+h) - P(u_j-h) ) / (2h)   , i=j
    !            0                                , otherwise
      
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Dmatrix
    real(dp), dimension(:), pointer :: p_Ddata,p_DdataMin,p_DdataMax
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: i
    real(DP) :: du1,du2
    
    ! Get the matrix/vector data
    call lsyssc_getbase_double (rmatrix,p_Dmatrix)
    call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

    call lsyssc_getbase_double (rfunction,p_Ddata)
    
    ! Lower/upper bound specified?
    nullify(p_DdataMin)
    nullify(p_DdataMax)
    
    if (present(rfunctionMin)) then
      call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
    end if

    if (present(rfunctionMax)) then
      call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
    end if
    
    ! We have 4 cases depending on the specified parameters.
    if (present(rfunctionMin) .and. present(rfunctionMax)) then
    
      ! Both bounds present
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax*p_DdataMax(i),max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) - dh))
        du2 = -min(dwMax*p_DdataMax(i),max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) + dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    else if (present(rfunctionMin)) then
    
      ! Lower bound present
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax,max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) - dh))
        du2 = -min(dwMax,max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) + dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    else if (present(rfunctionMax)) then
    
      ! Upper bound present
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax*p_DdataMax(i),max(dwMin,dwFct*p_Ddata(i) - dh))
        du2 = -min(dwMax*p_DdataMax(i),max(dwMin,dwFct*p_Ddata(i) + dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    else
    
      ! No bounds present. Constant bounds.
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax,max(dwMin,dwFct*p_Ddata(i) - dh))
        du2 = -min(dwMax,max(dwMin,dwFct*p_Ddata(i) + dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    end if

  end subroutine

end module

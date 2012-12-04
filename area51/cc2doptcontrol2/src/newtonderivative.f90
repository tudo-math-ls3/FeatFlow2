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
!# 1.) nwder_minMaxProjByCubature
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using a cubature approach.
!#
!# 2.) nwder_minMaxProjByAdaptCub
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using an adaptive cubature approach.
!#
!# 3.) nwder_minMaxProjByMass
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using a modification of the mass matrix.
!#
!# 4.) nwder_minMaxProjByApproxDer
!#     -> Calculates the Newton derivative of the minmax-Projection
!#        using a finite difference approach in each entry.
!#
!# 5.) nwder_rhsMinMaxProjByCubature
!#     -> Assemble a RHS vector of the minmax-Projection
!#        using a cubature approach.
!#
!# 6.) nwder_rhsMinMaxProjByMass
!#     -> Assemble a RHS vector of the minmax-Projection
!#        using a DOF based approach plus multiplication of a mass matrix.
!#
!# 7.) nwder_applyMinMaxProjByDof
!#     -> Applied the minmax-Projection using a DOF based approach.
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
  use linearformevaluation
  use bilinearformevaluation
  use matrixmodification
  
  implicit none
  
  private
  
  public :: nwder_minMaxProjByCubature
  public :: nwder_minMaxProjByAdaptCub
  public :: nwder_minMaxProjByMass
  public :: nwder_minMaxProjByApproxDer
  public :: nwder_rhsMinMaxProjByCubature
  public :: nwder_rhsMinMaxProjByMass
  public :: nwder_applyMinMaxProjByDof
  
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
    type(t_vectorScalar), pointer :: p_rvectorShift
    real(dp), dimension(:,:), allocatable :: Dfunc, DfuncMin, DfuncMax, Dshift
    integer(I32) :: celement
    real(DP) :: dwFct, dweight, dwMin, dwMax, dwShift
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
    nullify(p_rvectorShift)
    
    if (associated(rcollection%p_rvectorQuickAccess2)) then
      p_rvectorMin => rcollection%p_rvectorQuickAccess2%RvectorBlock(1)
    end if

    if (associated(rcollection%p_rvectorQuickAccess3)) then
      p_rvectorMax => rcollection%p_rvectorQuickAccess3%RvectorBlock(1)
    end if
  
    if (associated(rcollection%p_rvectorQuickAccess3)) then
      p_rvectorShift => rcollection%p_rvectorQuickAccess4%RvectorBlock(1)
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
    
    ! If a 'shift' function is present, shift the function.
    if (associated(p_rvectorShift)) then
      allocate(Dshift(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))    

      call fevl_evaluate_sim (p_rvectorShift, &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dshift)
          
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .lt. dwMax*DfuncMax(ipt,iel))) then
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .lt. dwMax)) then
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .gt. dwMin) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .lt. dwMax*DfuncMax(ipt,iel))) then
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .gt. dwMin) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel)) .lt. dwMax)) then
              Dcoefficients(1,ipt,iel) = dwFct*dweight
            else
              Dcoefficients(1,ipt,iel) = 0.0_DP
            end if
          end do
        end do
      
      end if
    
      deallocate (Dshift)

    else
        
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift) .lt. dwMax*DfuncMax(ipt,iel))) then
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift) .gt. dwMin*DfuncMin(ipt,iel)) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift) .lt. dwMax)) then
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift) .gt. dwMin) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift) .lt. dwMax*DfuncMax(ipt,iel))) then
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
            if (((dwFct*Dfunc(ipt,iel) + dwShift) .gt. dwMin) .and. &
                ((dwFct*Dfunc(ipt,iel) + dwShift) .lt. dwMax)) then
              Dcoefficients(1,ipt,iel) = dwFct*dweight
            else
              Dcoefficients(1,ipt,iel) = 0.0_DP
            end if
          end do
        end do
      
      end if

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
    ! This routine returns the coefficients for the Newton derivative,
    ! similar to coeff_MinMaxProj. However, elements which cross
    ! the border of the active set are associated to the active set.
    ! For these elements, this routine collects the elements,
    ! which can then be assembled in a second step with a different
    ! cubature formula.
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
    type(t_vectorScalar), pointer :: p_rvectorShift
    real(dp), dimension(:,:), allocatable :: Dfunc, DfuncMin, DfuncMax, Dshift
    integer(I32) :: celement
    real(DP) :: dweight, dwMin, dwMax, dwFct, dwShift
    integer :: iel, ipt, nptsInactive
    integer, dimension(:), pointer :: p_IelementList
    
    ! Get the bounds and the multiplier from the collection
    dwFct = rcollection%DquickAccess(1)
    dweight = rcollection%DquickAccess(2)
    dwMin = rcollection%DquickAccess(3)
    dwMax = rcollection%DquickAccess(4)
    dwShift = rcollection%DquickAccess(5)
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)
    
    ! Get the element buffer
    call storage_getbase_int (rcollection%IquickAccess(1),p_IelementList)

    ! Lower/upper bound specified?
    nullify(p_rvectorMin)
    nullify(p_rvectorMax)
    nullify(p_rvectorShift)
    
    if (associated(rcollection%p_rvectorQuickAccess2)) then
      p_rvectorMin => rcollection%p_rvectorQuickAccess2%RvectorBlock(1)
    end if

    if (associated(rcollection%p_rvectorQuickAccess3)) then
      p_rvectorMax => rcollection%p_rvectorQuickAccess3%RvectorBlock(1)
    end if

    if (associated(rcollection%p_rvectorQuickAccess3)) then
      p_rvectorShift => rcollection%p_rvectorQuickAccess4%RvectorBlock(1)
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
            rcollection%IquickAccess(2) = rcollection%IquickAccess(2) + 1
            p_IelementList(rcollection%IquickAccess(2)) = &
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
            rcollection%IquickAccess(2) = rcollection%IquickAccess(2) + 1
            p_IelementList(rcollection%IquickAccess(2)) = &
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
            rcollection%IquickAccess(2) = rcollection%IquickAccess(2) + 1
            p_IelementList(rcollection%IquickAccess(2)) = &
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
        
        ! All points inactive? Ok, return the operator.
        ! Partially active? Remember the element and return 0.
        ! Element will be reassembled later.
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
            rcollection%IquickAccess(2) = rcollection%IquickAccess(2) + 1
            p_IelementList(rcollection%IquickAccess(2)) = &
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

  subroutine nwder_minMaxProjByCubature (dweight,rmatrix,rcubatureInfo,&
      dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax,&
      dwShift,rfunctionShift)
  
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> 
  !     min( dwmin*rfunctionMin , max(dwmax*rfunctionMax, dwFct*rfunction + dwShift*rfunctionShift))
  ! If rfunctionMin/rfunctionMax/rfunctionShift are not specified, 
  ! they are assumed to be =1.
  ! The operator is added to rmatrix.
  !
  ! NOTE: This assumes rfunctionMin and rfunctionMax to be
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
  
  ! OPTIONAL: Weight for the shift function.
  ! If not present, =0 is assumed.
  real(DP), optional :: dwShift
  
  ! OPTIONAL: Shift function
  type(t_vectorScalar), intent(in), optional :: rfunctionShift
!</input>

!<inputoutput>
  ! The matrix which receives the operator.
  type(t_matrixScalar), intent(inout) :: rmatrix
!<inputoutput>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_bilinearForm) :: rform
    type(t_vectorBlock), target :: rvecFct, rvecMin, rvecMax, rvecShift
    
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

    if (present(rfunctionShift)) then
      call lsysbl_createVecFromScalar (rfunctionShift,rvecShift)
      rcollection%p_rvectorQuickAccess4 => rvecShift
    else
      nullify(rcollection%p_rvectorQuickAccess4)
    end if
    
    rcollection%DquickAccess(1) = dwFct
    rcollection%DquickAccess(2) = dweight
    rcollection%DquickAccess(3) = dwMin
    rcollection%DquickAccess(4) = dwMax
    rcollection%DquickAccess(5) = 0.0_DP
    if (present(dwShift)) then
      rcollection%DquickAccess(5) = dwShift
    end if

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

    if (present(rfunctionShift)) then
      call lsysbl_releaseVector (rvecShift)
    end if
    
  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine nwder_minMaxProjByAdaptCub (dweight,rmatrix,rcubatureInfo,&
      rcubatureInfoAdapt,dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax)
  
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> min( dwmin*rfunctionMin , max(dwmax*rfunctionMax, dwFct*rfunction))
  ! If rfunctionMin/rfunctionMax are not specified, they are assumed
  ! to be =1.
  ! The operator is added to rmatrix. The calculation is applied
  ! in two steps, using adaptive cubature.
  !
  ! NOTE: This assumes rfunctionMin and rfunctionMax to be
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

  ! Cubature info structure that defines how to apply cubature
  ! in the active/inactive set
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo

  ! Cubature info structure that defines how to apply cubature
  ! on the border of the active set
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfoAdapt
  
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
    integer :: nelements, ielemHandle
    type(t_bilfMatrixAssembly) :: rmatrixAssembly
    integer(I32) :: ccubType, celement
    integer, dimension(:), pointer :: p_IelementList
    
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

    ! Shift function currently not supported. Set to zero.
    nullify(rcollection%p_rvectorQuickAccess4)
    rcollection%DquickAccess(5) = 0.0

    ! Create an array that saves all elements on the border of the active set.
    nelements = rmatrix%p_rspatialDiscrTrial%p_rtriangulation%NEL
    call storage_new ('', 'Ielements', nelements, ST_INT, ielemHandle, &
        ST_NEWBLOCK_NOINIT)

    ! The IquickAccess(1) element saves the handle of the element list.
    ! IquickAccess(2) saves how many elements are collected.
    rcollection%IquickAccess(1) = ielemHandle
    rcollection%IquickAccess(2) = 0

    ! Call the assembly routine to calculate the operator.
    ! Simultaneously, collect the elements on the border of the
    ! active set. Do not calculate anything on the border.
    call bilf_buildMatrixScalar (rform,.false.,rmatrix,&
        rcubatureInfo,coeff_MinMaxProjColl,rcollection)
        
    ! In a second step, assemble a submesh matrix on the elements in the list
    ! with the extended cubature formula.
    ! Note: Up to now, this works only for uniform meshes!
    if (rcollection%IquickAccess(3) .gt. 0) then
      
      ! Get the underlying element
      celement = rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement
      
      ! Get the cubature formula.
      call spdiscr_getStdDiscrInfo(1,rcubatureInfoAdapt,&
          rmatrix%p_rspatialDiscrTrial,ccubature=ccubType)
      
      ! Assemble the matrix, now using the standard assembly callback
      ! routine coeff_MinMaxProj.
      call storage_getbase_int(ielemhandle,p_IelementList)
      call bilf_initAssembly(rmatrixAssembly,rform,celement,celement,ccubType)
      call bilf_assembleSubmeshMatrix9(rmatrixAssembly,rmatrix,&
          p_IelementList(1:rcollection%IquickAccess(2)),coeff_MinMaxProj,rcollection)
      call bilf_doneAssembly(rmatrixAssembly)
      
    end if


    ! Release memory
    call storage_free (ielemHandle)

    if (present(rfunctionMin)) then
      call lsysbl_releaseVector (rvecMin)
    end if

    if (present(rfunctionMax)) then
      call lsysbl_releaseVector (rvecMax)
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine nwder_minMaxProjByMass (rmassMatrix,dweight,rmatrix,bclear,&
      dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax,&
      dwShift, rfunctionShift)
      
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> 
  !      min( dwmin*rfunctionMin  max(dwmax*rfunctionMax, dwFct*rfunction + dwShift*rfunctionShift))
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
  
  ! If set to TRUE, rmatrix is cleared in advance.
  logical, intent(in) :: bclear
  
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

  ! OPTIONAL: Weight for the shift function.
  ! If not present, =0 is assumed.
  real(DP), optional :: dwShift
  
  ! OPTIONAL: Shift function
  type(t_vectorScalar), intent(in), optional :: rfunctionShift
!</input>
  
!<inputoutput>
  ! Matrix which receives the operator.
  type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>  
  
!</subroutine>
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata,p_DdataMin,p_DdataMax,p_DdataShift
    integer, dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du, dushift
    type(t_matrixScalar), target :: rmassCopy
    type(t_matrixScalar), pointer :: p_rmatrix
    
    ! Duplicate the mass matrix or create a new one --
    ! depending on whether to "clear" the matrix or not.
    ! "Clear" is implemented here by overwriting the entries with those
    ! of the mass matrix.
    if (bclear) then
      call lsyssc_duplicateMatrix (rmassMatrix,rmatrix,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      p_rmatrix => rmatrix
      
      ! Mark rmassCopy as "not used".
      rmassCopy%NA = 0
    else
      call lsyssc_duplicateMatrix (rmassMatrix,rmassCopy,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      p_rmatrix => rmassCopy
    end if
    
    ! Get the vector data
    call lsyssc_getbase_double (rfunction,p_Ddata)
    
    ! Lower/upper bound specified?
    nullify(p_DdataMin)
    nullify(p_DdataMax)
    nullify(p_DdataShift)
    
    if (present(rfunctionMin)) then
      call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
    end if

    if (present(rfunctionMax)) then
      call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
    end if

    ! Shift function =0 by default    
    dushift = 0.0_DP
    if (present(dwShift)) then
      duShift = dwShift
    end if

    if (present(rfunctionShift)) then
      call lsyssc_getbase_double (rfunctionShift,p_DdataShift)
    end if

    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rfunction%NEQ))
    
    if (associated(p_DdataShift)) then

      ! We have to take care of 4 cases, depending on which parameters
      ! are specified.
      if (associated(p_DdataMin) .and. associated(p_DdataMax)) then
        
        ! Both variable bounds specified
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift*p_DdataShift(i)
          if ((du .le. dwMin*p_DdataMin(i)) .or. (du .ge. dwMax*p_DdataMax(i))) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do

      else if (associated(p_DdataMin)) then
        
        ! Minimum specified
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift*p_DdataShift(i)
          if ((du .le. dwMin*p_DdataMin(i)) .or. (du .ge. dwMax)) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do
      
      else if (associated(p_DdataMax)) then
      
        ! Maximum specified
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift*p_DdataShift(i)
          if ((du .le. dwMin) .or. (du .ge. dwMax*p_DdataMax(i))) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do
      
      else
      
        ! None specified. Constant bounds.
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift*p_DdataShift(i)
          if ((du .le. dwMin) .or. (du .ge. dwMax)) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do
       
      end if
      
      if (nviolate .gt. 0) then
        ! Filter the matrix
        call mmod_replaceLinesByZero (p_rmatrix,p_Idofs(1:nviolate))
      end if
      
      ! If we created a temporary mass matrix, sum up to the original
      ! one. Otherwise we are done after scaling the entries.
      if (.not. bclear) then
        ! Sum up. p_rmatrix and rmassCopy coincide.
        call lsyssc_matrixLinearComb (rmassCopy,rmatrix,dwFct*dweight,1.0_DP,&
          .false.,.false.,.true.,.true.)
          
        ! Release the copy
        call lsyssc_releaseMatrix (rmassCopy)
      else
        call lsyssc_scaleMatrix (rmatrix,dwFct*dweight)
      end if

    else
      ! We have to take care of 4 cases, depending on which parameters
      ! are specified.
      if (associated(p_DdataMin) .and. associated(p_DdataMax)) then
        
        ! Both variable bounds specified
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift
          if ((du .le. dwMin*p_DdataMin(i)) .or. (du .ge. dwMax*p_DdataMax(i))) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do

      else if (associated(p_DdataMin)) then
        
        ! Minimum specified
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift
          if ((du .le. dwMin*p_DdataMin(i)) .or. (du .ge. dwMax)) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do
      
      else if (associated(p_DdataMax)) then
      
        ! Maximum specified
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift
          if ((du .le. dwMin) .or. (du .ge. dwMax*p_DdataMax(i))) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do
      
      else
      
        ! None specified. Constant bounds.
        nviolate = 0
        do i=1,rfunction%NEQ
          du = dwFct*p_Ddata(i) + duShift
          if ((du .le. dwMin) .or. (du .ge. dwMax)) then
            nviolate = nviolate + 1
            p_Idofs(nviolate) = i
          end if
        end do
       
      end if
      
      if (nviolate .gt. 0) then
        ! Filter the matrix
        call mmod_replaceLinesByZero (p_rmatrix,p_Idofs(1:nviolate))
      end if
      
      ! If we created a temporary mass matrix, sum up to the original
      ! one. Otherwise we are done after scaling the entries.
      if (.not. bclear) then
        ! Sum up. p_rmatrix and rmassCopy coincide.
        call lsyssc_matrixLinearComb (rmassCopy,rmatrix,dwFct*dweight,1.0_DP,&
          .false.,.false.,.true.,.true.)
          
        ! Release the copy
        call lsyssc_releaseMatrix (rmassCopy)
      else
        call lsyssc_scaleMatrix (rmatrix,dwFct*dweight)
      end if
    
    end if

    ! Release memory
    deallocate(p_Idofs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine nwder_minMaxProjByApproxDer (dweight,rmatrix,dh,&
      dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax)
      
!<description>
  ! Assembles the Newton derivative of the operator
  !   rfunction -> min( dwmin*rfunctionMin  max(dwmax*rfunctionMax, dwFct*rfunction))
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
        
        du1 = -min(dwMax*p_DdataMax(i),max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) + dh))
        du2 = -min(dwMax*p_DdataMax(i),max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) - dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    else if (present(rfunctionMin)) then
    
      ! Lower bound present
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax,max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) + dh))
        du2 = -min(dwMax,max(dwMin*p_DdataMin(i),dwFct*p_Ddata(i) - dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    else if (present(rfunctionMax)) then
    
      ! Upper bound present
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax*p_DdataMax(i),max(dwMin,dwFct*p_Ddata(i) + dh))
        du2 = -min(dwMax*p_DdataMax(i),max(dwMin,dwFct*p_Ddata(i) - dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    else
    
      ! No bounds present. Constant bounds.
    
      ! Loop through the diagonal entries
      do i=1,rmatrix%NEQ

        ! Calculate the diagonal entry, that's it
        
        du1 = -min(dwMax,max(dwMin,dwFct*p_Ddata(i) + dh))
        du2 = -min(dwMax,max(dwMin,dwFct*p_Ddata(i) - dh))
        
        p_Dmatrix(p_Kdiagonal(i)) = dweight*(du1-du2)/(2.0_DP*dh)
      
      end do

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_rhsMinMaxProj (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! For distributed control.
    ! Coefficients in the bilinear form of the projection operator.
    ! Variable constraints, given by a FEM function.
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
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
      ! local variables
    type(t_vectorBlock), pointer :: p_rvector,p_rvectorUmin,p_rvectorUmax,p_rvectorShift
    type(t_vectorBlock), pointer :: p_rvectorViol
    real(dp), dimension(:,:), allocatable :: Dfunc,Dshift,Dviol
    real(dp), dimension(:), allocatable :: Dumin,Dumax
    integer, dimension(:), allocatable :: Ielements
    integer(I32) :: celement
    integer :: ipt, iel
    real(DP) :: dwFct,dweight,dwMin,dwMax, dwShift, dwViol

    ! Get constant data    
    dwFct = rcollection%DquickAccess(1)
    dweight = rcollection%DquickAccess(2)
    dwMin = rcollection%DquickAccess(3)
    dwMax = rcollection%DquickAccess(4)
    dwShift = rcollection%DquickAccess(5)
    dwViol = rcollection%DquickAccess(6)
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1
    p_rvectorUmin => rcollection%p_rvectorQuickAccess2
    p_rvectorUmax => rcollection%p_rvectorQuickAccess3
    p_rvectorShift => rcollection%p_rvectorQuickAccess4
    p_rvectorViol => rcollection%p_rvectorQuickAccess5

    ! Allocate memory for the function values in the cubature points.
    ! Function value, minimum and maximum.
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))

    ! Allocate temp memory for element hints
    allocate(Ielements(ubound(Dcoefficients,2)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector we evaluate!
    
    celement = rdomainIntSubset%celement
    
    call fevl_evaluate_sim (p_rvector%RvectorBlock(1), &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc(:,:))

    ! Evaluate the 'violation' function
    allocate(Dviol(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    call fevl_evaluate_sim (p_rvectorViol%RvectorBlock(1), &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dviol(:,:))
    
    ! Now check the function values lambda.
    ! Return the projected control in every cubature point.
    
    if ((dwShift .ne. 0.0_DP) .and. associated (p_rvectorShift)) then
    
      ! Evaluate the 'shift' function
      allocate(Dshift(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
      call fevl_evaluate_sim (p_rvectorShift%RvectorBlock(1), &
          rdomainIntSubset%p_revalElementSet, &
          celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dshift(:,:))

      ! Nonconstant upper/lower bound given?
      ! Choose the appropriate loop...
      if (associated(p_rvectorUmin) .and. associated(p_rvectorUmax)) then
        
        ! Allocate temp memory for min/max values for u.
        allocate(Dumin(ubound(Dcoefficients,3)))
        allocate(Dumax(ubound(Dcoefficients,3)))

        do iel = 1,ubound(Dcoefficients,3)
        
          ! Evaluate min and max value of the control in the cubature points
          ! on the current element.
          Ielements(:) = rdomainIntSubset%p_Ielements(iel)

          call fevl_evaluate (DER_FUNC, Dumin, p_rvectorUmin%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          call fevl_evaluate (DER_FUNC, Dumax, p_rvectorUmax%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
          
            call projection (dwViol*Dviol(ipt,iel) + dwShift*Dshift(ipt,iel),&
                dwMin*Dumin(ipt), dwMax*Dumax(ipt),&
                dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel),&
                Dcoefficients(1,ipt,iel))

          end do
        end do

        deallocate(Dumin)
        deallocate(Dumax)
      
      else if (associated(p_rvectorUmin)) then
        
        ! Allocate temp memory for min/max values for u.
        allocate(Dumin(ubound(Dcoefficients,3)))

        do iel = 1,ubound(Dcoefficients,3)
        
          ! Evaluate min and max value of the control in the cubature points
          ! on the current element.
          Ielements(:) = rdomainIntSubset%p_Ielements(iel)

          call fevl_evaluate (DER_FUNC, Dumin, p_rvectorUmin%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift*Dshift(ipt,iel),&
                dwMin*Dumin(ipt), dwMax,&
                dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel),&
                Dcoefficients(1,ipt,iel))
          end do
        end do

        deallocate(Dumin)
      
      else if (associated(p_rvectorUmax)) then
        
        ! Allocate temp memory for min/max values for u.
        allocate(Dumax(ubound(Dcoefficients,3)))

        do iel = 1,ubound(Dcoefficients,3)
        
          ! Evaluate min and max value of the control in the cubature points
          ! on the current element.
          Ielements(:) = rdomainIntSubset%p_Ielements(iel)

          call fevl_evaluate (DER_FUNC, Dumax, p_rvectorUmax%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift*Dshift(ipt,iel),&
                dwMin, dwMax*Dumax(ipt),&
                dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel),&
                Dcoefficients(1,ipt,iel))
          end do
        end do
        
        deallocate(Dumax)

      else

        ! Constant bounds
        do iel = 1,ubound(Dcoefficients,3)
          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift*Dshift(ipt,iel),&
                dwMin, dwMax,&
                dwFct*Dfunc(ipt,iel) + dwShift*Dshift(ipt,iel),&
                Dcoefficients(1,ipt,iel))
          end do
        end do

      end if

      deallocate (DShift)

    else
    
      ! Nonconstant upper/lower bound given?
      ! Choose the appropriate loop...
      if (associated(p_rvectorUmin) .and. associated(p_rvectorUmax)) then
        
        ! Allocate temp memory for min/max values for u.
        allocate(Dumin(ubound(Dcoefficients,3)))
        allocate(Dumax(ubound(Dcoefficients,3)))

        do iel = 1,ubound(Dcoefficients,3)
        
          ! Evaluate min and max value of the control in the cubature points
          ! on the current element.
          Ielements(:) = rdomainIntSubset%p_Ielements(iel)

          call fevl_evaluate (DER_FUNC, Dumin, p_rvectorUmin%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          call fevl_evaluate (DER_FUNC, Dumax, p_rvectorUmax%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift,&
                dwMin*Dumin(ipt), dwMax*Dumax(ipt),&
                dwFct*Dfunc(ipt,iel) + dwShift,&
                Dcoefficients(1,ipt,iel))
          end do
        end do

        deallocate(Dumin)
        deallocate(Dumax)
      
      else if (associated(p_rvectorUmin)) then
        
        ! Allocate temp memory for min/max values for u.
        allocate(Dumin(ubound(Dcoefficients,3)))

        do iel = 1,ubound(Dcoefficients,3)
        
          ! Evaluate min and max value of the control in the cubature points
          ! on the current element.
          Ielements(:) = rdomainIntSubset%p_Ielements(iel)

          call fevl_evaluate (DER_FUNC, Dumin, p_rvectorUmin%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift,&
                dwMin*Dumin(ipt), dwMax,&
                dwFct*Dfunc(ipt,iel) + dwShift,&
                Dcoefficients(1,ipt,iel))
          end do
        end do

        deallocate(Dumin)
      
      else if (associated(p_rvectorUmax)) then
        
        ! Allocate temp memory for min/max values for u.
        allocate(Dumax(ubound(Dcoefficients,3)))

        do iel = 1,ubound(Dcoefficients,3)
        
          ! Evaluate min and max value of the control in the cubature points
          ! on the current element.
          Ielements(:) = rdomainIntSubset%p_Ielements(iel)

          call fevl_evaluate (DER_FUNC, Dumax, p_rvectorUmax%RvectorBlock(1), &
              Dpoints(:,:,iel), IelementsHint=Ielements)

          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift,&
                dwMin, dwMax*Dumax(ipt),&
                dwFct*Dfunc(ipt,iel) + dwShift,&
                Dcoefficients(1,ipt,iel))
          end do
        end do
        
        deallocate(Dumax)

      else

        ! Constant bounds
        do iel = 1,ubound(Dcoefficients,3)
          ! Calculate the projection in the cubature points
          do ipt = 1,ubound(Dcoefficients,2)
            call projection (dwViol*Dviol(ipt,iel) + dwShift,&
                dwMin, dwMax,&
                dwFct*Dfunc(ipt,iel) + dwShift,&
                Dcoefficients(1,ipt,iel))
          end do
        end do

      end if
      
    end if
    
    ! Release memory
    deallocate(Ielements)
    deallocate(Dfunc)
  
  contains

    subroutine projection (dtest,dmin,dmax,din,dout)
    
    ! The actual projection.
    ! Returns dmin/dmax where dtest violates the bounds,
    ! or din where it does not.
    
    real(DP), intent(in) :: dtest,dmin,dmax,din
    real(DP), intent(out) :: dout
    
      if (dtest .lt. dmin) then
        dout = dmin
      else if (dtest .gt. dmax) then
        dout = dmax
      else
        dout = din
      end if
    
    end subroutine

  end subroutine
  
! ***************************************************************************

!<subroutine>

  subroutine nwder_rhsMinMaxProjByCubature (dweight,rvector,rcubatureInfo,&
      dwFctViolate,rfunctionViolate,dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax,&
      dwShift, rfunctionShift)
  
!<description>
  ! Assembles the linear form
  !   ( dweight*rfunction , test)
  ! for the projection function
  !      dwmax*rfunctionMax + dwShift*rfunctionShift
  !            where dwFctViolate*rfunctionViolate + dwShift*rfunctionShift > dwmax*rfunctionMax
  !      dwmin*rfunctionMin + dwShift*rfunctionShift
  !            where dwFctViolate*rfunctionViolate + dwShift*rfunctionShift < dwmin*rfunctionMin
  !      dwFct*rfunction
  !            elsewhere
  ! If rfunctionMin/rfunctionMax are not specified, they are assumed
  ! to be =1.
  ! The operator is added to rvector.
  !
  ! NOTE: This assumes rfunctionMin and rfunctionMax to be
  ! discretised in the same space as rfunction!
!</description>

!<input>
  ! Weight in front of the operator when being added to rmatrix.
  real(DP), intent(in) :: dweight
  
  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFctViolate
  
  ! An FE function specifying where the projection should be
  ! applied. Whereever dwFctViolate*rfunctionViolate violates the bounds,
  ! the bounds are assumed.
  type(t_vectorScalar), intent(in) :: rfunctionViolate

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

  ! OPTIONAL: Weight for the shift function.
  ! If not present, =0 is assumed.
  real(DP), optional :: dwShift
  
  ! OPTIONAL: Shift function
  type(t_vectorScalar), intent(in), optional :: rfunctionShift
!</input>

!<inputoutput>
  ! The vector which receives the projection.
  type(t_vectorScalar), intent(inout) :: rvector
!<inputoutput>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type (t_linearForm) :: rlinform
    type(t_vectorBlock), target :: rvecFct, rvecMin, rvecMax, rvecShift, rvecViol
    
    ! Prepare a linearform-assembly.
    rlinform%itermCount = 1
    rlinform%Dcoefficients(1) = 1.0_DP
    rlinform%Idescriptors(1) = DER_FUNC

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
    rcollection%DquickAccess(5) = 0.0_DP
    
    if (present(dwShift)) then
      rcollection%DquickAccess(5) = dwShift
    end if

    rcollection%DquickAccess(6) = dwFctViolate

    if (present(rfunctionShift)) then
      call lsysbl_createVecFromScalar (rfunctionShift,rvecShift)
      rcollection%p_rvectorQuickAccess4 => rvecShift
    else
      nullify(rcollection%p_rvectorQuickAccess4)
    end if

    call lsysbl_createVecFromScalar (rfunctionViolate,rvecViol)
    rcollection%p_rvectorQuickAccess5 => rvecViol
    
    ! Assemble the vector
    call linf_buildVectorScalar(rlinform,.false.,rvector,rcubatureInfo,&
        coeff_rhsMinMaxProj,rcollection)

    ! Release memory
    if (present(rfunctionMin)) then
      call lsysbl_releaseVector (rvecMin)
    end if

    if (present(rfunctionMax)) then
      call lsysbl_releaseVector (rvecMax)
    end if

    call lsysbl_releaseVector (rvecViol)
    
  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine nwder_rhsMinMaxProjByMass (dweight,rvector,rmassMatrix,rvectorTemp,&
      dwFctViolate,rfunctionViolate,dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax,&
      dwShift,rfunctionShift)
  
!<description>
  ! Calculates 
  !   rvector = rvector + dweight * MassMatrix * rfunction
  ! for the projection function
  !   rfunction -> 
  !      dwmax*rfunctionMax + dwShift*rfunctionShift
  !            where dwFctViolate*rfunctionViolate + dwShift*rfunctionShift > dwmax*rfunctionMax
  !      dwmin*rfunctionMin + dwShift*rfunctionShift
  !            where dwFctViolate*rfunctionViolate + dwShift*rfunctionShift < dwmin*rfunctionMin
  !      dwFct*rfunction
  !            elsewhere
  ! If rfunctionMin/rfunctionMax/rfunctionShift are not specified, they are assumed
  ! to be =1.
  ! The operator is added to rvector.
  ! The calculation is done based on the degrees of freedom in
  ! rfunction, which corresponds to the operator
  ! nwder_MinMaxProjByMass.
  !
  ! NOTE: This assumes rfunctionMin and rfunctionMax to be
  ! discretised in the same space as rfunction!
!</description>

!<input>
  ! Weight in front of the operator when being added to rmatrix.
  real(DP), intent(in) :: dweight
  
  ! The mass matrix
  type(t_matrixScalar), intent(in) :: rmassMatrix
  
  ! A temporary vector.
  type(t_vectorScalar), intent(inout) :: rvectorTemp
  
  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFct
  
  ! An FE function
  type(t_vectorScalar), intent(in) :: rfunction

  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFctViolate
  
  ! An FE function specifying where the projection should be
  ! applied. Whereever dwFctViolate*rfunctionViolate violates the bounds,
  ! the bounds are assumed.
  type(t_vectorScalar), intent(in) :: rfunctionViolate
  
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

  ! OPTIONAL: Weight for the shift function.
  ! If not present, =0 is assumed.
  real(DP), optional :: dwShift
  
  ! OPTIONAL: Shift function
  type(t_vectorScalar), intent(in), optional :: rfunctionShift
!</input>

!<inputoutput>
  ! The vector which receives the projection.
  type(t_vectorScalar), intent(inout) :: rvector
!<inputoutput>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_vectorBlock), target :: rvecFct, rvecMin, rvecMax
    real(DP), dimension(:), pointer :: p_DdataIn,p_DdataOut,p_DdataMin,p_DdataMax
    real(DP), dimension(:), pointer :: p_DdataShift,p_DdataViol,p_DdataRhs
    integer :: i
    real(DP) :: duShift

    ! Get the source and target arrays.
    call lsyssc_getbase_double (rfunction,p_DdataIn)
    call lsyssc_getbase_double (rvectorTemp,p_DdataOut)
    call lsyssc_getbase_double (rfunctionViolate,p_DdataViol)
    call lsyssc_getbase_double (rvector,p_DdataRhs)
    
    duShift = 0.0_DP
    if (present(dwShift)) duShift = dwShift

    if ((duShift .ne. 0.0_DP) .and. present(rfunctionShift)) then
    
      call lsyssc_getbase_double (rfunctionShift,p_DdataShift)
    
      ! How to compute? Nonconstant min/max available?
      if (present(rfunctionMin) .and. present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector, componentwise
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin*p_DdataMin(i), dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do
        
      else if (present(rfunctionMin)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin*p_DdataMin(i), dwMax,&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do

      else if (present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin, dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do

      else
      
        ! Constant bounds
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin, dwMax,&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do

      end if

    else
    
      ! How to compute? Nonconstant min/max available?
      if (present(rfunctionMin) .and. present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin*p_DdataMin(i), dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do
        
      else if (present(rfunctionMin)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin*p_DdataMin(i), dwMax,&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do

      else if (present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin, dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do

      else
      
        ! Constant bounds
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin, dwMax,&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do

      end if
      
    end if
    
    ! Multiply by the mass matrix and sum up to rvector.
    call lsyssc_scalarMatVec (&
        rmassMatrix,rvectorTemp,rvector,dweight,1.0_DP)
  
  contains
  
    subroutine projection (dtest,dmin,dmax,din,dout)
    
    ! The actual projection.
    ! Returns dmin/dmax where dtest violates the bounds,
    ! or din where it does not.
    
    real(DP), intent(in) :: dtest,dmin,dmax,din
    real(DP), intent(out) :: dout
    
      if (dtest .lt. dmin) then
        dout = dmin
      else if (dtest .gt. dmax) then
        dout = dmax
      else
        dout = din
      end if
    
    end subroutine
  
  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine nwder_applyMinMaxProjByDof (dweight,rvector,&
      dwFctViolate,rfunctionViolate,dwFct,rfunction,dwMin,dwMax,rfunctionMin,rfunctionMax,&
      dwShift,rfunctionShift)
  
!<description>
  ! Calculates 
  !   rvector = dweight * rfunction
  ! for the projection function
  !   rfunction -> 
  !      dwmax*rfunctionMax + dwShift*rfunctionShift
  !            where dwFctViolate*rfunctionViolate + dwShift*rfunctionShift > dwmax*rfunctionMax
  !      dwmin*rfunctionMin + dwShift*rfunctionShift
  !            where dwFctViolate*rfunctionViolate + dwShift*rfunctionShift < dwmin*rfunctionMin
  !      dwFct*rfunction
  !            elsewhere
  ! If rfunctionMin/rfunctionMax are not specified, they are assumed
  ! to be =1.
  ! The operator is added to rvector.
  ! The calculation is done based on the degrees of freedom in
  ! rfunction, which corresponds to the operator
  ! nwder_MinMaxProjByMass.
  !
  ! NOTE: This assumes rfunctionMin and rfunctionMax to be
  ! discretised in the same space as rfunction!
!</description>

!<input>
  ! Weight in front of the operator when being added to rmatrix.
  real(DP), intent(in) :: dweight
  
  ! Weight for the function rfunction.
  real(DP), intent(in) :: dwFctViolate
  
  ! An FE function specifying where the projection should be
  ! applied. Whereever dwFctViolate*rfunctionViolate violates the bounds,
  ! the bounds are assumed.
  type(t_vectorScalar), intent(in) :: rfunctionViolate
  
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

  ! OPTIONAL: Weight for the shift function.
  ! If not present, =0 is assumed.
  real(DP), optional :: dwShift
  
  ! OPTIONAL: Shift function
  type(t_vectorScalar), intent(in), optional :: rfunctionShift
!</input>

!<inputoutput>
  ! The vector which receives the projection.
  type(t_vectorScalar), intent(inout) :: rvector
!<inputoutput>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_vectorBlock), target :: rvecFct, rvecMin, rvecMax,rvecShift
    real(DP), dimension(:), pointer :: p_DdataIn,p_DdataOut,p_DdataMin,p_DdataMax
    real(DP), dimension(:), pointer :: p_DdataShift,p_DdataViol
    integer :: i
    real(DP) :: dushift

    ! Get the source and target arrays.
    call lsyssc_getbase_double (rfunction,p_DdataIn)
    call lsyssc_getbase_double (rvector,p_DdataOut)
    call lsyssc_getbase_double (rfunctionViolate,p_DdataViol)

    duShift = 0.0_DP
    if (present(dwShift)) duShift = dwShift

    if ((duShift .ne. 0.0_DP) .and. present(rfunctionShift)) then
    
      call lsyssc_getbase_double (rfunctionShift,p_DdataShift)
    
      ! How to compute? Nonconstant min/max available?
      if (present(rfunctionMin) .and. present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector, componentwise
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin*p_DdataMin(i), dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do
        
      else if (present(rfunctionMin)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin*p_DdataMin(i), dwMax,&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do

      else if (present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin, dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do

      else
      
        ! Constant bounds
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift*p_DdataShift(i),&
              dwMin, dwMax,&
              dwFct*p_DdataIn(i) + duShift*p_DdataShift(i),&
              p_DdataOut(i))
        end do

      end if

    else
    
      ! How to compute? Nonconstant min/max available?
      if (present(rfunctionMin) .and. present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin*p_DdataMin(i), dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do
        
      else if (present(rfunctionMin)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMin,p_DdataMin)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin*p_DdataMin(i), dwMax,&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do

      else if (present(rfunctionMax)) then
      
        ! Get the bounding functions
        call lsyssc_getbase_double (rfunctionMax,p_DdataMax)
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin, dwMax*p_DdataMax(i),&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do

      else
      
        ! Constant bounds
      
        ! Restrict the vector
        do i=1,rvector%NEQ
          call projection (dwFctViolate*p_DdataViol(i) + duShift,&
              dwMin, dwMax,&
              dwFct*p_DdataIn(i) + duShift,&
              p_DdataOut(i))
        end do

      end if
      
    end if
    
  contains
    
    subroutine projection (dtest,dmin,dmax,din,dout)
    
    ! The actual projection.
    ! Returns dmin/dmax where dtest violates the bounds,
    ! or din where it does not.
    
    real(DP), intent(in) :: dtest,dmin,dmax,din
    real(DP), intent(out) :: dout
    
      if (dtest .lt. dmin) then
        dout = dmin
      else if (dtest .gt. dmax) then
        dout = dmax
      else
        dout = din
      end if
    
    end subroutine

  end subroutine

end module

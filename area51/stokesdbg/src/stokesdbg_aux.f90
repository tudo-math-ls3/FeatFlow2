module stokesdbg_aux

use fsystem
use storage
use genoutput
use cubature
use derivatives
use boundary
use triangulation
use domainintegration
use element
use spatialdiscretisation
use linearsystemscalar
use scalarpde
use collection, only: t_collection
use dofmapping
use transformation
use elementpreprocessing
use feevaluation
use feevaluation2
use blockmatassembly
use blockmatassemblybase
use blockmatassemblystdop

implicit none

contains

  ! ***********************************************************************************************

  subroutine stdbg_aux_funcZeroBC2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                     cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  integer, dimension(:), intent(in)                           :: Icomponents
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  integer, intent(in)                                         :: ielement
  integer, intent(in)                                         :: cinfoNeeded
  integer, intent(in)                                         :: iwhere
  real(DP), intent(in)                                        :: dwhere
  type(t_collection), intent(inout), optional                 :: rcollection
  real(DP), dimension(:), intent(out)                         :: Dvalues

    Dvalues = 0.0_DP

  end subroutine

  ! ***********************************************************************************************

  subroutine stdbg_aux_funcParProfileBC2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                     cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  integer, dimension(:), intent(in)                           :: Icomponents
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  integer, intent(in)                                         :: ielement
  integer, intent(in)                                         :: cinfoNeeded
  integer, intent(in)                                         :: iwhere
  real(DP), intent(in)                                        :: dwhere
  type(t_collection), intent(inout), optional                 :: rcollection
  real(DP), dimension(:), intent(out)                         :: Dvalues

  real(DP) :: y
  
    ! clear the output array; for vector-valued problems, this will set the Y-velocity to zero
    Dvalues = 0.0_DP
    if ((dwhere .ge. 3.0_DP) .and. (dwhere .le. 4.0_DP)) then
      y = 4.0_DP-dwhere
      Dvalues(1) = y*(1.0_DP-y)
    end if
  end subroutine


  ! ***********************************************************************************************

  subroutine stdbg_aux_funcRhsOne (rdiscretisation,rform,nelements,npointsPerElement,Dpoints, &
                         IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_linearForm), intent(in)                              :: rform
  integer, intent(in)                                         :: nelements
  integer, intent(in)                                         :: npointsPerElement
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints
  integer, dimension(:,:), intent(in) :: IdofsTest
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset
  type(t_collection), intent(inout), optional      :: rcollection
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  
    Dcoefficients = 1.0_DP
  
  end subroutine

  !****************************************************************************

  subroutine stdbg_aux_calcDiv2D (derror, rvector, rcubInfo)

!<description>
  ! This routine calculates the divergence of a 2D vector field.
!</description>

!<input> 
  ! The FE solution vector.
  type(t_vectorScalar), dimension(:), intent(in), target :: Rvector
  
  type(t_scalarCubatureInfo), intent(in) :: rcubInfo

!<output>
  ! The L2-norm of the divergence.
  real(DP), intent(out) :: derror
!</output>

!</subroutine>

  type(t_fev2Vectors) :: rvectorList

    derror = 0.0_DP
    call fev2_addVectorToEvalList(rvectorList,rvector(1),1)
    call fev2_addVectorToEvalList(rvectorList,rvector(2),1)
    call bma_buildIntegral (derror,BMA_CALC_STANDARD, bma_fcalc_divergenceL2norm, &
        revalVectors=rvectorList, rcubatureInfo=rcubInfo)
    call fev2_releaseVectorList(rvectorList)
    derror = sqrt(derror)

  end subroutine

!  old divergence computation routine follows
!
!  !****************************************************************************
!
!  subroutine stdbg_aux_calcDiv2D2 (derror, rvector)
!
!!<description>
!  ! This routine calculates the divergence of a 2D vector field.
!!</description>
!
!!<input> 
!  ! The FE solution vector.
!  type(t_vectorScalar), dimension(:), intent(in), target :: Rvector
!
!!<output>
!  ! The L2-norm of the divergence.
!  real(DP), intent(out) :: derror
!!</output>
!
!!</subroutine>
!
!  ! A discretisation structure specifying how to compute the error.
!  type(t_spatialDiscretisation), pointer :: p_rdiscretisation
!
!    ! local variables
!    integer :: ielementDistr, ICUBP, NVE, NCOEFF
!    integer :: IEL, IELmax, IELset, IELGlobal
!    real(DP) :: OM,dt
!    
!    ! Array to tell the element which derivatives to calculate
!    logical, dimension(EL_MAXNDER) :: Bder
!    
!    ! For every cubature point on the reference element,
!    ! the corresponding cubature weight
!    real(DP), dimension(:), allocatable :: Domega
!    
!    ! number of cubature points on the reference element
!    integer :: ncubp
!    
!    ! Number of local degees of freedom for test functions
!    integer :: indofTrial
!    
!    ! The triangulation structure - to shorten some things...
!    type(t_triangulation), pointer :: p_rtriangulation
!    
!    ! A pointer to an element-number list
!    integer, dimension(:), pointer :: p_IelementList
!    
!    ! An array receiving the coordinates of cubature points on
!    ! the reference element for all elements in a set.
!    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef
!
!    ! Arrays for saving Jacobian determinants and matrices
!    real(DP), dimension(:,:), pointer :: p_Ddetj
!    
!    ! Current element distribution
!    type(t_elementDistribution), pointer :: p_relementDistribution
!    
!    ! Number of elements in the current element distribution
!    integer :: NEL
!
!    ! Pointer to the values of the function that are computed by the callback routine.
!    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
!    
!    ! Number of elements in a block. Normally =PPERR_NELEMSIM,
!    ! except if there are less elements in the discretisation.
!    integer :: nelementsPerBlock
!    
!    ! A t_domainIntSubset structure that is used for storing information
!    ! and passing it to callback routines.
!    type(t_domainIntSubset) :: rintSubset
!    type(t_evalElementSet) :: revalElementSet
!    
!    ! An allocateable array accepting the DOF`s of a set of elements.
!    integer, dimension(:,:), allocatable, target :: IdofsTrial
!  
!    ! Type of transformation from the reference to the real element 
!    integer(I32) :: ctrafoType, ccubature
!    
!    ! Element evaluation tag; collects some information necessary for evaluating
!    ! the elements.
!    integer(I32) :: cevaluationTag
!
!    ! Pointer to the element-wise error
!    real(DP), dimension(:), pointer :: p_Derror
!
!
!    ! Which derivatives of basis functions are needed?
!    ! Check the descriptors of the bilinear form and set BDER
!    ! according to these.
!
!    Bder = .false.
!    Bder(DER_DERIV_X) = .true.
!    Bder(DER_DERIV_Y) = .true.
!
!    ! Fetch the discretisation
!    p_rdiscretisation => Rvector(1)%p_rspatialDiscr
!
!    ! Get a pointer to the triangulation - for easier access.
!    p_rtriangulation => p_rdiscretisation%p_rtriangulation
!    
!    ! For saving some memory in smaller discretisations, we calculate
!    ! the number of elements per block. For smaller triangulations,
!    ! this is NEL. If there are too many elements, it is at most
!    ! PPERR_NELEMSIM. This is only used for allocating some arrays.
!    nelementsPerBlock = min(1000, p_rtriangulation%NEL)
!    
!    ! Set the current error to 0 and add the error contributions of each element
!    ! to that.
!    Derror = 0.0_DP
!    
!    ! Now loop over the different element distributions (=combinations
!    ! of trial and test functions) in the discretisation.
!
!    do ielementDistr = 1, p_rdiscretisation%inumFESpaces
!    
!      ! Activate the current element distribution
!      p_relementDistribution => p_rdiscretisation%RelementDistr(ielementDistr)
!    
!      ! Cancel if this element distribution is empty.
!      if (p_relementDistribution%NEL .eq. 0) cycle
!
!      ! Get the number of local DOF`s for trial functions
!      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
!      
!      ! Get the number of corner vertices of the element
!      NVE = elem_igetNVE(p_relementDistribution%celement)
!      
!      ! Get from the trial element space the type of coordinate system
!      ! that is used there:
!      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
!
!      ccubature = p_relementDistribution%ccubTypeEval
!
!      ! Get the number of cubature points for the cubature formula
!      ncubp = cub_igetNumPts(ccubature)
!      
!      ! Allocate two arrays for the points and the weights
!      allocate(Domega(ncubp))
!      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
!      
!      ! Get the cubature formula
!      call cub_getCubature(ccubature, p_DcubPtsRef, Domega)
!      
!      ! Allocate memory for the DOF`s of all the elements.
!      allocate(IdofsTrial(indofTrial, nelementsPerBlock))
!
!      ! Allocate memory for the coefficients
!      allocate(Dcoefficients(ncubp, nelementsPerBlock, 2))
!    
!      ! Initialisation of the element set.
!      call elprep_init(revalElementSet)
!
!      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
!      ! the elements later. All of them can be combined with OR, what will give
!      ! a combined evaluation tag. 
!      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
!
!      ! Make sure that we have determinants.
!      cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)
!
!      ! p_IelementList must point to our set of elements in the discretisation
!      ! with that combination of trial functions
!      call storage_getbase_int (p_relementDistribution%h_IelementList, &
!                                p_IelementList)
!                     
!      ! Get the number of elements there.
!      NEL = p_relementDistribution%NEL
!    
!      ! Loop over the elements - blockwise.
!      do IELset = 1, NEL, 1000
!      
!        ! We always handle LINF_NELEMSIM elements simultaneously.
!        ! How many elements have we actually here?
!        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
!        ! elements simultaneously.
!        
!        IELmax = min(NEL,IELset-1+1000)
!      
!        ! Calculate the global DOF`s into IdofsTrial.
!        !
!        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
!        ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
!        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
!                                     IdofsTrial)
!                                     
!        ! Prepare the call to the evaluation routine of the analytic function.    
!        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
!        rintSubset%ielementStartIdx = IELset
!        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
!        rintSubset%p_IdofsTrial => IdofsTrial
!        rintSubset%celement = p_relementDistribution%celement
!    
!        ! Calculate all information that is necessary to evaluate the finite element
!        ! on all cells of our subset. This includes the coordinates of the points
!        ! on the cells.
!        call elprep_prepareSetForEvaluation (revalElementSet,&
!            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
!            ctrafoType, p_DcubPtsRef(:,1:ncubp))
!        p_Ddetj => revalElementSet%p_Ddetj
!
!        ! In the next loop, we do not have to evaluate the coordinates
!        ! on the reference elements anymore.
!        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
!
!        ! Calculate the X/Y-derivative of the FE function in the
!        ! cubature points: u_h(x,y).
!        call fevl_evaluate_sim3 (Rvector(1), revalElementSet,&
!                p_relementDistribution%celement, IdofsTrial, DER_DERIV_X,&
!                Dcoefficients(:,1:IELmax-IELset+1,1))
!
!        call fevl_evaluate_sim3 (Rvector(2), revalElementSet,&
!                p_relementDistribution%celement, IdofsTrial, DER_DERIV_Y,&
!                Dcoefficients(:,1:IELmax-IELset+1,2))
!
!        do IEL=1,IELmax-IELset+1
!              
!          ! Loop over all cubature points on the current element
!          do icubp = 1, ncubp
!                
!            ! calculate the current weighting factor in the cubature formula
!            ! in that cubature point.
!            !
!            ! Take the absolut value of the determinant of the mapping.
!            ! In 2D, the determinant is always positive, whereas in 3D,
!            ! the determinant might be negative -- that is normal!
!                
!            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
!
!            ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
!            derror = derror + OM * (Dcoefficients(icubp,IEL,1)+Dcoefficients(icubp,IEL,2))**2
!
!          end do ! ICUBP 
!              
!        end do ! IEL
!
!        ! Release the temporary domain integration structure again
!        call domint_doneIntegration (rintSubset)
!    
!      end do ! IELset
!      
!      ! Release memory
!      call elprep_releaseElementSet(revalElementSet)
!
!      deallocate(p_DcubPtsRef)
!      deallocate(Dcoefficients)
!      deallocate(IdofsTrial)
!      deallocate(Domega)
!
!    end do ! ielementDistr
!
!    ! derror is ||error||^2, so take the square root at last.
!    derror = sqrt(derror)
!
!  end subroutine 

end module

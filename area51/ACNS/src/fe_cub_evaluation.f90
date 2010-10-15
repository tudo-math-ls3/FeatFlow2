module fe_cub_evaluation

  use fsystem
  use storage
  use genoutput
  use cubature
  use basicgeometry
  use element
  use derivatives
  use transformation
  use triangulation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use scalarpde
  use spatialdiscretisation
  use domainintegration
  use elementpreprocessing
  use feevaluation
  use collection
  use dofmapping

  IMPLICIT NONE

!<constantblock description = "Identifiers for the type of error to be computed.">

  ! $L_2$-error/norm
  integer, parameter :: PPERR_L2ERROR = 1
  
  ! $H_1$-error/norm
  integer, parameter :: PPERR_H1ERROR = 2

  ! $L_1$-error/norm
  integer, parameter :: PPERR_L1ERROR = 3
  
  
!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! Number of elements to handle simultaneously when building vectors
  integer :: PPERR_NELEMSIM   = 1000
contains


!*******************************************************************

  subroutine fun_integral2d_conf (rvectorScalar,cerrortype,derror,&
                                  rdiscretisation,ffunctionReference,rcollection)


!<description>
  ! This routine calculates the integration of a given finite element function
  ! in rvectorscalar to a given analytical callback function ffunctionReference.
  ! ffunctionReference will provided the function values at the cubature points. 
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  
  ! type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN)                      :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information for callback routines.
  type(t_collection), intent(INOUT), optional      :: rcollection

  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFuncSc.inc'
  optional :: ffunctionReference
!</input>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer(I32) :: IEL, IELmax, IELset
    real(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
  
    ! type of transformation from the reference to the real element 
    integer(I32) :: ctrafotype
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

! Mcai, we first scale rvectorScalar to be 0
!    call lsyssc_scaleVector(rvectorScalar, 0.0_DP)

    Bder = .false.
    select case (cerrortype)
    case (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC) = .true.
    case (PPERR_H1ERROR) 
      Bder(DER_DERIV_X) = .true.
      Bder(DER_DERIV_Y) = .true.
    case DEFAULT
      call output_line('Unknown error type identifier!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
      call sys_halt()
    end select
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistribution%ccubtypeEval, ncubp, Dxi, Domega)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafotype = elem_igetTrafotype(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafotype),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp,nelementsPerBlock,4))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      if (present(ffunctionReference)) then
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_realPOINTS)
      end if
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafotype, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        select case (cerrortype)
        
        case (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC,rdiscretisation, &
                        int(IELmax-IELset+1),ncubp,&
                        revalElementSet%p_Dpointsreal,&
                        IdofsTrial,rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1_I32,1),rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1_I32,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1_I32,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx
          
          do IEL=1,IELmax-IELset+1
          
            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp
            
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.

              OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
              
              ! L1-error is:   int_... abs(u-u_h) dx
              
              derror = derror + &
                       OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))

            end do ! ICUBP 

          end do ! IEL

        case DEFAULT
          call output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
          call sys_halt()
        end select
        
        ! Release the temporary domain integration structure again
        call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)

    end do ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    if ((cerrortype .eq. PPERR_L2ERROR) .or.&
        (cerrortype .eq. PPERR_H1ERROR)) then
      derror = sqrt(derror)
    end if

  end subroutine
!*****************************************************************************

 ! ***************************************************************************
! The following subroutine is to calculate  
! \sum_{E} {\sum{\Phi_i b_i}^3-\sum{\Phi_i b_i}}
! Question: is it same as 
! (\sum_{E} {\sum{\Phi_i b_i}})^3-(\sum_{E}[\sum{\Phi_i b_i}])?

!<subroutine>

  subroutine fevl_evaluate_cubicpoly (rvectorScalar, &
                            rdomainIntSubset, idertype, Dvalues, iterm)
                                      
!<description>
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on multiple elements in a 
  ! discretisation.
  ! rdomainIntSubset must specify all information about where and how
  ! to evaluate; e.g. the coordinates of the evaluation points are
  ! to be found here.
  !
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating 'simultaneously' on some $Q_1$ and some $P_1$ elements
  ! is not allowed e.g.. 
  !
  ! The interface of this routine is designed to be called in callback
  ! functions during linearform and bilinearform evaluation.
  ! The target array Dvalues provides a shape which is compatible
  ! to the callback interface. The variable iterm specifies the
  ! subarray in Dvalues where values are written to; i.e.
  ! the result of the evaluation is written to Dvalues(iterm,:,:).
!</description>

!<input>
  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being evaluated.
  type(t_domainIntSubset), intent(IN)            :: rdomainIntSubset

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)               :: rvectorScalar
  
  ! type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: idertype
  
  ! Number of the subarray in Dvalues where the result is written to.
  ! The routine writes the resulting values to Dvalues(iterm,:,:).
  integer, intent(in)                            :: iterm
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:,:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  real(DP) :: dval
  integer :: iel,ipoint,ibas
  integer(I32) :: ieltyp
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  integer, dimension(:,:), pointer :: p_IdofsTrial
  
  npoints = ubound(Dvalues,2)
  nelements = ubound(Dvalues,3)
  
  ! What to evaluate?
  Bder = .false.
  Bder(idertype) = .true.
  
  ! Get the pointer to the trail DOF's.
  p_IdofsTrial => rdomainIntSubset%p_IdofsTrial
  
  ! Get the currently active element
  ieltyp = rvectorScalar%p_rspatialDiscr%RelementDistr( &
      rdomainIntSubset%ielementDistribution)%celement
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(ieltyp)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(ieltyp),npoints,nelements))
  
  ! Evaluate the basis functions
  call elem_generic_sim2 (ieltyp, rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
  
  if (rvectorScalar%cdatatype .eq. ST_DOUBLE) then
  
    ! Get the data array from the vector
    call lsyssc_getbase_double(rvectorScalar,p_Ddata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    do iel=1,nelements
      do ipoint = 1,npoints
        ! Calculate the value in the point
        dval = 0.0_DP
        do ibas = 1,indofTrial
          dval = dval + &
                 p_Ddata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,idertype,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(iterm,ipoint,iel) = (dval**2-1.0_DP)*dval
      end do

    end do
    
  else if (rvectorScalar%cdatatype .eq. ST_SINGLE) then
  
    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    do iel=1,nelements
      do ipoint = 1,npoints
        ! Calculate the value in the point
        dval = 0.0_DP
        do ibas = 1,indofTrial
          dval = dval + &
                 p_Fdata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,idertype,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(iterm,ipoint,iel) = (dval**2-1.0_DP)*dval
      end do
    end do
    
  else
    call output_line('Unsupported vector precision!',&
      OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate_sim')
    call sys_halt()
  end if
  
  ! Release memory, finish
  deallocate(DbasTrial)

  end subroutine

!**********************************************************************************************


end module

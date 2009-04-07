!##############################################################################
!# ****************************************************************************
!# <name> feevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to evaluate finite element functions
!# in given points.
!#
!# The module contains the following routines:
!#
!# 1.) fevl_evaluate
!#     -> Evaluate a FE function in an arbitrary set of points.
!#        Determines automatically the elements that contain the points.
!#
!# 3.) fevl_evaluate_mult
!#     -> Evaluate a FE function simultaneously in multiple points on
!#        one element. The element containing the point must be given.
!#
!# 3.) fevl_evaluate_sim
!#     -> Evaluate a FE function simultaneously in multiple points on
!#        multiple elements. Elements containing the points must be given.
!#
!# 4.) fevl_getVectorMagnitude
!#     -> Calculate the maximum norm of a given FEM vector field.
!# </purpose>
!##############################################################################

module feevaluation

  use fsystem
  use linearsystemscalar
  use linearsystemblock
  use triasearch
  use element
  use elementpreprocessing
  use domainintegration
  use derivatives
  use spdiscprojection
  
  implicit none

!<constants>

!<constantblock description="Constants for the cnonmeshPoints parameter.">

  ! Points outside of the domain not allowed.
  integer, parameter :: FEVL_NONMESHPTS_NONE   = 0

  ! Points outside of the domain allowed, try to find an element nearby to evaluate.
  integer, parameter :: FEVL_NONMESHPTS_NEARBY = 1

  ! Points outside of the domain allowed, assume 0 as value there.
  integer, parameter :: FEVL_NONMESHPTS_ZERO   = 2

!</constantblock>

!</constants>

  ! There are two functions fevl_evaluate_mult which do the same -- with
  ! different calling conventions and different complexities.
  interface fevl_evaluate_mult
    module procedure fevl_evaluate_mult1
    module procedure fevl_evaluate_mult2
  end interface

  ! There are two functions fevl_evaluate_sim which do the same -- with
  ! different calling conventions and different complexities.
  interface fevl_evaluate_sim
    module procedure fevl_evaluate_sim1
    module procedure fevl_evaluate_sim2
    module procedure fevl_evaluate_sim3
    module procedure fevl_evaluate_sim4
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate (iderType, Dvalues, rvectorScalar, Dpoints, &
      Ielements, IelementsHint, cnonmeshPoints)
                                      
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)              :: rvectorScalar
  
  ! A list of points where to evaluate.
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(:,:), intent(IN) :: Dpoints
  
  ! OPTIONAL: A list of elements containing the points Dpoints.
  ! If this is not specified, the element numbers containing the points
  ! are determined automatically.
  integer, dimension(:), intent(IN), optional :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements
  ! containing the points. This is ignored if Ielements is specified!
  integer, dimension(:), intent(IN), optional :: IelementsHint
  
  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
  ! if a point is located outside of the domain. May happen e.g. in
  ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default 
  ! parameter if cnonmeshPoints is not specified. 
  integer, intent(IN), optional :: cnonmeshPoints
  
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(DP), dimension(:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: cnonmesh
    integer :: ipoint,indof,nve,ibas
    integer(I32) :: celement
    integer :: iel
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval
    
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    
    ! Transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint
    
    ! Values of basis functions and DOF's
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    integer, dimension(EL_MAXNBAS) :: Idofs
    
    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! Evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag

    ! Ok, slow but general.
    
    p_RelementDistribution => rvectorScalar%p_rspatialDiscr%RelementDistr
    
    ! For uniform discretisations, we get the element type in advance...
    if (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      
      ! Element type
      celement = rvectorScalar%p_rspatialDiscr%&
          RelementDistr(1)%celement

      ! Get the number of local DOF's for trial and test functions
      indof = elem_igetNDofLoc(celement)
      
      ! Number of vertices on the element
      nve = elem_igetNVE(celement)
      
      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)
      
      ! Get the element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)
      
      nullify(p_IelementDistr)
    else
      call storage_getbase_int (&
          rvectorScalar%p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
    end if
    
    ! Get the data vector
    select case (rvectorScalar%cdataType)
    case (ST_DOUBLE) 
      call lsyssc_getbase_double(rvectorScalar,p_Ddata)
    case (ST_SINGLE)
      call lsyssc_getbase_single(rvectorScalar,p_Fdata)
    case DEFAULT
      call output_line ('Unsupported vector precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
      call sys_halt()
    end select
    
    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.
    
    cnonmesh = FEVL_NONMESHPTS_NONE
    if (present(cnonmeshPoints)) cnonmesh = cnonmeshPoints
    
    ! We loop over all points.

    iel = 1

    do ipoint = 1,ubound(Dpoints,2)
    
      ! Get the element number that contains the point. 
      if (present(Ielements)) then
        ! Either we have that...
        iel = Ielements (ipoint)
      else
        ! Or we have to find the element. Probably, the caller gave us a
        ! hint where we can start the search...
        if (present(IelementsHint)) then
          iel = IelementsHint (ipoint)
        end if
        
        ! Otherwise, we use iel from the previous iteration as starting point
        ! to search for the element. Use raytracing search to find the element
        ! containing the point.
        call tsrch_getElem_raytrace2D (&
          Dpoints(:,ipoint),rvectorScalar%p_rspatialDiscr%p_rtriangulation,iel)

        ! Ok, not found... Brute force search
        if (iel .eq. 0) then
          call tsrch_getElem_BruteForce (Dpoints(:,ipoint), &
            rvectorScalar%p_rspatialDiscr%p_rtriangulation,iel)
        end if
        
        if (iel .eq. 0) then
        
          ! We really have not luck here... Are nonmesh-points allowed?
          if (cnonmesh .eq. FEVL_NONMESHPTS_NEARBY) then
          
            ! Yes. Find the closest element!
            call tsrch_getNearestElem_BruteForce (Dpoints(:,ipoint), &
              rvectorScalar%p_rspatialDiscr%p_rtriangulation,iel)
              
            ! The evaluation routine then evaluates the FE function
            ! outside of the element...
            
          else if (cnonmesh .eq. FEVL_NONMESHPTS_ZERO) then
            
            ! Assume zero here.
            Dvalues(ipoint) = 0.0_DP
            
            cycle
            
          end if
          
        end if

      end if

      if (iel .eq. 0) then
        call output_line ('Point '//trim(sys_siL(ipoint,10))//' not found!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
        cycle
      end if
    
      ! Get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! Get the number of local DOF's for trial and test functions
        indof = elem_igetNDofLoc(celement)
        
        ! Number of vertices on the element
        nve = elem_igetNVE(celement)
        
        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)
        
        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)
      end if
        
      ! Calculate the global DOF's on that element into IdofsTest.
      call dof_locGlobMapping (rvectorScalar%p_rspatialDiscr, &
          iel,Idofs)
     
      ! Get the element shape information
      call elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
          rvectorScalar%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)
          
      ! Calculate the transformation of the point to the reference element
      call trafo_calcRefCoords (ctrafoType,revalElement%Dcoords,&
          Dpoints(:,ipoint),DparPoint)

      ! Now calculate everything else what is necessary for the element    
      call elprep_prepareForEvaluation (revalElement, &
          iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
          rvectorScalar%p_rspatialDiscr%p_rtriangulation, iel, &
          ctrafoType, DparPoint, Dpoints(:,ipoint))

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)
      
      ! Combine the basis functions to get the function value.
      dval = 0.0_DP
      if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !          
        ! Calculate the value in the point
        do ibas = 1,indof
          dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
        end do
      
      else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !
        ! Calculate the value in the point
        do ibas = 1,indof
          dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
        end do
        
      end if

      ! Save the value in the point
      Dvalues(ipoint) = dval
      
    end do ! ipoint

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_mult1 (iderType, Dvalues, rvectorScalar, ielement, &
      DpointsRef, Dpoints)
                                      
!<description>
  ! This is a rather general finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints on one element ielement. 
  ! The values of the FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)              :: rvectorScalar
  
  ! The element number containing the points in Dpoints.
  integer, intent(IN) :: ielement
  
  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! Either Dpoints or DpointsRef must be specified!
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(:,:), intent(IN), optional :: DpointsRef
  
  ! OPTIONAL: A list of points where to evaluate. All points must be inside
  ! of element ielement. 
  ! If not specified, the coordinates are automatically calculated.
  ! Either Dpoints or DpointsRef must be specified!
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(:,:), intent(IN), optional :: Dpoints
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(DP), dimension(:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: ipoint,indof,nve,ibas,npoints
    integer(I32) :: celement
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval
    
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    
    ! Transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint
    
    ! Values of basis functions and DOF's
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    integer, dimension(EL_MAXNBAS) :: Idofs
    
    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! Evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag

    ! Ok, slow but general.
    
    ! Are points given?
    if ((.not. present(Dpoints)) .and. (.not. present(DpointsRef))) then
      call output_line ('Evaluation points not specified!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate_mult')  
    end if
    
    if (present(DpointsRef)) then
      npoints = ubound(DpointsRef,2)
    else
      npoints = ubound(Dpoints,2)
    end if
    
    p_RelementDistribution => rvectorScalar%p_rspatialDiscr%RelementDistr
    
    ! For uniform discretisations, we get the element type in advance...
    if (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      
      ! Element type
      celement = p_RelementDistribution(1)%celement

    else
      call storage_getbase_int (rvectorScalar%p_rspatialDiscr%h_IelementDistr,&
          p_IelementDistr)

      celement = p_RelementDistribution(p_IelementDistr(ielement))%celement
    end if

    ! Get the number of local DOF's for trial and test functions
    indof = elem_igetNDofLoc(celement)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(celement)
    
    ! Type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(celement)
    
    ! Get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(celement)

    ! Get the data vector
    select case (rvectorScalar%cdataType)
    case (ST_DOUBLE) 
      call lsyssc_getbase_double(rvectorScalar,p_Ddata)
    case (ST_SINGLE)
      call lsyssc_getbase_single(rvectorScalar,p_Fdata)
    case DEFAULT
      call output_line ('Unsupported vector precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
      call sys_halt()
    end select
    
    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.
    
    ! Calculate the global DOF's on that element into IdofsTest.
    call dof_locGlobMapping (rvectorScalar%p_rspatialDiscr, &
        ielement,Idofs)
        
    ! Get the element shape information
    call elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
        rvectorScalar%p_rspatialDiscr%p_rtriangulation, ielement, ctrafoType)
  
    ! We loop over all points
    do ipoint = 1,npoints
    
      if (present(DpointsRef)) then
        DparPoint(1:ubound(DpointsRef,1)) = DpointsRef(:,ipoint)
      else
        ! Calculate the transformation of the point to the reference element
        call trafo_calcRefCoords (ctrafoType,revalElement%Dcoords,Dpoints(:,ipoint),DparPoint)
      end if

      ! Now calculate everything else what is necessary for the element.
      ! Don't calculate the shape of the cell again since we did this in advance above.
      if (present(Dpoints)) then
        call elprep_prepareForEvaluation (revalElement, &
            iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
            rvectorScalar%p_rspatialDiscr%p_rtriangulation, ielement, &
            ctrafoType, DparPoint, Dpoints(:,ipoint))
      else
        call elprep_prepareForEvaluation (revalElement, &
            iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
            rvectorScalar%p_rspatialDiscr%p_rtriangulation, ielement, &
            ctrafoType, DparPoint)
      end if

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)
      
      dval = 0.0_DP
      if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !          
        ! Calculate the value in the point
        do ibas = 1,indof
          dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
        end do
      
      else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !
        ! Calculate the value in the point
        do ibas = 1,indof
          dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
        end do
        
      end if

      ! Save the value in the point
      Dvalues(ipoint) = dval
      
    end do ! ipoint

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_mult2 (rvectorScalar, Dcoords, Djac, Ddetj, &
                  celement, IdofsTrial, npoints,  Dpoints, iderType,&
                  Dvalues, ItwistIndex)
                                      
!<description>
  ! DEPRECATED!
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on one elements in a 
  ! discretisation.
  ! The caller must provide all necessary information for the evaluation in the
  ! parameters.
!</description>

!<input>
  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)               :: rvectorScalar
  
  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those 
  ! elements that we are concerning here.
  integer(I32), intent(IN)                       :: celement

  ! A list of the corner vertices of the element.
  ! array [1..NDIM2D,1..TRIA_MAXNVE2D] of double
  real(DP), dimension(:,:), intent(IN)           :: Dcoords
  
  ! The Jacobian matrix of the mapping between the reference and the
  ! real element, for all points on the element.
  ! array [1..TRAFO_NJACENTRIES,1..npointsPerElement]
  real(DP), dimension(:,:),intent(IN)            :: Djac
  
  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to the real element.
  ! array [1..npointsPerElement]
  real(DP), dimension(:), intent(IN)             :: Ddetj
  
  ! An array accepting the DOF's on the element in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF's in trial space)
  integer, dimension(:), intent(IN)              :: IdofsTrial
  
  ! Number of points on the element where to evalate the function
  integer, intent(IN) :: npoints
  
  ! Array with coordinates of the points where to evaluate.
  ! DIMENSION(NDIM2D,npoints).
  ! The coordinates are expected 
  ! - on the reference element, if celement identifies a parametric element
  ! - on the real element, if celement identifies a nonparametric element
  ! It's assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:), intent(IN)           :: Dpoints

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType

  ! OPTIONAL: Twist index bitfield of the element. Defines for every edge its
  ! orientation.
  ! Can be omitted if the element does not need this information.
  integer(I32), intent(IN), optional             :: itwistIndex

!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints).
  real(DP), dimension(:), intent(OUT)            :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:), allocatable :: DbasTrial
  integer :: indofTrial
  real(DP) :: dval
  integer :: ipoint,ibas
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  
  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints))
  
  ! Evaluate the basis functions
  call elem_generic_mult (celement, Dcoords, Djac, Ddetj, &
                         Bder, DbasTrial, npoints, Dpoints, itwistIndex)  
  
  if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
  
    ! Get the data array from the vector
    call lsyssc_getbase_double(rvectorScalar,p_Ddata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    do ipoint = 1,npoints
      ! Calculate the value in the point
      dval = 0.0_DP
      do ibas = 1,indofTrial
        dval = dval + &
            p_Ddata(IdofsTrial(ibas)) * DbasTrial(ibas,iderType,ipoint)
      end do
      ! Save the value in the point
      Dvalues(ipoint) = dval
    end do
    
  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
  
    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    do ipoint = 1,npoints
      ! Calculate the value in the point
      dval = 0.0_DP
      do ibas = 1,indofTrial
        dval = dval + &
            p_Fdata(IdofsTrial(ibas)) * DbasTrial(ibas,iderType,ipoint)
      end do
      ! Save the value in the point
      Dvalues(ipoint) = dval
    end do
    
  else
    call output_line('Unsupported vector precision!',&
      OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate_mult2')
    call sys_halt()
  end if
  
  ! Release memory, finish
  deallocate(DbasTrial)
  
  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE fevl_evaluate_sim (iderType, Dvalues, rvectorScalar, Dpoints, &
!      Ielements, DpointsRef)
!                                      
!!<description>
!  ! This is a rather general finite element evaluation
!  ! routine. It allows to evaluate a general (scalar) FE function specified
!  ! by rvectorScalar in a set of points Dpoints on a set of elements Ielements. 
!  ! The values of the FE function are written into Dvalues.
!!</description>
!
!!<input>
!  ! Type of function value to evaluate. One of the DER_xxxx constants,
!  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
!  INTEGER, INTENT(IN)                            :: iderType
!
!  ! The scalar solution vector that is to be evaluated.
!  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
!  
!  ! A list of points where to evaluate. All points must be inside
!  ! of element ielement.
!  ! DIMENSION(1..ndim,1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!  
!  ! A list of elements containing the points in Dpoints.
!  ! All elements in this list must be of the same type!!!
!  INTEGER, DIMENSION(:), INTENT(IN) :: Ielements
!  
!  ! OPTIONAL: Coordinates of the points on the reference element.
!  ! If not specified, the coordinates are automatically calculated.
!  ! DIMENSION(1..ndim,1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:,:), INTENT(IN), TARGET, OPTIONAL :: DpointsRef
!!</input>
!
!!<output>
!  ! Values of the FE function at the points specified by Dpoints.
!  ! DIMENSION(1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dvalues
!!</output>
!
!!</subroutine>
!
!    ! local variables
!    LOGICAL :: bnonpar
!    INTEGER :: ipoint,celement,indof,nve,ibas,iel,ndim,ntwistsize
!    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementDistr
!    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
!    REAL(DP) :: dval
!    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DpointsRef
!    
!    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
!    
!    ! Triangulation information
!    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!    INTEGER, DIMENSION(:,:), POINTER :: p_IverticesAtElement
!    
!    ! Transformation
!    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE :: Djac
!    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Ddetj
!    INTEGER(I32) :: ctrafoType
!    
!    ! Values of basis functions and DOF's
!    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Dbas
!    INTEGER, DIMENSION(:,:), ALLOCATABLE :: Idofs
!    
!    ! Coordinates of the corners of one element
!    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE :: Dcoord
!
!    ! List of element distributions in the discretisation structure
!    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: p_RelementDistribution
!
!    ! Twist index array to define the orientation of faces/edges
!    INTEGER(I32), DIMENSION(:,:), ALLOCATABLE :: ItwistIndex
!
!    ! Ok, slow but general.
!    
!    ! Get triangulation information
!    CALL storage_getbase_double2d (&
!        rvectorScalar%p_rspatialDiscr%p_rtriangulation%h_DvertexCoords,&
!        p_DvertexCoords)
!    CALL storage_getbase_int2d (&
!        rvectorScalar%p_rspatialDiscr%p_rtriangulation%h_IverticesAtElement,&
!        p_IverticesAtElement)
!    
!    p_RelementDistribution => rvectorScalar%p_rspatialDiscr%RelementDistr
!
!    ! For uniform discretisations, we get the element type in advance...
!    IF (rvectorScalar%p_rspatialDiscr%ccomplexity .EQ. SPDISC_UNIFORM) THEN
!      
!      ! Element type
!      celement = p_RelementDistribution(1)%celement
!
!      ! Get the number of local DOF's for trial and test functions
!      indof = elem_igetNDofLoc(celement)
!      
!      ! Number of vertices on the element
!      nve = elem_igetNVE(celement)
!      
!      ! Type of transformation from/to the reference element
!      ctrafoType = elem_igetTrafoType(celement)
!      
!      ! Element nonparametric?
!      bnonpar = elem_isNonparametric(celement)
!      
!      NULLIFY(p_IelementDistr)
!    ELSE
!      CALL storage_getbase_int (rvectorScalar%p_rspatialDiscr%h_IelementDistr,&
!          p_IelementDistr)
!    END IF
!    
!    ! Get the data vector
!    SELECT CASE (rvectorScalar%cdataType)
!    CASE (ST_DOUBLE) 
!      CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
!    CASE (ST_SINGLE)
!      CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
!    CASE DEFAULT
!      CALL output_line ('Unsupported vector precision!',&
!          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
!      CALL sys_halt()
!    END SELECT
!    
!    ! What to evaluate?
!    Bder = .FALSE.
!    Bder(iderType) = .TRUE.
!    
!    ! Get the type of the element ielement
!    IF (ASSOCIATED(p_IelementDistr)) THEN
!      ! As all elements have the same type, we get the element
!      ! characteristics by checking the first element.
!      celement = p_RelementDistribution(p_IelementDistr(Ielements(1)))%celement
!
!      ! Get the number of local DOF's for trial and test functions
!      indof = elem_igetNDofLoc(celement)
!      
!      ! Number of vertices on the element
!      nve = elem_igetNVE(celement)
!      
!      ! Type of transformation from/to the reference element
!      ctrafoType = elem_igetTrafoType(celement)
!      
!      ! Element nonparametric?
!      bnonpar = elem_isNonparametric(celement)
!      
!    END IF
!      
!    ! Calculate the global DOF's on that element into IdofsTest.
!    ALLOCATE(Idofs(indof,SIZE(Ielements)))
!    CALL dof_locGlobMapping_mult(rvectorScalar%p_rspatialDiscr, &
!        Ielements, .FALSE.,Idofs)
!        
!    ! Get the coordinates forming the elements
!    ALLOCATE(Dcoord(UBOUND(p_DvertexCoords,1),NVE,SIZE(Ielements)))
!    DO iel = 1,SIZE(Ielements)
!      Dcoord(:,1:nve,iel) = &
!        p_DvertexCoords(:,p_IverticesAtElement(:,Ielements(iel)))
!    END DO
!    
!    ! Get the coordinates of all points on the reference element
!    IF (PRESENT(DpointsRef)) THEN
!      p_DpointsRef => DpointsRef
!    ELSE
!      ALLOCATE(p_DpointsRef(UBOUND(Dpoints,1),UBOUND(Dpoints,2),&
!                            UBOUND(Dpoints,3)))
!      ! Calculate the transformation of the point to the reference element
!      DO iel = 1,SIZE(Ielements)
!        DO ipoint = 1,UBOUND(Dpoints,2)
!          CALL trafo_calcRefCoords (ctrafoType,Dcoord(:,:,iel),&
!              Dpoints(:,ipoint,iel),p_DpointsRef(:,ipoint,iel))
!        END DO
!      END DO
!    END IF
!    
!    ! Calculate the transformation of all the points from the reference
!    ! to the real element(s).
!    ndim = UBOUND(Dcoord,1)
!    ALLOCATE(Djac(ndim*ndim,UBOUND(Dpoints,2),UBOUND(Dpoints,3)))
!    ALLOCATE(Ddetj(UBOUND(Dpoints,2),UBOUND(Dpoints,3)))
!    CALL trafo_calctrafo_sim (elem_igetTrafoType(celement),SIZE(Ielements),&
!        UBOUND(Dpoints,2),Dcoord,&
!        p_DpointsRef,Djac,Ddetj)    
!  
!    ! Does the element need twist indices?
!    ntwistsize = elem_getTwistIndexSize(celement)
!
!    ! If necessary, calculate the twist index array. The element may need it.
!    ! We always allocate so that we have something we can pass to the element...
!    ALLOCATE(ItwistIndex(MAX(1,ntwistSize),SIZE(Ielements)))
!    IF (ntwistSize .NE. 0) THEN
!      CALL trafo_calcTwistIndices(&
!        rvectorScalar%p_rspatialDiscr%p_rtriangulation,Ielements,ItwistIndex)
!    END IF
!
!    ! Calculate the values of the basis functions in the given points.
!    ALLOCATE(Dbas(indof,&
!             elem_getMaxDerivative(celement),&
!             UBOUND(Dpoints,2), UBOUND(Dpoints,3)))
!    IF (bnonpar) THEN
!      CALL elem_generic_sim (celement, Dcoord, Djac, Ddetj, &
!                           Bder, Dbas, UBOUND(Dpoints,2), UBOUND(Dpoints,3), &
!                           Dpoints,ItwistIndex)    
!    ELSE
!      CALL elem_generic_sim (celement, Dcoord, Djac, Ddetj, &
!                           Bder, Dbas, UBOUND(Dpoints,2), UBOUND(Dpoints,3), &
!                           p_DpointsRef,ItwistIndex)    
!    END IF
!  
!    ! Calculate the desired values. We loop over all points and all elements
!    IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
!      DO iel = 1, UBOUND(Dpoints,3)
!        DO ipoint = 1,UBOUND(Dpoints,2)
!      
!          dval = 0.0_DP
!        
!          ! Now that we have the basis functions, we want to have the function values.
!          ! We get them by multiplying the FE-coefficients with the values of the
!          ! basis functions and summing up.
!          !          
!          ! Calculate the value in the point
!          DO ibas = 1,indof
!            dval = dval + p_Ddata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
!          END DO
!        
!          ! Save the value in the point
!          Dvalues(ipoint,iel) = dval
!        
!        END DO ! ipoint
!      END DO ! iel
!
!    ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
!    
!      DO iel = 1, UBOUND(Dpoints,3)
!        DO ipoint = 1,UBOUND(Dpoints,2)
!          
!          ! Now that we have the basis functions, we want to have the function values.
!          ! We get them by multiplying the FE-coefficients with the values of the
!          ! basis functions and summing up.
!          !
!          ! Calculate the value in the point
!          DO ibas = 1,indof
!            dval = dval + p_Fdata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
!          END DO
!        
!          ! Save the value in the point
!          Dvalues(ipoint,iel) = dval
!        
!        END DO ! ipoint
!      END DO ! iel
!      
!    END IF
!    
!    ! Release allocated memory
!    DEALLOCATE(ItwistIndex)
!    DEALLOCATE(Dbas)
!    DEALLOCATE(Ddetj)
!    DEALLOCATE(Djac)
!    IF (.NOT. PRESENT(DpointsRef)) THEN
!      DEALLOCATE(p_DpointsRef)
!    END IF
!    DEALLOCATE(Dcoord)
!    DEALLOCATE(Idofs)
!
!  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_sim1 (iderType, Dvalues, rvectorScalar, Dpoints, &
      Ielements, DpointsRef)
                                      
!<description>
  ! This is a rather general finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints on a set of elements Ielements. 
  ! The values of the FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)              :: rvectorScalar
  
  ! A list of points where to evaluate. All points must be inside
  ! of element ielement.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints
  
  ! A list of elements containing the points in Dpoints.
  ! All elements in this list must be of the same type!!!
  integer, dimension(:), intent(IN) :: Ielements
  
  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  real(DP), dimension(:,:,:), intent(IN), target, optional :: DpointsRef
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(1..npoints,1..nelements)
  real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    logical :: bnonpar
    integer :: ipoint,indof,nve,ibas,iel
    integer(I32) :: celement
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval
    
    real(DP), dimension(:,:,:), pointer :: p_DpointsRef
    
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! Transformation
    integer(I32) :: ctrafoType
    
    ! Values of basis functions and DOF's
    real(DP), dimension(:,:,:,:), allocatable :: Dbas
    integer, dimension(:,:), allocatable :: Idofs

    ! Element evaluation set that collects element specific information
    ! during the evaluation
    type(t_evalElementSet)  :: revalElementSet
    integer(I32) :: cevaluationTag
    
    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! Ok, slow but general.
    
    ! Get triangulation information
    p_rtriangulation => rvectorScalar%p_rspatialDiscr%p_rtriangulation
    
    p_RelementDistribution => rvectorScalar%p_rspatialDiscr%RelementDistr

    ! For uniform discretisations, we get the element type in advance...
    if (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      
      ! Element type
      celement = p_RelementDistribution(1)%celement

      ! Get the number of local DOF's for trial and test functions
      indof = elem_igetNDofLoc(celement)
      
      ! Number of vertices on the element
      nve = elem_igetNVE(celement)
      
      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)
      
      ! Element nonparametric?
      bnonpar = elem_isNonparametric(celement)
      
      nullify(p_IelementDistr)
    else
      call storage_getbase_int (rvectorScalar%p_rspatialDiscr%h_IelementDistr,&
          p_IelementDistr)
    end if
    
    ! Get the data vector
    select case (rvectorScalar%cdataType)
    case (ST_DOUBLE) 
      call lsyssc_getbase_double(rvectorScalar,p_Ddata)
    case (ST_SINGLE)
      call lsyssc_getbase_single(rvectorScalar,p_Fdata)
    case DEFAULT
      call output_line ('Unsupported vector precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
      call sys_halt()
    end select
    
    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.
    
    ! Get the type of the element ielement
    if (associated(p_IelementDistr)) then
      ! As all elements have the same type, we get the element
      ! characteristics by checking the first element.
      celement = p_RelementDistribution(p_IelementDistr(Ielements(1)))%celement

      ! Get the number of local DOF's for trial and test functions
      indof = elem_igetNDofLoc(celement)
      
      ! Number of vertices on the element
      nve = elem_igetNVE(celement)
      
      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)
      
      ! Element nonparametric?
      bnonpar = elem_isNonparametric(celement)
      
    end if
      
    ! Calculate the global DOF's on that element into IdofsTest.
    allocate(Idofs(indof,size(Ielements)))
    call dof_locGlobMapping_mult(rvectorScalar%p_rspatialDiscr, &
        Ielements, Idofs)
    
    ! Initialisation of the element set.    
    call elprep_init(revalElementSet)   
        
    ! Get the coordinates of the corners of the elements
    call elprep_prepareSetForEvaluation (revalElementSet,&
        EL_EVLTAG_COORDS, p_rtriangulation, Ielements, ctrafoType,&
        DpointsRef=DpointsRef,DpointsReal=Dpoints)

    ! Get the coordinates of all points on the reference element
    if (present(DpointsRef)) then
      p_DpointsRef => DpointsRef
    else
      allocate(p_DpointsRef(ubound(Dpoints,1),ubound(Dpoints,2),&
               ubound(Dpoints,3)))
      ! Calculate the transformation of the point to the reference element
      do iel = 1,size(Ielements)
        do ipoint = 1,ubound(Dpoints,2)
          call trafo_calcRefCoords (ctrafoType,revalElementSet%p_Dcoords(:,:,iel),&
              Dpoints(:,ipoint,iel),p_DpointsRef(:,ipoint,iel))
        end do
      end do
    end if

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(celement)
    
    ! Don't create coordinates on the reference/real element; we do this manually!
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REALPOINTS))

    ! Don't calculate element shape information, we have that already.
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_COORDS))
                    
    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.

    ! Prepare the element set for the evaluation
    call elprep_prepareSetForEvaluation (revalElementSet,&
        cevaluationTag, p_rtriangulation, Ielements, ctrafoType,&
        DpointsRef=p_DpointsRef,DpointsReal=Dpoints)
    
    ! Calculate the values of the basis functions in the given points.
    allocate(Dbas(indof,&
             elem_getMaxDerivative(celement),&
             ubound(Dpoints,2), size(Ielements)))
    call elem_generic_sim2 (celement, revalElementSet, Bder, Dbas)
             
    ! Calculate the desired values. We loop over all points and all elements
    if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
      do iel = 1, size(Ielements)
        do ipoint = 1,ubound(Dpoints,2)
      
          dval = 0.0_DP
        
          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.
          !          
          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
          end do
        
          ! Save the value in the point
          Dvalues(ipoint,iel) = dval
        
        end do ! ipoint
      end do ! iel

    else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
    
      do iel = 1, size(Ielements)
        do ipoint = 1,ubound(Dpoints,2)
          
          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.
          !
          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
          end do
        
          ! Save the value in the point
          Dvalues(ipoint,iel) = dval
        
        end do ! ipoint
      end do ! iel
      
    end if
    
    ! Release allocated memory
    ! Remove the reference to DpointsRef again
    deallocate(Dbas)

    call elprep_releaseElementSet(revalElementSet)
    deallocate(Idofs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_sim2 (rvectorScalar, Dcoords, Djac, Ddetj, &
                  celement, IdofsTrial, npoints,  nelements, Dpoints, iderType,&
                  Dvalues,ItwistIndexEdges)
                                      
!<description>
  ! DEPRECATED!
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on multiple elements in a 
  ! discretisation.
  ! The caller must provide all necessary information for the evaluation in the
  ! parameters.
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating 'simultaneously' on some $Q_1$ and some $P_1$ elements
  ! is not allowed e.g.. 
!</description>

!<input>
  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)              :: rvectorScalar
  
  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those 
  ! elements that we are concerning here.
  integer(I32), intent(IN)                      :: celement

  ! A list of the corner vertices of all elements in progress.
  ! array [1..NDIM2D,1..TRIA_MAXNVE2D,1..Number of elements] of double
  real(DP), dimension(:,:,:), intent(IN)        :: Dcoords
  
  ! The Jacobian matrix of the mapping between the reference and each
  ! real element, for all points on all elements in progress.
  ! array [1..TRAFO_NJACENTRIES,1..npointsPerElement,1..Number of elements]
  real(DP), dimension(:,:,:),intent(IN)         :: Djac
  
  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to each real element in progress.
  ! array [1..npointsPerElement,1..Number of elements]
  real(DP), dimension(:,:), intent(IN)          :: Ddetj
  
  ! An array accepting the DOF's on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF's in trial space,nelements)
  integer, dimension(:,:), intent(IN)           :: IdofsTrial
  
  ! Number of points on every element where to evalate the function
  integer, intent(IN) :: npoints
  
  ! Number of elements, the function is evaluated at
  integer, intent(IN)  :: nelements
  
  ! Array with coordinates of the points where to evaluate.
  ! DIMENSION(NDIM2D,npoints,nelements).
  ! The coordinates are expected 
  ! - on the reference element, if celement identifies a parametric element
  ! - on the real element, if celement identifies a nonparametric element
  ! It's assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType

  ! OPTIONAL: List of twist indices. Defines for every element the
  ! orgientation of the edges.
  ! Can be omitted if the element does not need it.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  integer(I32), dimension(:), intent(IN), optional :: ItwistIndexEdges
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial
  real(DP) :: dval
  integer :: iel,ipoint,ibas
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  
  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))
  
  ! Evaluate the basis functions
  call elem_generic_sim (celement, Dcoords, Djac, Ddetj, &
                         Bder, DbasTrial, npoints, nelements, Dpoints,ItwistIndexEdges)  
  
  if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
  
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
                 p_Ddata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(ipoint,iel) = dval
      end do
    end do
    
  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
  
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
                 p_Fdata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(ipoint,iel) = dval
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

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  celement, IdofsTrial, iderType, Dvalues)
                                      
!<description>
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on multiple elements in a 
  ! discretisation.
  ! revalElementScalar must specify all information about where and how
  ! to evaluate; e.g. the coordinates of the evaluation points are
  ! to be found here. The routine will then evaluate the element celement
  ! in these points.
  !
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating 'simultaneously' on some $Q_1$ and some $P_1$ elements
  ! is not allowed e.g.. 
!</description>

!<input>
  ! Element evaluation set that contains all information necessary
  ! for the evaluation (coordinates of the points, transformation,...)
  type(t_evalElementSet), intent(IN) :: revalElementSet

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(IN)              :: rvectorScalar
  
  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those 
  ! elements that we are concerning here.
  integer(I32), intent(IN)                      :: celement

  ! An array accepting the DOF's on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF's in trial space,nelements)
  integer, dimension(:,:), intent(IN) :: IdofsTrial
  
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  real(DP) :: dval
  integer :: iel,ipoint,ibas
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  
  npoints = ubound(Dvalues,1)
  nelements = ubound(Dvalues,2)
  
  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))
  
  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, revalElementSet, Bder, DbasTrial)
  
  if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
  
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
                 p_Ddata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(ipoint,iel) = dval
      end do
    end do
    
  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
  
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
                 p_Fdata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(ipoint,iel) = dval
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

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_sim4 (rvectorScalar, &
                                 rdomainIntSubset, iderType, Dvalues, iterm)
                                      
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
  
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(IN)                            :: iderType
  
  ! Number of the subarray in Dvalues where the result is written to.
  ! The routine writes the resulting values to Dvalues(iterm,:,:).
  integer, intent(in)                            :: iterm
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(#possible terms,npoints,nelements).
  real(DP), dimension(:,:,:), intent(OUT) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  real(DP) :: dval
  integer :: iel,ipoint,ibas
  integer(I32) :: celement
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  integer, dimension(:,:), pointer :: p_IdofsTrial
  
  npoints = ubound(Dvalues,2)
  nelements = ubound(Dvalues,3)
  
  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.
  
  ! Get the currently active element
  celement = rvectorScalar%p_rspatialDiscr%RelementDistr( &
      rdomainIntSubset%ielementDistribution)%celement
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))
  
  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)
  
  ! Get the pointer to the trail DOF's.
  ! If the IdofsTrial in the domain subset fits to our current element,
  ! take that. Otherwise, we have to compute the actual DOF's.
  if (rdomainIntSubset%celement .eq. celement) then
    p_IdofsTrial => rdomainIntSubset%p_IdofsTrial
  else
    allocate (p_IdofsTrial(indofTrial,nelements))
    call dof_locGlobMapping_mult(rvectorScalar%p_rspatialDiscr, &
        rdomainIntSubset%p_Ielements, p_IdofsTrial)
  end if
  
  if (rvectorScalar%cdataType .eq. ST_DOUBLE) then
  
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
                 p_Ddata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(iterm,ipoint,iel) = dval
      end do
    end do
    
  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then
  
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
                 p_Fdata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        end do
        ! Save the value in the point
        Dvalues(iterm,ipoint,iel) = dval
      end do
    end do
    
  else
    call output_line('Unsupported vector precision!',&
      OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate_sim')
    call sys_halt()
  end if
  
  ! Release memory, finish
  deallocate(DbasTrial)
  
  if (rdomainIntSubset%celement .ne. celement) then
    deallocate (p_IdofsTrial)
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_getVectorMagnitude (rvector,dumax)
                                      
!<description>
  ! Routine to calculate the vector magnitude of a FE vector field.
!</description>

!<input>
  ! A given finite element vector field. May be e.g. a velocity field.
  type(t_vectorBlock), intent(IN)            :: rvector
!</input>

!<output>
  ! The maximum vector field length. E.g. if rvector is a velocity field,
  ! this returns the maximum velocity.
  real(DP) :: dumax
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dvalues,p_Dvalues2
    integer :: i,j,celement
    logical :: bevaluate

    ! There are some special cases that we can handle directly.
    ! If all FE spaces are primal spaces (P1,Q1,...), we can just
    ! calculate dumax using lsysbl_getVectorMagnitude.
    ! Otherwise, we project the vector into a Q1 space
    ! and get the maximum value from the values in the vertices.
    bevaluate = .true.
    
    do i=2,rvector%nblocks
      ! All vectors must have the same shape
      if (rvector%RvectorBlock(i)%p_rspatialDiscr%inumFESpaces .ne. &
          rvector%RvectorBlock(1)%p_rspatialDiscr%inumFESpaces) then
         bevaluate = .false.
         exit
      end if
      
      do j=1,rvector%RvectorBlock(i)%p_rspatialDiscr%inumFESpaces
        celement = rvector%RvectorBlock(i)%p_rspatialDiscr%RelementDistr(j)%celement
        if (celement .ne. &
            rvector%RvectorBlock(1)%p_rspatialDiscr%RelementDistr(j)%celement) then
           bevaluate = .false.
           exit
        end if
        
        select case (elem_getPrimaryElement(celement))
        case (EL_P0_1D,EL_P1_1D,EL_P2_1D,&
              EL_P0_2D,EL_P1_2D,EL_P2_2D,EL_P3_2D,EL_P1T_2D,&
              EL_Q0_2D,EL_Q1_2D,EL_Q2_2D,EL_Q3_2D,EL_Q1T_2D,&
              EL_P0_3D,EL_P1_3D,EL_P2_3D,EL_Q0_3D,EL_Q1_3D,EL_Q2_3D,EL_Q1T_3D)
        case default
           bevaluate = .false.
           exit
        end select
        
      end do
    end do
    
    if (.not. bevaluate) then
    
      ! Evaluate directly
      call lsysbl_getVectorMagnitude (rvector,dumax=dumax)
    
    else
    
      ! Calculate the VecMag value by projection into the vertices.
      ! We project each component separately and sum up in p_Dvalues2.
      call spdp_projectToVertices (rvector%RvectorBlock(1), p_Dvalues, DER_FUNC)
      
      ! If there is more than one block...
      if (rvector%nblocks .gt. 1) then
       
        ! Calculate val=sqrt(val^2 + val2^2 + val3^2 + ...)
       
        do i=1,size(p_Dvalues)
          p_Dvalues(i) = p_Dvalues(i)**2
        end do
      
        do i=2,rvector%nblocks
          call spdp_projectToVertices (rvector%RvectorBlock(i), p_Dvalues2, DER_FUNC)

          do j=1,size(p_Dvalues)
            p_Dvalues(j) = p_Dvalues(j) + p_Dvalues2(j)**2
          end do
        end do

        do i=1,size(p_Dvalues)
          p_Dvalues(i) = sqrt(p_Dvalues(i))
        end do

      end if
      
      ! Get dumax
      dumax = lalg_norm(p_Dvalues,LINALG_NORMMAX,n=size(p_Dvalues))
      
      deallocate(p_Dvalues,p_Dvalues2)
      
    end if

  end subroutine

end module

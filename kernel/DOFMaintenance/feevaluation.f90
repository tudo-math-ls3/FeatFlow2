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
!#
!# 5.) fevl_evaluateBdr1d
!#     -> Evaluate a FE funtion in the unique boundary point in 1D.
!#
!# 6.) fevl_evaluateBdr2d
!#     -> Evaluate a FE funtion in an arbitrary set of points on the
!#        boundary in 2D. Determines automatically the elements that
!#        contain the points.
!#
!# </purpose>
!##############################################################################

module feevaluation

!$use omp_lib
  use basicgeometry
  use boundary
  use boundaryaux
  use collection
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use element
  use elementbase
  use elementpreprocessing
  use elementpreprocessing
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use perfconfig
  use spatialdiscretisation
  use spdiscprojection
  use storage
  use transformation
  use triangulation
  use triasearch
  use extstdassemblyinfo

  implicit none

  private

!<constants>

!<constantblock description="Constants for the cnonmeshPoints parameter.">

  ! Points outside of the domain not allowed.
  integer, parameter, public :: FEVL_NONMESHPTS_NONE   = 0

  ! Points outside of the domain allowed, try to find an element nearby to evaluate.
  integer, parameter, public :: FEVL_NONMESHPTS_NEARBY = 1

  ! Points outside of the domain allowed, assume 0 as value there.
  integer, parameter, public :: FEVL_NONMESHPTS_ZERO   = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! An evaluation assembly structure that saves crucial data during the assembly.
  type t_fevlAssembly

    ! Current cubature settings
    type(t_stdCubatureData) :: rcubatureInfo

    ! Current FE evaluation settings
    type(t_stdFEBasisEvalData) :: rfeBasisEvalData

  end type

!</typeblock>

  public :: t_fevlAssembly

!</types>


  ! There are two functions fevl_evaluate, one for scalar functions, one for
  ! multivariate ones.
  interface fevl_evaluate
    module procedure fevl_evaluate1
    module procedure fevl_evaluate2
    module procedure fevl_evaluate3
  end interface

  ! There are multiple functions fevl_evaluate_mult which do the same --
  ! with different calling conventions and different complexities.
  interface fevl_evaluate_mult
    module procedure fevl_evaluate_mult1
    module procedure fevl_evaluate_mult2
  end interface

  public :: fevl_evaluate_mult
  public :: fevl_evaluate_mult1
  public :: fevl_evaluate_mult2

  ! There are multiple functions fevl_evaluate_sim which do the same --
  ! with different calling conventions and different complexities.
  interface fevl_evaluate_sim
    module procedure fevl_evaluate_sim1
    module procedure fevl_evaluate_sim2 ! DEPRECATED
    module procedure fevl_evaluate_sim3
    module procedure fevl_evaluate_sim4
    module procedure fevl_evaluate_sim5
    module procedure fevl_evaluate_generic_sim1
    module procedure fevl_evaluate_generic_sim3
    module procedure fevl_evaluate_generic_sim4
    module procedure fevl_evaluate_generic_sim5
  end interface

  public :: fevl_evaluate
  public :: fevl_evaluate_sim
  public :: fevl_evaluate_sim1
  public :: fevl_evaluate_sim2 ! DEPRECATED
  public :: fevl_evaluate_sim3
  public :: fevl_evaluate_sim4
  public :: fevl_evaluate_sim5
  public :: fevl_evaluate_generic_sim1
  public :: fevl_evaluate_generic_sim3
  public :: fevl_evaluate_generic_sim4
  public :: fevl_evaluate_generic_sim5

  ! There are multiple functions fevl_evaluateBdr1d which do the same --
  ! with different calling conventions and different complexities.
  interface fevl_evaluateBdr1d
    module procedure fevl_evaluateBdr1d1
    module procedure fevl_evaluateBdr1d2
    module procedure fevl_evaluateBdr1d3
    module procedure fevl_evaluateBdr1d4
  end interface

  public :: fevl_evaluateBdr1d
  public :: fevl_evaluateBdr1d1
  public :: fevl_evaluateBdr1d2
  public :: fevl_evaluateBdr1d3
  public :: fevl_evaluateBdr1d4

  ! There are multiple functions fevl_evaluateBdr2d which do the same --
  ! with different calling conventions and different complexities.
  interface fevl_evaluateBdr2d
    module procedure fevl_evaluateBdr2d1
    module procedure fevl_evaluateBdr2d2
    module procedure fevl_evaluateBdr2d3
    module procedure fevl_evaluateBdr2d4
  end interface

  public :: fevl_evaluateBdr2d
  public :: fevl_evaluateBdr2d1
  public :: fevl_evaluateBdr2d2
  public :: fevl_evaluateBdr2d3
  public :: fevl_evaluateBdr2d4

  public :: fevl_getVectorMagnitude

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate1 (iderType, Dvalues, rvectorScalar, Dpoints, &
      Ielements, IelementsHint, cnonmeshPoints)

!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues. The routine is called via
  ! the interface fevl_evaluate(...).
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A list of points where to evaluate.
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(:,:), intent(in) :: Dpoints

  ! OPTIONAL: A list of elements containing the points Dpoints.
  ! If this is not specified, the element numbers containing the points
  ! are determined automatically.
  integer, dimension(:), intent(in), optional :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements
  ! containing the points. This is ignored if Ielements is specified!
  integer, dimension(:), intent(in), optional :: IelementsHint

  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
  ! if a point is located outside of the domain. May happen e.g. in
  ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default
  ! parameter if cnonmeshPoints is not specified.
  integer, intent(in), optional :: cnonmeshPoints

!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: cnonmesh
    integer :: ipoint,indof,nve,ibas
    integer(I32) :: celement
    integer :: iel,iresult,iellast,ivar
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! Transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint

    ! Values of basis functions and DOF`s
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
      celement = rvectorScalar%p_rspatialDiscr%RelementDistr(1)%celement

      ! Get the number of local DOF`s for trial and test functions
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
      call output_line ("Unsupported vector precision!", OU_CLASS_ERROR, OU_MODE_STD, &
                        "fevl_evaluate1")
      call sys_halt()
    end select

    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.

    cnonmesh = FEVL_NONMESHPTS_NONE
    if (present(cnonmeshPoints)) cnonmesh = cnonmeshPoints

    ! We loop over all points.

    iel = 1

    !$omp parallel do default(shared) firstprivate(iel)&
    !$omp private(Dbas,DparPoint,Idofs,celement,cevaluationTag,ctrafoType,&
    !$omp         dval,ibas,iellast,indof,iresult,nve,revalElement)
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
          Dpoints(:,ipoint),rvectorScalar%p_rspatialDiscr%p_rtriangulation,iel,&
          iresult=iresult,ilastElement=iellast)

        if (iresult .eq. -2) then
          ! Not found, too many iterations. Probably the domain is too big.
          ! Search again starting from the last found element, specify
          ! 2*sqrt(NEL) as maximum number of steps (should work for rather
          ! isotropic meshes) and try again.
          iel = iellast
          call tsrch_getElem_raytrace2D (&
            Dpoints(:,ipoint),rvectorScalar%p_rspatialDiscr%p_rtriangulation,iel,&
            iresult=iresult,ilastElement=iellast,&
            imaxIterations=max(100,2*int(sqrt(real(rvectorScalar%p_rspatialDiscr%p_rtriangulation%NEL,DP)))))
        end if

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
        call output_line ("Point "//trim(sys_siL(ipoint,10))//" not found!", &
                          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate1")
        cycle
      end if

      ! Get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! Get the number of local DOF`s for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! Number of vertices on the element
        nve = elem_igetNVE(celement)

        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)
      end if

      ! Calculate the global DOF`s on that element into IdofsTest.
      call dof_locGlobMapping (rvectorScalar%p_rspatialDiscr, iel, Idofs)

      ! Get the element shape information
      call elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
          rvectorScalar%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)

      ! Calculate the transformation of the point to the reference element
      call trafo_calcRefCoords (ctrafoType, revalElement%Dcoords,&
          Dpoints(:,ipoint), DparPoint)

      ! Now calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
          iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
          rvectorScalar%p_rspatialDiscr%p_rtriangulation, iel, &
          ctrafoType, DparPoint, Dpoints(:,ipoint))

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! Combine the basis functions to get the function value.
      if (rvectorScalar%cdataType .eq. ST_DOUBLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the
          ! function values.  We get them by multiplying the
          ! FE-coefficients with the values of the basis functions and
          ! summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(ipoint) = dval

        else

          do ivar = 1, rvectorScalar%NVAR
            ! Now that we have the basis functions, we want to have the
            ! function values.  We get them by multiplying the
            ! FE-coefficients with the values of the basis functions and
            ! summing up.

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Ddata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval
          end do

        end if

      else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the
          ! function values.  We get them by multiplying the
          ! FE-coefficients with the values of the basis functions and
          ! summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

        else

          do ivar = 1, rvectorScalar%NVAR
            ! Now that we have the basis functions, we want to have the
            ! function values.  We get them by multiplying the
            ! FE-coefficients with the values of the basis functions and
            ! summing up.

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Fdata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval
          end do

        end if

      end if

    end do ! ipoint

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine fevl_evaluate2 (CderType, Dvalues, rvectorBlock, Dpoints, Ielements, &
                             IelementsHint, cnonmeshPoints, iblockMin, iblockMax)

!<description>
  ! This subroutine is an extension of fevl_evaluate1(). It simultaneously evaluates
  ! several types of function values (given by CderType) for all components of a
  ! multivariate FE function (given by rvectorBlock). Thus, the evaluation of all the
  ! element basis functions etc. has to be performed only once. The routine is called via
  ! the interface fevl_evaluate(...).
  ! E.g., instead of,
  !
  !   real(DP), dimension(nblocks,npoints) :: Dvalues
  !   call fevl_evaluate(DER_FUNC, Dvalues(1,:), rsol%RvectorBlock(1), Dpoints)
  !   call fevl_evaluate(DER_FUNC, Dvalues(2,:), rsol%RvectorBlock(2), Dpoints)
  !   call fevl_evaluate(DER_DERIV_X, DderivX(1,:), rsol%RvectorBlock(1), Dpoints)
  !   call fevl_evaluate(DER_DERIV_X, DderivX(2,:), rsol%RvectorBlock(2), Dpoints)
  !   call fevl_evaluate(DER_DERIV_Y, DderivY(1,:), rsol%RvectorBlock(1), Dpoints)
  !   call fevl_evaluate(DER_DERIV_Y, DderivY(2,:), rsol%RvectorBlock(2), Dpoints)
  !
  ! one can use
  !
  !   integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)
  !   real(DP), dimension(nblocks,3,npoints) :: Dvalues
  !   call fevl_evaluate(CderType, Dvalues, rsol, Dpoints)
  !
!</description>

!<input>
  ! array of type of function values to evaluate (DER_FUNC, DER_DERIV_X etc.)
  integer, dimension(:), intent(in) :: CderType

  ! block solution vector representing the FE function that is to be evaluated
  ! It is assumed that all components correspond to the same spatial discretisation and
  ! use the same element type (i.e., only that of the first component is inquired)
  type(t_vectorBlock), intent(in) :: rvectorBlock

  ! list of points where to evaluate. DIMENSION(ndim,npoints)
  real(DP), dimension(:,:), intent(in) :: Dpoints

  ! OPTIONAL: A list of elements containing the points Dpoints.
  ! If this is not specified, the element numbers containing the points are determined
  ! automatically.
  integer, dimension(:), intent(in), optional :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements containing
  ! the points. This is ignored if Ielements is specified!
  integer, dimension(:), intent(in), optional :: IelementsHint

  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens if a point is
  ! located outside of the domain. May happen e.g. in nonconvex domains.
  ! FEVL_NONMESHPTS_NONE is the default parameter if cnonmeshPoints is not specified.
  integer, intent(in), optional :: cnonmeshPoints

  ! OPTIONAL: For the case that not all components of the vector are to be processed, the
  ! user can provide these two parameters so that only the components from iblockMin to
  ! iblockMax are processed. This can be necessary, e.g., in case of a saddle point
  ! problem where the last component of the solution vector (=pressure) is discretised
  ! in a different way than the first components (=velocity/displacement).
  ! If iblockMin (iblockMax) is not present, the minimal (maximal) block number is
  ! assumed to be 1 (rvectorBlock%nblocks).
  integer, intent(in), optional :: iblockMin, iblockMax
!</input>

!<output>
  ! values of the FE function at the points specified by Dpoints.
  ! DIMENSION(rvectorBlock%nblocks, size(CderType), npoints)
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: cnonmesh
    integer :: ipoint, indof, nve, ibas, ider, iblock, iblMin, iblMax
    integer(I32) :: celement
    integer :: iel,iresult,iellast
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval

    ! number of equations in one scalar component of the block vector
    integer :: neqsc

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint

    ! values of basis functions and DOF`s
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    integer, dimension(EL_MAXNBAS) :: Idofs

    ! pointer to the list of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! pointer to the spatial discretisation structure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    ! evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag

    if (present(iblockMin)) then
      iblMin = iblockMin
    else
      iblMin = 1
    endif
    if (present(iblockMax)) then
      iblMax = iblockMax
    else
      iblMax = rvectorBlock%nblocks
    endif

    ! set shortcuts
    p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(iblMin)
    p_RelementDistribution => p_rspatialDiscr%RelementDistr

    if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      ! for uniform discretisations, get the element type in advance

      ! element type
      celement = p_RelementDistribution(1)%celement

      ! number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! number of vertices on the element
      nve = elem_igetNVE(celement)

      ! type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)

      nullify(p_IelementDistr)
    else
      ! non-uniform discretisations
      call storage_getbase_int (p_rspatialDiscr%h_IelementDistr, p_IelementDistr)
    end if

    ! set pointer to the data vector
    select case (rvectorBlock%cdataType)
    case (ST_DOUBLE)
      call lsysbl_getbase_double(rvectorBlock, p_Ddata)
    case (ST_SINGLE)
      call lsysbl_getbase_single(rvectorBlock, p_Fdata)
    case DEFAULT
      call output_line ("Unsupported vector precision!", OU_CLASS_ERROR, OU_MODE_STD, &
                        "fevl_evaluate2")
      call sys_halt()
    end select

    ! inquire which types of function values are to be evaluated
    Bder = .false.
    do ider = 1, size(CderType)
      Bder(CderType(ider)) = .true.
    end do

    cnonmesh = FEVL_NONMESHPTS_NONE
    if (present(cnonmeshPoints)) then
      cnonmesh = cnonmeshPoints
    endif

    iel = 1

    ! loop over all points
    do ipoint = 1,ubound(Dpoints,2)

      ! get the element number that contains the point
      if (present(Ielements)) then
        ! either it is provided...
        iel = Ielements (ipoint)
      else
        ! ...or the element has to be searched

        if (present(IelementsHint)) then
          ! if available, use the hint provided by the user
          iel = IelementsHint (ipoint)
        end if
        ! otherwise, we use iel from the previous iteration as inital guess

        ! use raytracing search to find the element containing the point
        call tsrch_getElem_raytrace2D(Dpoints(:,ipoint), p_rspatialDiscr%p_rtriangulation, &
          iel, iresult=iresult,ilastElement=iellast)

        if (iresult .eq. -2) then
          ! Not found, too many iterations. Probably the domain is too big.
          ! Search again starting from the last found element, specify
          ! 2*sqrt(NEL) as maximum number of steps (should work for rather
          ! isotropic meshes) and try again.
          iel = iellast
          call tsrch_getElem_raytrace2D (&
            Dpoints(:,ipoint), p_rspatialDiscr%p_rtriangulation,iel,&
            iresult=iresult,ilastElement=iellast,&
            imaxIterations=max(100,2*int(sqrt(real(p_rspatialDiscr%p_rtriangulation%NEL,DP)))))
        end if

        if (iel .eq. 0) then
          ! if not found, use brute force search
          call tsrch_getElem_BruteForce(Dpoints(:,ipoint), &
                                        p_rspatialDiscr%p_rtriangulation, iel)
        end if

        if (iel .eq. 0) then
          ! if element is still not found, inquire if nonmesh-points are allowed
          if (cnonmesh .eq. FEVL_NONMESHPTS_NEARBY) then
            ! if yes then find the closest element
            call tsrch_getNearestElem_BruteForce(Dpoints(:,ipoint), &
                   p_rspatialDiscr%p_rtriangulation, iel)
            ! The evaluation routine then computes the FE function outside of the element.

          else if (cnonmesh .eq. FEVL_NONMESHPTS_ZERO) then
            ! if no then set the value to zero
            Dvalues(iblMin:iblMax,:,ipoint) = 0.0_DP
            ! go to next iteration
            cycle
          end if
        end if
      end if

      if (iel .eq. 0) then
        ! report when element has not been found and go to next iteration
        call output_line ("Point " // trim(sys_siL(ipoint,10)) // " not found!", &
                          OU_CLASS_ERROR, OU_MODE_STD, "fevl_evaluate2")
        Dvalues(iblMin:iblMax,:,ipoint) = 0.0_DP
        cycle
      end if

      ! get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! get the number of local DOFs for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! number of vertices on the element
        nve = elem_igetNVE(celement)

        ! type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)
      end if

      ! calculate the global DOFs on that element into IdofsTest.
      call dof_locGlobMapping (p_rspatialDiscr, iel, Idofs)

      ! get the element shape information
      call elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
             p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)

      ! calculate the transformation of the point to the reference element
      call trafo_calcRefCoords (ctrafoType, revalElement%Dcoords, &
                                Dpoints(:,ipoint), DparPoint)

      ! calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
             iand(cevaluationTag,not(EL_EVLTAG_COORDS)), p_rspatialDiscr%p_rtriangulation, &
             iel, ctrafoType, DparPoint, Dpoints(:,ipoint))

      ! call the element to calculate the values of the basis functions in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! number of equations in one scalar component
      neqsc = rvectorBlock%RvectorBlock(iblMin)%neq

      ! calculate function values by multiplying the FE-coefficients with the values of
      ! the basis functions and summing up
      Dvalues(iblMin:iblMax,:,ipoint) = 0.0_DP
      if (rvectorBlock%cdataType .eq. ST_DOUBLE) then
        do ider = 1,size(CderType)
          do iblock = iblMin,iblMax
            dval = 0.0_DP
            do ibas = 1,indof
              dval = dval +   p_Ddata((iblock-1)*neqsc + Idofs(ibas)) &
                            * Dbas(ibas,CderType(ider))
            end do
            Dvalues(iblock, ider, ipoint) = dval
          end do
        end do
      else if (rvectorBlock%cdataType .eq. ST_SINGLE) then
        do ider = 1,size(CderType)
          do iblock = iblMin,iblMax
            dval = 0.0_DP
            do ibas = 1,indof
              dval = dval +   p_Fdata((iblock-1)*neqsc + Idofs(ibas)) &
                            * Dbas(ibas,CderType(ider))
            end do
            Dvalues(iblock, ider, ipoint) = dval
          end do
        end do
      end if

    end do ! ipoint

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate3 (CderType, Dvalues, rvector, Dpoints, Ielements, &
      IelementsHint, cnonmeshPoints)

!<description>
  ! This subroutine is an extension of fevl_evaluate1(). It simultaneously evaluates
  ! several types of function values (given by CderType) for all components of a
  ! FE function (given by rvector). Thus, the evaluation of all the
  ! element basis functions etc. has to be performed only once. The routine is called via
  ! the interface fevl_evaluate(...).
  ! E.g., instead of,
  !
  !   real(DP), dimension(3,npoints) :: Dvalues
  !   call fevl_evaluate(DER_FUNC, Dvalues(1,:), rsol, Dpoints)
  !   call fevl_evaluate(DER_DERIV_X, DderivX(1,:), rsol, Dpoints)
  !   call fevl_evaluate(DER_DERIV_Y, DderivY(1,:), rsol, Dpoints)
  !
  ! one can use
  !
  !   integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)
  !   real(DP), dimension(3,npoints) :: Dvalues
  !   call fevl_evaluate(CderType, Dvalues, rsol, Dpoints)
  !
!</description>

!<input>
  ! array of type of function values to evaluate (DER_FUNC, DER_DERIV_X etc.)
  integer, dimension(:), intent(in) :: CderType

  ! Solution vector representing the FE function that is to be evaluated
  ! It is assumed that all components correspond to the same spatial discretisation and
  ! use the same element type (i.e., only that of the first component is inquired)
  type(t_vectorScalar), intent(in) :: rvector

  ! list of points where to evaluate. DIMENSION(ndim,npoints)
  real(DP), dimension(:,:), intent(in) :: Dpoints

  ! OPTIONAL: A list of elements containing the points Dpoints.
  ! If this is not specified, the element numbers containing the points are determined
  ! automatically.
  integer, dimension(:), intent(in), optional :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements containing
  ! the points. This is ignored if Ielements is specified!
  integer, dimension(:), intent(in), optional :: IelementsHint

  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens if a point is
  ! located outside of the domain. May happen e.g. in nonconvex domains.
  ! FEVL_NONMESHPTS_NONE is the default parameter if cnonmeshPoints is not specified.
  integer, intent(in), optional :: cnonmeshPoints
!</input>

!<output>
  ! values of the FE function at the points specified by Dpoints.
  ! DIMENSION(size(CderType), npoints)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: cnonmesh
    integer :: ipoint, indof, nve, ibas, ider
    integer(I32) :: celement
    integer :: iel,iresult,iellast
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint

    ! values of basis functions and DOF`s
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    integer, dimension(EL_MAXNBAS) :: Idofs

    ! pointer to the list of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! pointer to the spatial discretisation structure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    ! evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag

    ! set shortcuts
    p_rspatialDiscr => rvector%p_rspatialDiscr
    p_RelementDistribution => p_rspatialDiscr%RelementDistr

    if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      ! for uniform discretisations, get the element type in advance

      ! element type
      celement = p_RelementDistribution(1)%celement

      ! number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! number of vertices on the element
      nve = elem_igetNVE(celement)

      ! type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)

      nullify(p_IelementDistr)
    else
      ! non-uniform discretisations
      call storage_getbase_int (p_rspatialDiscr%h_IelementDistr, p_IelementDistr)
    end if

    ! set pointer to the data vector
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rvector, p_Ddata)
    case (ST_SINGLE)
      call lsyssc_getbase_single(rvector, p_Fdata)
    case default
      call output_line ("Unsupported vector precision!", OU_CLASS_ERROR, OU_MODE_STD, &
                        "fevl_evaluate3")
      call sys_halt()
    end select

    ! inquire which types of function values are to be evaluated
    Bder = .false.
    do ider = 1, size(CderType)
      Bder(CderType(ider)) = .true.
    end do

    cnonmesh = FEVL_NONMESHPTS_NONE
    if (present(cnonmeshPoints)) then
      cnonmesh = cnonmeshPoints
    endif

    iel = 1

    ! loop over all points
    do ipoint = 1,ubound(Dpoints,2)

      ! get the element number that contains the point
      if (present(Ielements)) then
        ! either it is provided...
        iel = Ielements (ipoint)
      else
        ! ...or the element has to be searched

        if (present(IelementsHint)) then
          ! if available, use the hint provided by the user
          iel = IelementsHint (ipoint)
        end if
        ! otherwise, we use iel from the previous iteration as inital guess

        ! use raytracing search to find the element containing the point
        call tsrch_getElem_raytrace2D(Dpoints(:,ipoint), p_rspatialDiscr%p_rtriangulation, &
          iel, iresult=iresult,ilastElement=iellast)

        if (iresult .eq. -2) then
          ! Not found, too many iterations. Probably the domain is too big.
          ! Search again starting from the last found element, specify
          ! 2*sqrt(NEL) as maximum number of steps (should work for rather
          ! isotropic meshes) and try again.
          iel = iellast
          call tsrch_getElem_raytrace2D (&
            Dpoints(:,ipoint), p_rspatialDiscr%p_rtriangulation,iel,&
            iresult=iresult,ilastElement=iellast,&
            imaxIterations=max(100,2*int(sqrt(real(p_rspatialDiscr%p_rtriangulation%NEL,DP)))))
        end if

        if (iel .eq. 0) then
          ! if not found, use brute force search
          call tsrch_getElem_BruteForce(Dpoints(:,ipoint), &
                                        p_rspatialDiscr%p_rtriangulation, iel)
        end if

        if (iel .eq. 0) then
          ! if element is still not found, inquire if nonmesh-points are allowed
          if (cnonmesh .eq. FEVL_NONMESHPTS_NEARBY) then
            ! if yes then find the closest element
            call tsrch_getNearestElem_BruteForce(Dpoints(:,ipoint), &
                   p_rspatialDiscr%p_rtriangulation, iel)
            ! The evaluation routine then computes the FE function outside of the element.

          else if (cnonmesh .eq. FEVL_NONMESHPTS_ZERO) then
            ! if no then set the value to zero
            Dvalues(:,ipoint) = 0.0_DP
            ! go to next iteration
            cycle
          end if
        end if
      end if

      if (iel .eq. 0) then
        ! report when element has not been found and go to next iteration
        call output_line ("Point " // trim(sys_siL(ipoint,10)) // " not found!", &
                          OU_CLASS_ERROR, OU_MODE_STD, "fevl_evaluate3")
        Dvalues(:,ipoint) = 0.0_DP
        cycle
      end if

      ! get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! get the number of local DOFs for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! number of vertices on the element
        nve = elem_igetNVE(celement)

        ! type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)
      end if

      ! calculate the global DOFs on that element into IdofsTest.
      call dof_locGlobMapping (p_rspatialDiscr, iel, Idofs)

      ! get the element shape information
      call elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
             p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)

      ! calculate the transformation of the point to the reference element
      call trafo_calcRefCoords (ctrafoType, revalElement%Dcoords, &
                                Dpoints(:,ipoint), DparPoint)

      ! calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
             iand(cevaluationTag,not(EL_EVLTAG_COORDS)), p_rspatialDiscr%p_rtriangulation, &
             iel, ctrafoType, DparPoint, Dpoints(:,ipoint))

      ! call the element to calculate the values of the basis functions in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! calculate function values by multiplying the FE-coefficients with the values of
      ! the basis functions and summing up
      Dvalues(:,ipoint) = 0.0_DP
      if (rvector%cdataType .eq. ST_DOUBLE) then
        do ider = 1,size(CderType)
          dval = 0.0_DP
          do ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,CderType(ider))
          end do
          Dvalues(ider, ipoint) = dval
        end do
      else if (rvector%cdataType .eq. ST_SINGLE) then
        do ider = 1,size(CderType)
          dval = 0.0_DP
          do ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,CderType(ider))
          end do
          Dvalues(ider, ipoint) = dval
        end do
      end if

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
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! The element number containing the points in Dpoints.
  integer, intent(in) :: ielement

  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! Either Dpoints or DpointsRef must be specified!
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(:,:), intent(in), optional :: DpointsRef

  ! OPTIONAL: A list of points where to evaluate. All points must be inside
  ! of element ielement.
  ! If not specified, the coordinates are automatically calculated.
  ! Either Dpoints or DpointsRef must be specified!
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(:,:), intent(in), optional :: Dpoints
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: ipoint,indof,nve,ibas,npoints,ivar
    integer(I32) :: celement
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP) :: dval

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! Transformation
    integer(I32) :: ctrafoType
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint

    ! Values of basis functions and DOF`s
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
      call output_line ("Evaluation points not specified!", &
                        OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_mult")
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

    ! Get the number of local DOF`s for trial and test functions
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
      call output_line ("Unsupported vector precision!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate")
      call sys_halt()
    end select

    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.

    ! Calculate the global DOF`s on that element into IdofsTest.
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
      ! Do not calculate the shape of the cell again since we did this in advance above.
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

      if (rvectorScalar%cdataType .eq. ST_DOUBLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(ipoint) = dval

        else

          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.

          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Ddata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval

          end do

        end if

      else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(ipoint) = dval

        else

          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.

          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Fdata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval

          end do

        end if

      end if

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
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those
  ! elements that we are concerning here.
  integer(I32), intent(in) :: celement

  ! A list of the corner vertices of the element.
  ! array [1..NDIM2D,1..TRIA_MAXNVE2D] of double
  real(DP), dimension(:,:), intent(in) :: Dcoords

  ! The Jacobian matrix of the mapping between the reference and the
  ! real element, for all points on the element.
  ! array [1..TRAFO_NJACENTRIES,1..npointsPerElement]
  real(DP), dimension(:,:),intent(in) :: Djac

  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to the real element.
  ! array [1..npointsPerElement]
  real(DP), dimension(:), intent(in) :: Ddetj

  ! An array accepting the DOF`s on the element in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF`s in trial space)
  integer, dimension(:), intent(in) :: IdofsTrial

  ! Number of points on the element where to evalate the function
  integer, intent(in) :: npoints

  ! Array with coordinates of the points where to evaluate.
  ! DIMENSION(NDIM2D,npoints).
  ! The coordinates are expected
  ! - on the reference element, if celement identifies a parametric element
  ! - on the real element, if celement identifies a nonparametric element
  ! It is assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:), intent(in) :: Dpoints

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! OPTIONAL: Twist index bitfield of the element. Defines for every edge its
  ! orientation.
  ! Can be omitted if the element does not need this information.
  integer(I32), intent(in), optional :: itwistIndex

!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints).
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:), allocatable :: DbasTrial
  integer :: indofTrial
  real(DP) :: dval
  integer :: ipoint,ibas,ivar
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

#if WARN_DEPREC
  call output_line ("Using deprecated feature. Please update your code.", &
      OU_CLASS_WARNING,OU_MODE_STD,"fevl_evaluate_mult2")
#endif

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

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

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

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

      do ipoint = 1,npoints
        ! Calculate the values in the point
        do ivar = 1,rvectorScalar%NVAR

          dval = 0.0_DP
          do ibas = 1,indofTrial
            dval = dval + p_Ddata((IdofsTrial(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                          DbasTrial(ibas,iderType,ipoint)
          end do
          ! Save the value in the point
          Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval
        end do
      end do

    end if

  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

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

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

      do ipoint = 1,npoints
        ! Calculate the values in the point
        do ivar = 1,rvectorScalar%NVAR

          dval = 0.0_DP
          do ibas = 1,indofTrial
            dval = dval + p_Fdata((IdofsTrial(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                          DbasTrial(ibas,iderType,ipoint)
          end do
          ! Save the value in the point
          Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval
        end do
      end do

    end if

  else
    call output_line("Unsupported vector precision!",&
      OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_mult2")
    call sys_halt()
  end if

  ! Release memory, finish
  deallocate(DbasTrial)

  end subroutine

!  ! ***************************************************************************
!<!-- // hide from automatic documentation parser
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
!  INTEGER, INTENT(in)                            :: iderType
!
!  ! The scalar solution vector that is to be evaluated.
!  TYPE(t_vectorScalar), INTENT(in)              :: rvectorScalar
!
!  ! A list of points where to evaluate. All points must be inside
!  ! of element ielement.
!  ! DIMENSION(1..ndim,1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:,:), INTENT(in) :: Dpoints
!
!  ! A list of elements containing the points in Dpoints.
!  ! All elements in this list must be of the same type!!!
!  INTEGER, DIMENSION(:), INTENT(in) :: Ielements
!
!  ! OPTIONAL: Coordinates of the points on the reference element.
!  ! If not specified, the coordinates are automatically calculated.
!  ! DIMENSION(1..ndim,1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:,:), INTENT(in), TARGET, OPTIONAL :: DpointsRef
!!</input>
!
!!<output>
!  ! Values of the FE function at the points specified by Dpoints.
!  ! DIMENSION(1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:), INTENT(out) :: Dvalues
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
!    ! Values of basis functions and DOF`s
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
!      ! Get the number of local DOF`s for trial and test functions
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
!      CALL output_line ("Unsupported vector precision!",&
!          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate")
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
!      ! Get the number of local DOF`s for trial and test functions
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
!    ! Calculate the global DOF`s on that element into IdofsTest.
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
! // unhide from automatic documentation parser -->

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_sim1 (iderType, Dvalues, rvectorScalar, Dpoints, &
      Ielements, DpointsRef, rperfconfig)

!<description>
  ! This is a rather general finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints on a set of elements Ielements.
  ! The values of the FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A list of points where to evaluate. All points must be inside
  ! of element ielement.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! A list of elements containing the points in Dpoints.
  ! All elements in this list must be of the same type!!!
  integer, dimension(:), intent(in) :: Ielements

  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  real(DP), dimension(:,:,:), intent(in), target, optional :: DpointsRef

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(1..npoints,1..nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    logical :: bnonpar
    integer :: ipoint,indof,nve,ibas,iel,ivar
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

    ! Values of basis functions and DOF`s
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

      ! Get the number of local DOF`s for trial and test functions
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
      call output_line ("Unsupported vector precision!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_sim1")
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

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Element nonparametric?
      bnonpar = elem_isNonparametric(celement)

    end if

    ! Calculate the global DOF`s on that element into IdofsTest.
    allocate(Idofs(indof,size(Ielements)))
    call dof_locGlobMapping_mult(rvectorScalar%p_rspatialDiscr, &
        Ielements, Idofs)

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Get the coordinates of the corners of the elements
    call elprep_prepareSetForEvaluation (revalElementSet,&
        EL_EVLTAG_COORDS, p_rtriangulation, Ielements, ctrafoType,&
        DpointsRef=DpointsRef,DpointsReal=Dpoints,rperfconfig=rperfconfig)

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

    ! Do not create coordinates on the reference/real element; we do this manually!
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REALPOINTS))

    ! Do not calculate element shape information, we have that already.
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_COORDS))

    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.

    ! Prepare the element set for the evaluation
    call elprep_prepareSetForEvaluation (revalElementSet,&
        cevaluationTag, p_rtriangulation, Ielements, ctrafoType,&
        DpointsRef=p_DpointsRef,DpointsReal=Dpoints,rperfconfig=rperfconfig)

    ! Calculate the values of the basis functions in the given points.
    allocate(Dbas(indof,&
             elem_getMaxDerivative(celement),&
             ubound(Dpoints,2), size(Ielements)))
    call elem_generic_sim2 (celement, revalElementSet, Bder, Dbas)

    ! Calculate the desired values. We loop over all points and all elements
    if (rvectorScalar%cdataType .eq. ST_DOUBLE) then

      ! Check if vector is stored in interleaved format
      if (rvectorScalar%NVAR .eq. 1) then

        ! Now that we have the basis functions, we want to have the
        ! function values.  We get them by multiplying the
        ! FE-coefficients with the values of the basis functions and
        ! summing up.

        do iel = 1, size(Ielements)
          do ipoint = 1,ubound(Dpoints,2)

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Ddata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues(ipoint,iel) = dval

          end do ! ipoint
        end do ! iel

      else

        ! Now that we have the basis functions, we want to have the
        ! function values.  We get them by multiplying the
        ! FE-coefficients with the values of the basis functions and
        ! summing up.

        do iel = 1,size(Ielements)
          do ipoint = 1,ubound(Dpoints,2)
            do ivar = 1,rvectorScalar%NVAR

              dval = 0.0_DP

              ! Calculate the value in the point
              do ibas = 1,indof
                dval = dval + p_Ddata((Idofs(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                              Dbas(ibas,iderType,ipoint,iel)
              end do

              ! Save the value in the point
              Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

            end do ! ivar
          end do ! ipoint
        end do ! iel
      end if

    else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

      ! Check if vector is stored in interleaved format
      if (rvectorScalar%NVAR .eq. 1) then

        ! Now that we have the basis functions, we want to have the
        ! function values.  We get them by multiplying the
        ! FE-coefficients with the values of the basis functions and
        ! summing up.

        do iel = 1, size(Ielements)
          do ipoint = 1,ubound(Dpoints,2)

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Fdata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues(ipoint,iel) = dval

          end do ! ipoint
        end do ! iel

      else

        ! Now that we have the basis functions, we want to have the
        ! function values.  We get them by multiplying the
        ! FE-coefficients with the values of the basis functions and
        ! summing up.

        do iel = 1,size(Ielements)
          do ipoint = 1,ubound(Dpoints,2)
            do ivar = 1,rvectorScalar%NVAR

              dval = 0.0_DP

              ! Calculate the value in the point
              do ibas = 1,indof
                dval = dval + p_Fdata((Idofs(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                              Dbas(ibas,iderType,ipoint,iel)
              end do

              ! Save the value in the point
              Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

            end do ! ivar
          end do ! ipoint
        end do ! iel
      end if

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
                  celement, IdofsTrial, npoints, nelements, Dpoints, iderType,&
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
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
!</description>

!<input>
  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those
  ! elements that we are concerning here.
  integer(I32), intent(in) :: celement

  ! A list of the corner vertices of all elements in progress.
  ! array [1..NDIM2D,1..TRIA_MAXNVE2D,1..Number of elements] of double
  real(DP), dimension(:,:,:), intent(in) :: Dcoords

  ! The Jacobian matrix of the mapping between the reference and each
  ! real element, for all points on all elements in progress.
  ! array [1..TRAFO_NJACENTRIES,1..npointsPerElement,1..Number of elements]
  real(DP), dimension(:,:,:),intent(in) :: Djac

  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to each real element in progress.
  ! array [1..npointsPerElement,1..Number of elements]
  real(DP), dimension(:,:), intent(in) :: Ddetj

  ! An array accepting the DOF`s on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF`s in trial space,nelements)
  integer, dimension(:,:), intent(in) :: IdofsTrial

  ! Number of points on every element where to evalate the function
  integer, intent(in) :: npoints

  ! Number of elements, the function is evaluated at
  integer, intent(in)  :: nelements

  ! Array with coordinates of the points where to evaluate.
  ! DIMENSION(NDIM2D,npoints,nelements).
  ! The coordinates are expected
  ! - on the reference element, if celement identifies a parametric element
  ! - on the real element, if celement identifies a nonparametric element
  ! It is assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! OPTIONAL: List of twist indices. Defines for every element the
  ! orgientation of the edges.
  ! Can be omitted if the element does not need it.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  integer(I32), dimension(:), intent(in), optional :: ItwistIndexEdges
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial
  real(DP) :: dval
  integer :: iel,ipoint,ibas,ivar
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

#if WARN_DEPREC
  call output_line ("Using deprecated feature. Please update your code.", &
      OU_CLASS_WARNING,OU_MODE_STD,"fevl_evaluate_sim2")
#endif

  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.

  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))

  ! Evaluate the basis functions
  call elem_generic_sim (celement, Dcoords, Djac, Ddetj, &
      Bder, DbasTrial, npoints, nelements, Dpoints, ItwistIndexEdges)

  if (rvectorScalar%cdataType .eq. ST_DOUBLE) then

    ! Get the data array from the vector
    call lsyssc_getbase_double(rvectorScalar,p_Ddata)

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Ddata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Ddata((IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Fdata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Fdata((IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else
    call output_line("Unsupported vector precision!",&
      OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_sim2")
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
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
!</description>

!<input>
  ! Element evaluation set that contains all information necessary
  ! for the evaluation (coordinates of the points, transformation,...)
  type(t_evalElementSet), intent(in) :: revalElementSet

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those
  ! elements that we are concerning here.
  integer(I32), intent(in) :: celement

  ! An array accepting the DOF`s on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF`s in trial space,nelements)
  integer, dimension(:,:), intent(in) :: IdofsTrial

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  real(DP) :: dval
  integer :: iel,ipoint,ibas,ivar
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata


  npoints = revalElementSet%npointsPerElement
  nelements = revalElementSet%nelements

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

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP
          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Ddata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Ddata((IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Fdata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Fdata((IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else
    call output_line("Unsupported vector precision!",&
      OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_sim3")
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
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
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
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! Number of the subarray in Dvalues where the result is written to.
  ! The routine writes the resulting values to Dvalues(iterm,:,:).
  integer, intent(in) :: iterm
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(#possible terms,npoints,nelements).
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  real(DP) :: dval
  integer :: iel,ipoint,ibas,ivar
  integer(I32) :: celement
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  integer, dimension(:,:), pointer :: p_IdofsTrial


  npoints = rdomainIntSubset%npointsPerElement
  nelements = rdomainIntSubset%nelements

  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.

  ! Get the currently active element
  celement = rdomainIntSubset%celement

  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))

  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

  ! Get the pointer to the trial DOF`s.
  ! If the IdofsTrial in the domain subset fits to our current element,
  ! take that. Otherwise, we have to compute the actual DOF`s.
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

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Ddata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(iterm,ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Ddata((p_IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues(iterm,(ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Fdata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(iterm,ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the FE-coefficients
      ! with the values of the basis functions and summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Fdata((p_IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues(iterm,(ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else
    call output_line("Unsupported vector precision!",&
      OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_sim4")
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

  subroutine fevl_evaluate_sim5 (rvectorScalar, &
                                 rdomainIntSubset, iderType, Dvalues)

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
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
!</description>

!<input>
  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being evaluated.
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  real(DP) :: dval
  integer :: iel,ipoint,ibas,ivar
  integer(I32) :: celement
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  integer, dimension(:,:), pointer :: p_IdofsTrial


  npoints = rdomainIntSubset%npointsPerElement
  nelements = rdomainIntSubset%nelements

  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.

  ! Get the currently active element
  celement = rdomainIntSubset%celement

  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))

  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

  ! Get the pointer to the trial DOF`s.
  ! If the IdofsTrial in the domain subset fits to our current element,
  ! take that. Otherwise, we have to compute the actual DOF`s.
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

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Ddata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar = 1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Ddata((p_IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

    ! Get the data array from the vector
    call lsyssc_getbase_single(rvectorScalar,p_Fdata)

    ! Check if vector is stored in interleaved format
    if (rvectorScalar%NVAR .eq. 1) then

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indofTrial
            dval = dval + &
                   p_Fdata(p_IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
          end do

          ! Save the value in the point
          Dvalues(ipoint,iel) = dval

        end do ! ipoint
      end do ! iel

    else

      ! Now that we have the basis functions, we want to have the
      ! function values.  We get them by multiplying the
      ! FE-coefficients with the values of the basis functions and
      ! summing up.

      do iel=1,nelements
        do ipoint = 1,npoints
          do ivar =1,rvectorScalar%NVAR

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indofTrial
              dval = dval + &
                     p_Fdata((p_IdofsTrial(ibas,iel)-1)*rvectorScalar%NVAR+ivar) *&
                     DbasTrial(ibas,iderType,ipoint,iel)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar,iel) = dval

          end do ! ivar
        end do ! ipoint
      end do ! iel

    end if

  else
    call output_line("Unsupported vector precision!",&
      OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluate_sim5")
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

  subroutine fevl_evaluate_generic_sim1 (iderType, Dvalues, rspatialDiscr,&
      ffunctionCoefficient, Dpoints, Ielements, DpointsRef, rcollection, rperfconfig)

!<description>
  ! This is a rather general finite element evaluation routine. It
  ! allows to evaluate a general (scalar) FE function in a set of
  ! points Dpoints on a set of elements Ielements. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! Spatial discretisation underlying the FE function
  type(t_spatialdiscretisation), intent(in) :: rspatialDiscr

  ! Callback function used to evaluate the FE function
  include "intf_fecoefficient_sim.inc"

  ! A list of points where to evaluate. All points must be inside
  ! of element ielement.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! A list of elements containing the points in Dpoints.
  ! All elements in this list must be of the same type!!!
  integer, dimension(:), intent(in) :: Ielements

  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  real(DP), dimension(:,:,:), intent(in), target, optional :: DpointsRef

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(1..npoints,1..nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    logical :: bnonpar
    integer :: indof,nve,iel,ipoint
    integer(I32) :: celement
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder

    real(DP), dimension(:,:,:), pointer :: p_DpointsRef

    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation

    ! Transformation
    integer(I32) :: ctrafoType

    ! Values of basis functions and DOF`s
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
    p_rtriangulation => rspatialDiscr%p_rtriangulation

    p_RelementDistribution => rspatialDiscr%RelementDistr

    ! For uniform discretisations, we get the element type in advance...
    if (rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then

      ! Element type
      celement = p_RelementDistribution(1)%celement

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Element nonparametric?
      bnonpar = elem_isNonparametric(celement)

      nullify(p_IelementDistr)
    else
      call storage_getbase_int (rspatialDiscr%h_IelementDistr, p_IelementDistr)
    end if

    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.

    ! Get the type of the element ielement
    if (associated(p_IelementDistr)) then
      ! As all elements have the same type, we get the element
      ! characteristics by checking the first element.
      celement = p_RelementDistribution(p_IelementDistr(Ielements(1)))%celement

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Element nonparametric?
      bnonpar = elem_isNonparametric(celement)

    end if

    ! Calculate the global DOF`s on that element into IdofsTest.
    allocate(Idofs(indof,size(Ielements)))
    call dof_locGlobMapping_mult(rspatialDiscr, Ielements, Idofs)

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Get the coordinates of the corners of the elements
    call elprep_prepareSetForEvaluation (revalElementSet,&
        EL_EVLTAG_COORDS, p_rtriangulation, Ielements, ctrafoType,&
        DpointsRef=DpointsRef,DpointsReal=Dpoints,rperfconfig=rperfconfig)

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

    ! Do not create coordinates on the reference/real element; we do this manually!
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REALPOINTS))

    ! Do not calculate element shape information, we have that already.
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_COORDS))

    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.

    ! Prepare the element set for the evaluation
    call elprep_prepareSetForEvaluation (revalElementSet,&
        cevaluationTag, p_rtriangulation, Ielements, ctrafoType,&
        DpointsRef=p_DpointsRef,DpointsReal=Dpoints,rperfconfig=rperfconfig)

    ! Calculate the values of the basis functions in the given points.
    allocate(Dbas(indof,&
             elem_getMaxDerivative(celement),&
             ubound(Dpoints,2), size(Ielements)))
    call elem_generic_sim2 (celement, revalElementSet, Bder, Dbas)

    ! Calculate the desired values by the callback function.
    call ffunctionCoefficient(iderType, size(Ielements), ubound(Dpoints,2), indof,&
                              Idofs, Dpoints, Dbas, Dvalues, rcollection)

    ! Release allocated memory
    ! Remove the reference to DpointsRef again
    deallocate(Dbas)

    call elprep_releaseElementSet(revalElementSet)
    deallocate(Idofs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_generic_sim3 (revalElementSet, ffunctionCoefficient,&
      celement, IdofsTrial, iderType, Dvalues, rcollection)

!<description>
  ! This routine allows to evaluate a finite element solution
  ! simultaneously in multiple points on multiple elements in a
  ! discretisation.
  ! revalElementScalar must specify all information about where and how
  ! to evaluate; e.g. the coordinates of the evaluation points are
  ! to be found here. The routine will then evaluate the element celement
  ! in these points.
  !
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
!</description>

!<input>
  ! Element evaluation set that contains all information necessary
  ! for the evaluation (coordinates of the points, transformation,...)
  type(t_evalElementSet), intent(in) :: revalElementSet

  ! Callback function used to evaluate the FE function
  include "intf_fecoefficient_sim.inc"

  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. celement defines the type
  ! of FE trial function that was used for the discretisation on those
  ! elements that we are concerning here.
  integer(I32), intent(in) :: celement

  ! An array accepting the DOF`s on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF`s in trial space,nelements)
  integer, dimension(:,:), intent(in) :: IdofsTrial

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements


  npoints = revalElementSet%npointsPerElement
  nelements = revalElementSet%nelements

  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.

  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))

  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, revalElementSet, Bder, DbasTrial)

  ! Calculate the desired values by the callback function.
  call ffunctionCoefficient(iderType, nelements, npoints, indofTrial,&
                            IdofsTrial, revalElementSet%p_DpointsRef,&
                            DbasTrial, Dvalues, rcollection)

  ! Release memory, finish
  deallocate(DbasTrial)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_generic_sim4 (rspatialDiscr, ffunctionCoefficient, &
                                         rdomainIntSubset, iderType, Dvalues, &
                                         iterm, rcollection)

!<description>
  ! This routine allows to evaluate a finite element solution
  ! simultaneously in multiple points on multiple elements in a
  ! discretisation.
  ! rdomainIntSubset must specify all information about where and how
  ! to evaluate; e.g. the coordinates of the evaluation points are
  ! to be found here.
  !
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
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
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! Spatial discretisation underlying the FE function
  type(t_spatialdiscretisation), intent(in) :: rspatialDiscr

  ! Callback function used to evaluate the FE function
  include "intf_fecoefficient_sim.inc"

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! Number of the subarray in Dvalues where the result is written to.
  ! The routine writes the resulting values to Dvalues(iterm,:,:).
  integer, intent(in) :: iterm
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(#possible terms,npoints,nelements).
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  integer(I32) :: celement
  integer, dimension(:,:), pointer :: p_IdofsTrial


  npoints = rdomainIntSubset%npointsPerElement
  nelements = rdomainIntSubset%nelements

  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.

  ! Get the currently active element
  celement = rdomainIntSubset%celement

  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))

  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

  ! Get the pointer to the trial DOF`s.
  ! If the IdofsTrial in the domain subset fits to our current element,
  ! take that. Otherwise, we have to compute the actual DOF`s.
  if (rdomainIntSubset%celement .eq. celement) then
    p_IdofsTrial => rdomainIntSubset%p_IdofsTrial
  else
    allocate (p_IdofsTrial(indofTrial,nelements))
    call dof_locGlobMapping_mult(rspatialDiscr, rdomainIntSubset%p_Ielements, p_IdofsTrial)
  end if

  ! Calculate the desired values by the callback function.
  call ffunctionCoefficient(iderType, nelements, npoints, indofTrial,&
                            p_IdofsTrial, rdomainIntSubset%p_DcubPtsReal,&
                            DbasTrial, Dvalues(iterm,:,:), rcollection)

  ! Release memory, finish
  deallocate(DbasTrial)

  if (rdomainIntSubset%celement .ne. celement) then
    deallocate (p_IdofsTrial)
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluate_generic_sim5 (rspatialDiscr, ffunctionCoefficient, &
                                         rdomainIntSubset, iderType, Dvalues, &
                                         rcollection)

!<description>
  ! This routine allows to evaluate a finite element solution
  ! simultaneously in multiple points on multiple elements in a
  ! discretisation.
  ! rdomainIntSubset must specify all information about where and how
  ! to evaluate; e.g. the coordinates of the evaluation points are
  ! to be found here.
  !
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating "simultaneously" on some <tex>$Q_1$</tex> and some <tex>$P_1$</tex>
  ! elements is not allowed e.g..
!</description>

!<input>
  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being evaluated.
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! Spatial discretisation underlying the FE function
  type(t_spatialdiscretisation), intent(in) :: rspatialDiscr

  ! Callback function used to evaluate the FE function
  include "intf_fecoefficient_sim.inc"

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: DbasTrial
  integer :: indofTrial,npoints,nelements
  integer(I32) :: celement
  integer, dimension(:,:), pointer :: p_IdofsTrial


  npoints = rdomainIntSubset%npointsPerElement
  nelements = rdomainIntSubset%nelements

  ! What to evaluate?
  Bder = .false.
  Bder(iderType) = .true.

  ! Get the currently active element
  celement = rdomainIntSubset%celement

  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(celement)
  allocate(DbasTrial(indofTrial,elem_getMaxDerivative(celement),npoints,nelements))

  ! Evaluate the basis functions
  call elem_generic_sim2 (celement, rdomainIntSubset%p_revalElementSet, Bder, DbasTrial)

  ! Get the pointer to the trial DOF`s.
  ! If the IdofsTrial in the domain subset fits to our current element,
  ! take that. Otherwise, we have to compute the actual DOF`s.
  if (rdomainIntSubset%celement .eq. celement) then
    p_IdofsTrial => rdomainIntSubset%p_IdofsTrial
  else
    allocate (p_IdofsTrial(indofTrial,nelements))
    call dof_locGlobMapping_mult(rspatialDiscr, rdomainIntSubset%p_Ielements, p_IdofsTrial)
  end if

  ! Calculate the desired values by the callback function.
  call ffunctionCoefficient(iderType, nelements, npoints, indofTrial,&
                            p_IdofsTrial, rdomainIntSubset%p_DcubPtsReal,&
                            DbasTrial, Dvalues, rcollection)

  ! Release memory, finish
  deallocate(DbasTrial)

  if (rdomainIntSubset%celement .ne. celement) then
    deallocate (p_IdofsTrial)
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr1d1 (iderType, Dvalues, rvectorScalar, ibdc)

!<description>
  ! This is the most general (and completely slowest) finite element
  ! evaluation routine. It allows to evaluate a general (scalar) FE
  ! function specified by rvectorScalar at the unique boundary
  ! point. The values of the FE function are written into Dvalues. The
  ! routine is called via the interface fevl_evaluatBdr1d(...).
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! Number of the boundary component
  integer, intent(in) :: ibdc
!</input>

!<output>
  ! Values of the FE function at the points specified by DpointsPar.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialdiscretisation), pointer :: p_rspatialDiscr
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: NELbdc

    ! The vector must provide a discretisation structure
    if (associated(rvectorScalar%p_rspatialDiscr)) then
      p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
    else
      call output_line("Vector does not provide discretisation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr1d1")
      call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if (associated(p_rspatialDiscr%p_rtriangulation)) then
      p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    else
      call output_line("Discretisation does not provide triangulation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr1d1")
      call sys_halt()
    end if

    ! Get number of elements at boundary component
    NELbdc = bdraux_getNELAtBdrComp(ibdc, p_rtriangulation)

    ! Allocate temporal memory
    allocate(IelementList(NELbdc), IelementOrientation(NELbdc))

    ! Get list of elements at boundary region
    call bdraux_getElementsAtBdrComp(ibdc, p_rspatialDiscr, NELbdc,&
        IelementList, IelementOrientation)

    ! Call FE-evaluation routine
    call fevl_evaluateBdr1d2(iderType, Dvalues, rvectorScalar,&
        IelementList, IelementOrientation)

    ! Deallocate temporal memory
    deallocate(IelementList, IelementOrientation)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr1d2 (iderType, Dvalues, rvectorScalar, &
      IelementList, IelementOrientation)

!<description>
  ! This is the most general (and completely slowest) finite element
  ! evaluation routine. It allows to evaluate a general (scalar) FE
  ! function specified by rvectorScalar in the unique boundary
  ! point. The values of the FE function are written into Dvalues. The
  ! routine is called via the interface fevl_evaluatBdr1d(...).
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A list of elements adjacent to the boundary
  ! containing the points DpointsPar.
  integer, dimension(:), intent(in) :: IelementList

  ! The orientation of elements adjacent to the boundary region
  integer, dimension(:), intent(in) :: IelementOrientation
!</input>

!<output>
  ! Values of the FE function at the points specified by DpointsPar.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: icoordSystem
    integer :: indof,nve,ibas
    integer(I32) :: celement
    integer :: iel,idx,ivar
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP), dimension(NDIM1D+2) :: DpointRef
    real(DP) :: dval

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! Transformation
    integer(I32) :: ctrafoType

    ! Values of basis functions and DOF`s
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
      celement = rvectorScalar%p_rspatialDiscr%RelementDistr(1)%celement

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Get the element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)

      ! Type of coordinate system
      icoordSystem = elem_igetCoordSystem(celement)

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
      call output_line ("Unsupported vector precision!",&
          OU_CLASS_ERROR, OU_MODE_STD, "fevl_evaluateBdr1d2")
      call sys_halt()
    end select

    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.

    ! We loop over all elements.
    do idx = 1, size(IelementList)

      ! Get element number
      iel = IelementList(idx)

      ! Get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! Get the number of local DOF`s for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! Number of vertices on the element
        nve = elem_igetNVE(celement)

        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)

        ! Type of coordinate system
        icoordSystem = elem_igetCoordSystem(celement)
      end if

      ! Calculate the global DOF`s on that element into IdofsTest.
      call dof_locGlobMapping (rvectorScalar%p_rspatialDiscr, iel, Idofs)

      ! Compute coordinate of cubature point on reference element
      if (IelementOrientation(idx) .eq. 1) then
        DpointRef = -1.0
      else
        DpointRef = 1.0
      end if

      ! Now calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
          cevaluationTag, rvectorScalar%p_rspatialDiscr%p_rtriangulation,&
          iel, ctrafoType, DpointRef=DpointRef)

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! Combine the basis functions to get the function value.
      if (rvectorScalar%cdataType .eq. ST_DOUBLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the
          ! function values.  We get them by multiplying the
          ! FE-coefficients with the values of the basis functions and
          ! summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(idx) = dval

        else

          do ivar = 1, rvectorScalar%NVAR
            ! Now that we have the basis functions, we want to have the
            ! function values.  We get them by multiplying the
            ! FE-coefficients with the values of the basis functions and
            ! summing up.

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Ddata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((idx-1)*rvectorScalar%NVAR+ivar) = dval
          end do

        end if

      else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the
          ! function values.  We get them by multiplying the
          ! FE-coefficients with the values of the basis functions and
          ! summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(idx) = dval

        else

          do ivar = 1, rvectorScalar%NVAR
            ! Now that we have the basis functions, we want to have the
            ! function values.  We get them by multiplying the
            ! FE-coefficients with the values of the basis functions and
            ! summing up.

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Fdata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((idx-1)*rvectorScalar%NVAR+ivar) = dval
          end do

        end if

      end if

    end do ! idx

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr1d3 (CderType, Dvalues, rvectorBlock, ibdc,&
      iblockMin, iblockMax)

!<description>
  ! This subroutine is an extension of fevl_evaluateBdr1d1(). It
  ! simultaneously evaluates several types of function values (given
  ! by CderType) for all components of a multivariate FE function
  ! (given by rvectorBlock) in the unique boundary point. Thus, the
  ! evaluation of all the element basis functions etc. has to be
  ! performed only once. The routine is called via the interface
  ! fevl_evaluateBdr1D(...).
  ! E.g., instead of,
  !
  !   real(DP), dimension(nblocks,npoints) :: Dvalues
  !   call fevl_evaluateBdr1d(DER_FUNC, Dvalues(1,:), rsol%RvectorBlock(1), ibdc)
  !   call fevl_evaluateBdr1d(DER_FUNC, Dvalues(2,:), rsol%RvectorBlock(2), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_X, DderivX(1,:), rsol%RvectorBlock(1), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_X, DderivX(2,:), rsol%RvectorBlock(2), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_Y, DderivY(1,:), rsol%RvectorBlock(1), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_Y, DderivY(2,:), rsol%RvectorBlock(2), ibdc)
  !
  ! one can use
  !
  !   integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)
  !   real(DP), dimension(nblocks,3,npoints) :: Dvalues
  !   call fevl_evaluateBdr1d(CderType, Dvalues, rsol, ibdc)
  !
!</description>

!<input>
  ! array of type of function values to evaluate (DER_FUNC, DER_DERIV_X etc.)
  integer, dimension(:), intent(in) :: CderType

  ! block solution vector representing the FE function that is to be evaluated
  ! It is assumed that all components correspond to the same spatial discretisation and
  ! use the same element type (i.e., only that of the first component is inquired)
  type(t_vectorBlock), intent(in) :: rvectorBlock

  ! Number of the boundary component
  integer, intent(in) :: ibdc

  ! OPTIONAL: For the case that not all components of the vector are to be processed, the
  ! user can provide these two parameters so that only the components from iblockMin to
  ! iblockMax are processed. This can be necessary, e.g., in case of a saddle point
  ! problem where the last component of the solution vector (=pressure) is discretised
  ! in a different way than the first components (=velocity/displacement).
  ! If iblockMin (iblockMax) is not present, the minimal (maximal) block number is
  ! assumed to be 1 (rvectorBlock%nblocks).
  integer, intent(in), optional :: iblockMin, iblockMax
!</input>

!<output>
  ! values of the FE function at the points specified by Dpoints.
  ! DIMENSION(rvectorBlock%nblocks, size(CderType), npoints)
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialdiscretisation), pointer :: p_rspatialDiscr
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: NELbdc,iblMin

    if (present(iblockMin)) then
      iblMin = iblockMin
    else
      iblMin = 1
    endif

    ! The vector must provide a discretisation structure
    if (associated(rvectorBlock%p_rblockDiscr%RspatialDiscr)) then
      p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(iblMin)
    else
      call output_line("Vector does not provide discretisation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr1d3")
      call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if (associated(p_rspatialDiscr%p_rtriangulation)) then
      p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    else
      call output_line("Discretisation does not provide triangulation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr1d3")
      call sys_halt()
    end if

    ! Get number of elements at boundary component
    NELbdc = bdraux_getNELAtBdrComp(ibdc, p_rtriangulation)

    ! Allocate temporal memory
    allocate(IelementList(NELbdc), IelementOrientation(NELbdc))

    ! Get list of elements at boundary region
    call bdraux_getElementsAtBdrComp(ibdc, p_rspatialDiscr, NELbdc,&
        IelementList, IelementOrientation)

    ! Call FE-evaluation routine
    call fevl_evaluateBdr1d4(CderType, Dvalues, rvectorBlock,&
        IelementList, IelementOrientation, iblockMin, iblockMax)

    ! Deallocate temporal memory
    deallocate(IelementList, IelementOrientation)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr1d4 (CderType, Dvalues, rvectorBlock,&
      IelementList, IelementOrientation, iblockMin, iblockMax)

!<description>
  ! This is the most general (and completely slowest) finite element
  ! evaluation routine. It allows to evaluate a general (scalar) FE
  ! function specified by rvectorScalar in the unique boundary
  ! point. The values of the FE function are written into Dvalues. The
  ! routine is called via the interface fevl_evaluatBdr1d(...).
!</description>

!<description>

  ! This subroutine is an extension of fevl_evaluateBdr2d1(). It
  ! simultaneously evaluates several types of function values (given
  ! by CderType) for all components of a multivariate FE function
  ! (given by rvectorBlock) in a set of points on the boundary
  ! specified by its parameter values DpointsPar. Thus, the evaluation
  ! of all the element basis functions etc. has to be performed only
  ! once. The routine is called via the interface fevl_evaluate(...).
  ! E.g., instead of,
  !
  !   real(DP), dimension(nblocks,npoints) :: Dvalues
  !   call fevl_evaluateBdr1d(DER_FUNC, Dvalues(1,:), rsol%RvectorBlock(1), ibdc)
  !   call fevl_evaluateBdr1d(DER_FUNC, Dvalues(2,:), rsol%RvectorBlock(2), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_X, DderivX(1,:), rsol%RvectorBlock(1), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_X, DderivX(2,:), rsol%RvectorBlock(2), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_Y, DderivY(1,:), rsol%RvectorBlock(1), ibdc)
  !   call fevl_evaluateBdr1d(DER_DERIV_Y, DderivY(2,:), rsol%RvectorBlock(2), ibdc)
  !
  ! one can use
  !
  !   integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)
  !   real(DP), dimension(nblocks,3,npoints) :: Dvalues
  !   call fevl_evaluateBdr1d(CderType, Dvalues, rsol, ibdc)
  !
!</description>

!<input>
  ! array of type of function values to evaluate (DER_FUNC, DER_DERIV_X etc.)
  integer, dimension(:), intent(in) :: CderType

  ! block solution vector representing the FE function that is to be evaluated
  ! It is assumed that all components correspond to the same spatial discretisation and
  ! use the same element type (i.e., only that of the first component is inquired)
  type(t_vectorBlock), intent(in) :: rvectorBlock

  ! A list of elements adjacent to the boundary
  ! containing the points DpointsPar.
  integer, dimension(:), intent(in) :: IelementList

  ! The orientation of elements adjacent to the boundary region
  integer, dimension(:), intent(in) :: IelementOrientation

  ! OPTIONAL: For the case that not all components of the vector are to be processed, the
  ! user can provide these two parameters so that only the components from iblockMin to
  ! iblockMax are processed. This can be necessary, e.g., in case of a saddle point
  ! problem where the last component of the solution vector (=pressure) is discretised
  ! in a different way than the first components (=velocity/displacement).
  ! If iblockMin (iblockMax) is not present, the minimal (maximal) block number is
  ! assumed to be 1 (rvectorBlock%nblocks).
  integer, intent(in), optional :: iblockMin, iblockMax
!</input>

!<output>
  ! values of the FE function at the points specified by Dpoints.
  ! DIMENSION(rvectorBlock%nblocks, size(CderType), npoints)
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: icoordSystem
    integer :: ipoint,indof,nve,ibas,ider,iblock,iblMin,iblMax
    integer(I32) :: celement
    integer :: iel,idx
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP), dimension(NDIM1D+2) :: DpointRef
    real(DP) :: dval

    ! number of equations in one scalar component of the block vector
    integer :: neqsc

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! Transformation
    integer(I32) :: ctrafoType

    ! Values of basis functions and DOF`s
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    integer, dimension(EL_MAXNBAS) :: Idofs

    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! pointer to the spatial discretisation structure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    ! Evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag

    if (present(iblockMin)) then
      iblMin = iblockMin
    else
      iblMin = 1
    endif
    if (present(iblockMax)) then
      iblMax = iblockMax
    else
      iblMax = rvectorBlock%nblocks
    endif

    ! Ok, slow but general.

    p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(iblMin)
    p_RelementDistribution => p_rspatialDiscr%RelementDistr

    ! For uniform discretisations, we get the element type in advance...
    if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then

      ! Element type
      celement = p_rspatialDiscr%RelementDistr(1)%celement

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Get the element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)

      ! Type of coordinate system
      icoordSystem = elem_igetCoordSystem(celement)

      nullify(p_IelementDistr)
    else
      call storage_getbase_int (p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
    end if

    ! Get the data vector
    select case (rvectorBlock%cdataType)
    case (ST_DOUBLE)
      call lsysbl_getbase_double(rvectorBlock,p_Ddata)
    case (ST_SINGLE)
      call lsysbl_getbase_single(rvectorBlock,p_Fdata)
    case DEFAULT
      call output_line ("Unsupported vector precision!",&
          OU_CLASS_ERROR, OU_MODE_STD, "fevl_evaluateBdr1d4")
      call sys_halt()
    end select

    ! inquire which types of function values are to be evaluated
    Bder = .false.
    do ider = 1, size(CderType)
      Bder(CderType(ider)) = .true.
    end do

    ! We loop over all elements.
    do idx = 1, size(IelementList)

      ! Get element number
      iel = IelementList(idx)

      ! Get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! Get the number of local DOF`s for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! Number of vertices on the element
        nve = elem_igetNVE(celement)

        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)

        ! Type of coordinate system
        icoordSystem = elem_igetCoordSystem(celement)
      end if

      ! Calculate the global DOF`s on that element into IdofsTest.
      call dof_locGlobMapping (p_rspatialDiscr, iel, Idofs)

      ! Compute coordinate of cubature point on reference element
      if (IelementOrientation(idx) .eq. 1) then
        DpointRef = -1.0
      else
        DpointRef = 1.0
      end if

      ! Now calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
          cevaluationTag, p_rspatialDiscr%p_rtriangulation,&
          iel, ctrafoType, DpointRef=DpointRef)

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! number of equations in one scalar component
      neqsc = rvectorBlock%RvectorBlock(iblMin)%neq

      ! calculate function values by multiplying the FE-coefficients with the values of
      ! the basis functions and summing up
      Dvalues(iblMin:iblMax,:,ipoint) = 0.0_DP
      if (rvectorBlock%cdataType .eq. ST_DOUBLE) then
        do ider = 1,size(CderType)
          do iblock = iblMin,iblMax
            dval = 0.0_DP
            do ibas = 1,indof
              dval = dval +   p_Ddata((iblock-1)*neqsc + Idofs(ibas)) &
                            * Dbas(ibas,CderType(ider))
            end do
            Dvalues(iblock, ider, ipoint) = dval
          end do
        end do
      else if (rvectorBlock%cdataType .eq. ST_SINGLE) then
        do ider = 1,size(CderType)
          do iblock = iblMin,iblMax
            dval = 0.0_DP
            do ibas = 1,indof
              dval = dval +   p_Fdata((iblock-1)*neqsc + Idofs(ibas)) &
                            * Dbas(ibas,CderType(ider))
            end do
            Dvalues(iblock, ider, ipoint) = dval
          end do
        end do
      end if

    end do ! idx

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr2d1 (iderType, Dvalues, rvectorScalar, DpointsPar,&
      ibdc, cparType, rboundaryRegion)

!<description>
  ! This is the most general (and completely slowest) finite element
  ! evaluation routine. It allows to evaluate a general (scalar) FE
  ! function specified by rvectorScalar in a set of points on the
  ! boundary specified by its parameter values DpointsPar. The values
  ! of the FE function are written into Dvalues. The routine is called
  ! via the interface fevl_evaluatBdr2d(...).
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A list of parametrised points where to evaluate.
  ! DIMENSION(1..npoints)
  real(DP), dimension(:), intent(in) :: DpointsPar

  ! Number of the boundary component
  integer, intent(in) :: ibdc

  ! Type of parametrisation used for array DpointsPar
  integer, intent(in) :: cparType

  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to search for the elements (must be given in 0-1 parametrisation)
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion
!</input>

!<output>
  ! Values of the FE function at the points specified by DpointsPar.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialdiscretisation), pointer :: p_rspatialDiscr
    real(DP), dimension(:,:), pointer :: DedgePosition
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: NELbdc

    ! The vector must provide a discretisation structure
    if (associated(rvectorScalar%p_rspatialDiscr)) then
      p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
    else
      call output_line("Vector does not provide discretisation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d1")
      call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if (associated(p_rspatialDiscr%p_rtriangulation)) then
      p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    else
      call output_line("Discretisation does not provide triangulation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d1")
      call sys_halt()
    end if

    ! Do we have a boundary region?
    if (present(rboundaryRegion)) then
      ! Check if boundary region is given in 0-1 parametrisation
      if (rboundaryRegion%cparType .ne. BDR_PAR_01) then
        call output_line ("Boundary region must be given in 0-1 parametrisation!", &
            OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d1")
        call sys_halt()
      end if

      ! Get number of elements at boundary region
      NELbdc = bdraux_getNELAtRegion(rboundaryRegion, p_rtriangulation)

      ! Allocate temporal memory
      allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
      allocate(DedgePosition(2,NELbdc))

      ! Get list of elements at boundary region
      call bdraux_getElementsAtRegion(rboundaryRegion, p_rspatialDiscr,&
          NELbdc, IelementList, IelementOrientation, DedgePosition,&
          cparType=cparType)

      ! Call FE-evaluation routine
      call fevl_evaluateBdr2d2(iderType, Dvalues, rvectorScalar,&
          DpointsPar, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
          DedgePosition(:,1:NELbdc))

      ! Deallocate temporal memory
      deallocate(IelementList, IelementOrientation, DedgePosition)

    else
      ! Get number of elements at boundary component
      NELbdc = bdraux_getNELAtBdrComp(ibdc, p_rtriangulation)

      ! Allocate temporal memory
      allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
      allocate(DedgePosition(2,NELbdc))

      ! Get list of elements at boundary region
      call bdraux_getElementsAtBdrComp(ibdc, p_rspatialDiscr, NELbdc,&
          IelementList, IelementOrientation, DedgePosition, cparType=cparType)

      ! Call FE-evaluation routine
      call fevl_evaluateBdr2d2(iderType, Dvalues, rvectorScalar, DpointsPar,&
          IelementList, IelementOrientation, DedgePosition)

      ! Deallocate temporal memory
      deallocate(IelementList, IelementOrientation, DedgePosition)
    end if

  end subroutine fevl_evaluateBdr2d1

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr2d2 (iderType, Dvalues, rvectorScalar, DpointsPar,&
      IelementList, IelementOrientation, DedgePosition)

!<description>
  ! This is the most general (and completely slowest) finite element
  ! evaluation routine. It allows to evaluate a general (scalar) FE
  ! function specified by rvectorScalar in a set of points on the
  ! boundary specified by its parameter values DpointsPar. The values
  ! of the FE function are written into Dvalues. The routine is called
  ! via the interface fevl_evaluatBdr2d(...).
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  integer, intent(in) :: iderType

  ! The scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A list of parametrised points where to evaluate.
  ! DIMENSION(1..npoints)
  real(DP), dimension(:), intent(in) :: DpointsPar

  ! A list of elements adjacent to the boundary
  ! containing the points DpointsPar.
  integer, dimension(:), intent(in) :: IelementList

  ! The orientation of elements adjacent to the boundary region
  integer, dimension(:), intent(in) :: IelementOrientation

  ! The start- and end-parameter values of the edges on the boundary.
  ! The type of parametrisation of DpointsPar and DedgePosition must
  ! be the same which is not checked by this routine.
  real(DP), dimension(:,:), intent(in) :: DedgePosition
!</input>

!<output>
  ! Values of the FE function at the points specified by DpointsPar.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer(I32) :: icoordSystem
    integer :: ipoint,indof,nve,ibas
    integer(I32) :: celement
    integer :: iel,idx,idxFirst,idxLast,ivar
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP), dimension(1,NDIM2D+1) :: Dxi2d
    real(DP), dimension(NDIM2D+1) :: DpointsRef
    real(DP), dimension(1,1) :: Dxi1D
    real(DP) :: dval

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! Transformation
    integer(I32) :: ctrafoType

    ! Values of basis functions and DOF`s
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
      celement = rvectorScalar%p_rspatialDiscr%RelementDistr(1)%celement

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Get the element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)

      ! Type of coordinate system
      icoordSystem = elem_igetCoordSystem(celement)

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
      call output_line ("Unsupported vector precision!",&
          OU_CLASS_ERROR, OU_MODE_STD, "fevl_evaluateBdr2d2")
      call sys_halt()
    end select

    ! What to evaluate?
    Bder = .false.
    Bder(iderType) = .true.

    ! We loop over all points.
    do ipoint = 1,size(DpointsPar)

      iel = 0

      ! We have to find the element using binary search
      if (iel .le. 0) then
        idxFirst = 1; idxLast = size(IelementList)

        do while ((idxFirst .le. idxLast) .and. (iel .le. 0))
          idx = idxFirst + ((idxLast-idxFirst)/2)

          if (DpointsPar(ipoint) .lt. DedgePosition(1,idx)) then
            ! Continue binary seach in first part
            idxLast = idx-1
          elseif (DpointsPar(ipoint) .gt. DedgePosition(2,idx)) then
            ! Continue binary seach in second part
            idxFirst = idx+1
          else
            iel = IelementList(idx)
          end if
        end do
      end if

      if (iel .le. 0) then
        call output_line ("Point "//trim(sys_siL(ipoint,10))//" not found!", &
            OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d2")
        cycle
      end if

      ! Get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! Get the number of local DOF`s for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! Number of vertices on the element
        nve = elem_igetNVE(celement)

        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)

        ! Type of coordinate system
        icoordSystem = elem_igetCoordSystem(celement)
      end if

      ! Calculate the global DOF`s on that element into IdofsTest.
      call dof_locGlobMapping (rvectorScalar%p_rspatialDiscr, iel, Idofs)

      ! Transform the parameter value of the boundary point into
      ! reference coordinates along the boundary edge.
      call mprim_linearRescale(DpointsPar(ipoint), DedgePosition(1,idx),&
          DedgePosition(2,idx), -1.0_DP, 1.0_DP, Dxi1D)

      ! Map the 1D cubature point to the edge in 2D.
      call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(idx), &
          1, Dxi1D, Dxi2d)
      DpointsRef = Dxi2d(1,:)

      ! Now calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
          cevaluationTag, rvectorScalar%p_rspatialDiscr%p_rtriangulation,&
          iel, ctrafoType, DpointRef=DpointsRef)

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! Combine the basis functions to get the function value.
      if (rvectorScalar%cdataType .eq. ST_DOUBLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the
          ! function values.  We get them by multiplying the
          ! FE-coefficients with the values of the basis functions and
          ! summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(ipoint) = dval

        else

          do ivar = 1, rvectorScalar%NVAR
            ! Now that we have the basis functions, we want to have the
            ! function values.  We get them by multiplying the
            ! FE-coefficients with the values of the basis functions and
            ! summing up.

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Ddata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval
          end do

        end if

      else if (rvectorScalar%cdataType .eq. ST_SINGLE) then

        ! Check if vector is stored in interleaved format
        if (rvectorScalar%NVAR .eq. 1) then

          ! Now that we have the basis functions, we want to have the
          ! function values.  We get them by multiplying the
          ! FE-coefficients with the values of the basis functions and
          ! summing up.

          dval = 0.0_DP

          ! Calculate the value in the point
          do ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
          end do

          ! Save the value in the point
          Dvalues(ipoint) = dval

        else

          do ivar = 1, rvectorScalar%NVAR
            ! Now that we have the basis functions, we want to have the
            ! function values.  We get them by multiplying the
            ! FE-coefficients with the values of the basis functions and
            ! summing up.

            dval = 0.0_DP

            ! Calculate the value in the point
            do ibas = 1,indof
              dval = dval + p_Fdata((Idofs(ibas)-1)*rvectorScalar%NVAR+ivar) *&
                            Dbas(ibas,iderType)
            end do

            ! Save the value in the point
            Dvalues((ipoint-1)*rvectorScalar%NVAR+ivar) = dval
          end do

        end if

      end if

    end do ! ipoint

  end subroutine fevl_evaluateBdr2d2

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr2d3 (CderType, Dvalues, rvectorBlock, DpointsPar,&
      ibdc, cparType, rboundaryRegion, iblockMin, iblockMax)

!<description>

  ! This subroutine is an extension of fevl_evaluateBdr2d1(). It
  ! simultaneously evaluates several types of function values (given
  ! by CderType) for all components of a multivariate FE function
  ! (given by rvectorBlock) in a set of points on the boundary
  ! specified by its parameter values DpointsPar. Thus, the evaluation
  ! of all the element basis functions etc. has to be performed only
  ! once. The routine is called via the interface fevl_evaluateBdr2D(...).
  ! E.g., instead of,
  !
  !   real(DP), dimension(nblocks,npoints) :: Dvalues
  !   call fevl_evaluateBdr2d(DER_FUNC, Dvalues(1,:), rsol%RvectorBlock(1), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_FUNC, Dvalues(2,:), rsol%RvectorBlock(2), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_X, DderivX(1,:), rsol%RvectorBlock(1), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_X, DderivX(2,:), rsol%RvectorBlock(2), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_Y, DderivY(1,:), rsol%RvectorBlock(1), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_Y, DderivY(2,:), rsol%RvectorBlock(2), DpointsPar, ibdc, cparType)
  !
  ! one can use
  !
  !   integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)
  !   real(DP), dimension(nblocks,3,npoints) :: Dvalues
  !   call fevl_evaluateBdr2d(CderType, Dvalues, rsol, DpointsPar, ibdc, cparType)
  !
!</description>

!<input>
  ! array of type of function values to evaluate (DER_FUNC, DER_DERIV_X etc.)
  integer, dimension(:), intent(in) :: CderType

  ! block solution vector representing the FE function that is to be evaluated
  ! It is assumed that all components correspond to the same spatial discretisation and
  ! use the same element type (i.e., only that of the first component is inquired)
  type(t_vectorBlock), intent(in) :: rvectorBlock

  ! A list of parametrised points where to evaluate.
  ! DIMENSION(1..npoints)
  real(DP), dimension(:), intent(in) :: DpointsPar

  ! Number of the boundary component
  integer, intent(in) :: ibdc

  ! Type of parametrisation used for array DpointsPar
  integer, intent(in) :: cparType

  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to search for the elements (must be given in 0-1 parametrisation)
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion

  ! OPTIONAL: For the case that not all components of the vector are to be processed, the
  ! user can provide these two parameters so that only the components from iblockMin to
  ! iblockMax are processed. This can be necessary, e.g., in case of a saddle point
  ! problem where the last component of the solution vector (=pressure) is discretised
  ! in a different way than the first components (=velocity/displacement).
  ! If iblockMin (iblockMax) is not present, the minimal (maximal) block number is
  ! assumed to be 1 (rvectorBlock%nblocks).
  integer, intent(in), optional :: iblockMin, iblockMax
!</input>

!<output>
  ! values of the FE function at the points specified by Dpoints.
  ! DIMENSION(rvectorBlock%nblocks, size(CderType), npoints)
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialdiscretisation), pointer :: p_rspatialDiscr
    real(DP), dimension(:,:), pointer :: DedgePosition
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: NELbdc,iblMin

    if (present(iblockMin)) then
      iblMin = iblockMin
    else
      iblMin = 1
    endif

    ! The vector must provide a discretisation structure
    if (associated(rvectorBlock%p_rblockDiscr%RspatialDiscr)) then
      p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(iblMin)
    else
      call output_line("Vector does not provide discretisation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d3")
      call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if (associated(p_rspatialDiscr%p_rtriangulation)) then
      p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    else
      call output_line("Discretisation does not provide triangulation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d3")
      call sys_halt()
    end if

    ! Do we have a boundary region?
    if (present(rboundaryRegion)) then
      ! Check if boundary region is given in 0-1 parametrisation
      if (rboundaryRegion%cparType .ne. BDR_PAR_01) then
        call output_line ("Boundary region must be given in 0-1 parametrisation!", &
            OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d3")
        call sys_halt()
      end if

      ! Get number of elements at boundary region
      NELbdc = bdraux_getNELAtRegion(rboundaryRegion, p_rtriangulation)

      ! Allocate temporal memory
      allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
      allocate(DedgePosition(2,NELbdc))

      ! Get list of elements at boundary region
      call bdraux_getElementsAtRegion(rboundaryRegion, p_rspatialDiscr,&
          NELbdc, IelementList, IelementOrientation, DedgePosition,&
          cparType=cparType)

      ! Call FE-evaluation routine
      call fevl_evaluateBdr2d4(CderType, Dvalues, rvectorBlock,&
          DpointsPar, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
          DedgePosition(:,1:NELbdc), iblockMin, iblockMax)

      ! Deallocate temporal memory
      deallocate(IelementList, IelementOrientation, DedgePosition)

    else
      ! Get number of elements at boundary component
      NELbdc = bdraux_getNELAtBdrComp(ibdc, p_rtriangulation)

      ! Allocate temporal memory
      allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
      allocate(DedgePosition(2,NELbdc))

      ! Get list of elements at boundary region
      call bdraux_getElementsAtBdrComp(ibdc, p_rspatialDiscr, NELbdc,&
          IelementList, IelementOrientation, DedgePosition, cparType=cparType)

      ! Call FE-evaluation routine
      call fevl_evaluateBdr2d4(CderType, Dvalues, rvectorBlock, DpointsPar,&
          IelementList, IelementOrientation, DedgePosition, iblockMin, iblockMax)

      ! Deallocate temporal memory
      deallocate(IelementList, IelementOrientation, DedgePosition)
    end if

  end subroutine fevl_evaluateBdr2d3

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_evaluateBdr2d4 (CderType, Dvalues, rvectorBlock, DpointsPar,&
      IelementList, IelementOrientation, DedgePosition, iblockMin, iblockMax)

!<description>
  ! This is the most general (and completely slowest) finite element
  ! evaluation routine. It allows to evaluate a general (scalar) FE
  ! function specified by rvectorScalar in a set of points on the
  ! boundary specified by its parameter values DpointsPar. The values
  ! of the FE function are written into Dvalues. The routine is called
  ! via the interface fevl_evaluatBdr2d(...).
!</description>

!<description>

  ! This subroutine is an extension of fevl_evaluateBdr2d1(). It
  ! simultaneously evaluates several types of function values (given
  ! by CderType) for all components of a multivariate FE function
  ! (given by rvectorBlock) in a set of points on the boundary
  ! specified by its parameter values DpointsPar. Thus, the evaluation
  ! of all the element basis functions etc. has to be performed only
  ! once. The routine is called via the interface fevl_evaluate(...).
  ! E.g., instead of,
  !
  !   real(DP), dimension(nblocks,npoints) :: Dvalues
  !   call fevl_evaluateBdr2d(DER_FUNC, Dvalues(1,:), rsol%RvectorBlock(1), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_FUNC, Dvalues(2,:), rsol%RvectorBlock(2), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_X, DderivX(1,:), rsol%RvectorBlock(1), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_X, DderivX(2,:), rsol%RvectorBlock(2), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_Y, DderivY(1,:), rsol%RvectorBlock(1), DpointsPar, ibdc, cparType)
  !   call fevl_evaluateBdr2d(DER_DERIV_Y, DderivY(2,:), rsol%RvectorBlock(2), DpointsPar, ibdc, cparType)
  !
  ! one can use
  !
  !   integer, dimension(3) :: CderType = (/DER_FUNC, DER_DERIV_X, DER_DERIV_Y/)
  !   real(DP), dimension(nblocks,3,npoints) :: Dvalues
  !   call fevl_evaluateBdr2d(CderType, Dvalues, rsol, DpointsPar, ibdc, cparType)
  !
!</description>

!<input>
  ! array of type of function values to evaluate (DER_FUNC, DER_DERIV_X etc.)
  integer, dimension(:), intent(in) :: CderType

  ! block solution vector representing the FE function that is to be evaluated
  ! It is assumed that all components correspond to the same spatial discretisation and
  ! use the same element type (i.e., only that of the first component is inquired)
  type(t_vectorBlock), intent(in) :: rvectorBlock

  ! A list of parametrised points where to evaluate.
  ! DIMENSION(1..npoints)
  real(DP), dimension(:), intent(in) :: DpointsPar

  ! A list of elements adjacent to the boundary
  ! containing the points DpointsPar.
  integer, dimension(:), intent(in) :: IelementList

  ! The orientation of elements adjacent to the boundary region
  integer, dimension(:), intent(in) :: IelementOrientation

  ! The start- and end-parameter values of the edges on the boundary.
  ! The type of parametrisation of DpointsPar and DedgePosition must
  ! be the same which is not checked by this routine.
  real(DP), dimension(:,:), intent(in) :: DedgePosition

  ! OPTIONAL: For the case that not all components of the vector are to be processed, the
  ! user can provide these two parameters so that only the components from iblockMin to
  ! iblockMax are processed. This can be necessary, e.g., in case of a saddle point
  ! problem where the last component of the solution vector (=pressure) is discretised
  ! in a different way than the first components (=velocity/displacement).
  ! If iblockMin (iblockMax) is not present, the minimal (maximal) block number is
  ! assumed to be 1 (rvectorBlock%nblocks).
  integer, intent(in), optional :: iblockMin, iblockMax
!</input>

!<output>
  ! values of the FE function at the points specified by Dpoints.
  ! DIMENSION(rvectorBlock%nblocks, size(CderType), npoints)
  real(DP), dimension(:,:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer(I32) :: icoordSystem
    integer :: ipoint,indof,nve,ibas,ider,iblock,iblMin,iblMax
    integer(I32) :: celement
    integer :: iel,idx,idxFirst,idxLast
    integer, dimension(:), pointer :: p_IelementDistr
    logical, dimension(EL_MAXNDER) :: Bder
    real(DP), dimension(1,NDIM2D+1) :: Dxi2d
    real(DP), dimension(NDIM2D+1) :: DpointsRef
    real(DP), dimension(1,1) :: Dxi1D
    real(DP) :: dval

    ! number of equations in one scalar component of the block vector
    integer :: neqsc

    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata

    ! Transformation
    integer(I32) :: ctrafoType

    ! Values of basis functions and DOF`s
    real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    integer, dimension(EL_MAXNBAS) :: Idofs

    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! pointer to the spatial discretisation structure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    ! Evaluation structure and tag
    type(t_evalElement) :: revalElement
    integer(I32) :: cevaluationTag

    if (present(iblockMin)) then
      iblMin = iblockMin
    else
      iblMin = 1
    endif
    if (present(iblockMax)) then
      iblMax = iblockMax
    else
      iblMax = rvectorBlock%nblocks
    endif

    ! Ok, slow but general.

    p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(iblMin)
    p_RelementDistribution => p_rspatialDiscr%RelementDistr

    ! For uniform discretisations, we get the element type in advance...
    if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then

      ! Element type
      celement = p_rspatialDiscr%RelementDistr(1)%celement

      ! Get the number of local DOF`s for trial and test functions
      indof = elem_igetNDofLoc(celement)

      ! Number of vertices on the element
      nve = elem_igetNVE(celement)

      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(celement)

      ! Get the element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(celement)

      ! Type of coordinate system
      icoordSystem = elem_igetCoordSystem(celement)

      nullify(p_IelementDistr)
    else
      call storage_getbase_int (p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
    end if

    ! Get the data vector
    select case (rvectorBlock%cdataType)
    case (ST_DOUBLE)
      call lsysbl_getbase_double(rvectorBlock,p_Ddata)
    case (ST_SINGLE)
      call lsysbl_getbase_single(rvectorBlock,p_Fdata)
    case DEFAULT
      call output_line ("Unsupported vector precision!",&
          OU_CLASS_ERROR, OU_MODE_STD, "fevl_evaluateBdr2d4")
      call sys_halt()
    end select

    ! inquire which types of function values are to be evaluated
    Bder = .false.
    do ider = 1, size(CderType)
      Bder(CderType(ider)) = .true.
    end do

    ! We loop over all points.

    do ipoint = 1,size(DpointsPar)

      iel = 0

      ! We have to find the element using binary search
      if (iel .le. 0) then
        idxFirst = 1; idxLast = size(IelementList)

        do while ((idxFirst .le. idxLast) .and. (iel .le. 0))
          idx = idxFirst + ((idxLast-idxFirst)/2)

          if (DpointsPar(ipoint) .lt. DedgePosition(1,idx)) then
            ! Continue binary seach in first part
            idxLast = idx-1
          elseif (DpointsPar(ipoint) .gt. DedgePosition(2,idx)) then
            ! Continue binary seach in second part
            idxFirst = idx+1
          else
            iel = IelementList(idx)
          end if
        end do
      end if

      if (iel .le. 0) then
        call output_line ("Point "//trim(sys_siL(ipoint,10))//" not found!", &
            OU_CLASS_ERROR,OU_MODE_STD,"fevl_evaluateBdr2d3")
        cycle
      end if

      ! Get the type of the element iel
      if (associated(p_IelementDistr)) then
        celement = p_RelementDistribution(p_IelementDistr(iel))%celement

        ! Get the number of local DOF`s for trial and test functions
        indof = elem_igetNDofLoc(celement)

        ! Number of vertices on the element
        nve = elem_igetNVE(celement)

        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(celement)

        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(celement)

        ! Type of coordinate system
        icoordSystem = elem_igetCoordSystem(celement)
      end if

      ! Calculate the global DOF`s on that element into IdofsTest.
      call dof_locGlobMapping (p_rspatialDiscr, iel, Idofs)

      ! Transform the parameter value of the boundary point into
      ! reference coordinates along the boundary edge.
      call mprim_linearRescale(DpointsPar(ipoint), DedgePosition(1,idx),&
          DedgePosition(2,idx), -1.0_DP, 1.0_DP, Dxi1D)

      ! Map the 1D cubature point to the edge in 2D.
      call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(idx), &
          1, Dxi1D, Dxi2d); DpointsRef = Dxi2d(1,:)

      ! Now calculate everything else what is necessary for the element
      call elprep_prepareForEvaluation (revalElement, &
          cevaluationTag, p_rspatialDiscr%p_rtriangulation,&
          iel, ctrafoType, DpointRef=DpointsRef)

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      call elem_generic2 (celement, revalElement, Bder, Dbas)

      ! number of equations in one scalar component
      neqsc = rvectorBlock%RvectorBlock(iblMin)%neq

      ! calculate function values by multiplying the FE-coefficients with the values of
      ! the basis functions and summing up
      Dvalues(iblMin:iblMax,:,ipoint) = 0.0_DP
      if (rvectorBlock%cdataType .eq. ST_DOUBLE) then
        do ider = 1,size(CderType)
          do iblock = iblMin,iblMax
            dval = 0.0_DP
            do ibas = 1,indof
              dval = dval +   p_Ddata((iblock-1)*neqsc + Idofs(ibas)) &
                            * Dbas(ibas,CderType(ider))
            end do
            Dvalues(iblock, ider, ipoint) = dval
          end do
        end do
      else if (rvectorBlock%cdataType .eq. ST_SINGLE) then
        do ider = 1,size(CderType)
          do iblock = iblMin,iblMax
            dval = 0.0_DP
            do ibas = 1,indof
              dval = dval +   p_Fdata((iblock-1)*neqsc + Idofs(ibas)) &
                            * Dbas(ibas,CderType(ider))
            end do
            Dvalues(iblock, ider, ipoint) = dval
          end do
        end do
      end if

    end do ! ipoint

  end subroutine fevl_evaluateBdr2d4

  ! ***************************************************************************

!<subroutine>

  subroutine fevl_getVectorMagnitude (rvector,dumax,rperfconfig)

!<description>
  ! Routine to calculate the vector magnitude of a FE vector field.
!</description>

!<input>
  ! A given finite element vector field. May be e.g. a velocity field.
  type(t_vectorBlock), intent(in) :: rvector

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<output>
  ! The maximum vector field length. E.g. if rvector is a velocity field,
  ! this returns the maximum velocity.
  real(DP) :: dumax
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dvalues,p_Dvalues2
    integer :: i,j
    integer(I32) :: celement
    logical :: bevaluate

    nullify(p_Dvalues)
    nullify(p_Dvalues2)

    ! There are some special cases that we can handle directly.
    ! If all FE spaces are primal spaces (P1,Q1,...), we can just
    ! calculate dumax using lsysbl_getVectorMagnitude.
    ! Otherwise, we project the vector into a Q1 space
    ! and get the maximum value from the values in the vertices.
    bevaluate = .true.

    do i=1,rvector%nblocks
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

    if (bevaluate) then

      ! Evaluate directly
      call lsysbl_getVectorMagnitude (rvector,dumax=dumax)

    else

      ! Calculate the VecMag value by projection into the vertices.
      ! We project each component separately and sum up in p_Dvalues2.
      call spdp_projectToVertices (rvector%RvectorBlock(1), p_Dvalues,&
                                   DER_FUNC, rperfconfig)

      ! If there is more than one block...
      if (rvector%nblocks .gt. 1) then

        ! Calculate val=sqrt(val^2 + val2^2 + val3^2 + ...)

        do i=1,size(p_Dvalues)
          p_Dvalues(i) = p_Dvalues(i)**2
        end do

        do i=2,rvector%nblocks
          call spdp_projectToVertices (rvector%RvectorBlock(i), p_Dvalues2,&
                                       DER_FUNC, rperfconfig)

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

!#########################################################################
!# ***********************************************************************
!# <name> pprocerror </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating errors and norms
!# of finite element functions.
!#
!# The following routines can be found in this module:
!#
!# 1.) pperr_scalarVec
!#     -> Calculate <tex>$ L_1 $</tex>-error, $L_2$-error or $H_1$-error to an
!#        analytic reference function or the <tex>$ L_1 $</tex>-norm, $L_2$-norm
!#        or $H_1$-norm of a FE function:
!#   <tex>
!#   $$ \int_\Omega u-u_h dx , \qquad \int_\Omega \nabla u-\nabla u_h dx $$
!#   </tex>
!#
!# 2.) pperr_scalar
!#     -> Calculate, e.g., <tex>$ L_1 $</tex>-error, $L_2$-error or $H_1$-error to an
!#        analytic reference function or the <tex>$ L_1 $</tex>-norm, $L_2$-norm
!#        or $H_1$-norm of a FE function:
!#   <tex>
!#   $$ \int_\Omega w(u-u_h) dx , \qquad \int_\Omega w(\nabla u-\nabla u_h) dx $$
!#   </tex>
!#
!# 3.) pperr_scalarBoundary2d
!#     -> On a 2D boundary segment, calculate <tex>$ L_1 $</tex>-error, $L_2$-error
!#        or $H_1$-error to an analytic reference function or the
!#        <tex>$ L_1 $</tex>-norm, $L_2$-norm or $H_1$-norm of a FE function.
!#   <tex>
!#   $$ \int_\Gamma u-cu_h dx , \qquad \int_\Gamma \nabla u-c\nabla u_h dx $$
!#   </tex>
!#
!# 4.) pperr_scalarErrorEstimate
!#     -> Calculate error to two different scalar vectors of a  FE function:
!#   <tex>
!#   $$ \int_\Omega u_h-u_ref dx $$
!#   </tex>
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is
!#        some reference solution vector which is supposed to be a
!#        better approximation of the true solution.
!#
!# 5.) pperr_blockErrorEstimate
!#     -> Calculate error to two different block vectors of a  FE function:
!#   <tex>
!#   $$ \int_\Omega u_h-u_ref dx $$
!#   </tex>
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is
!#        some reference solution vector which is supposed to be a
!#        better approximation of the true solution.
!#
!# 6.) pperr_scalarStandardDeviation
!#     -> Calculate the standard deviation of a scalar vector
!#
!# 7.) pperr_blockStandardDeviation
!#     -> Calculate the standard deviation of a block vector
!#
!# 8.) pperr_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# </purpose>
!#########################################################################

module pprocerror

!$use omp_lib
  use basicgeometry
  use boundary
  use boundaryaux
  use collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use feevaluation
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use perfconfig
  use scalarpde
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation
  use extstdassemblyinfo

  implicit none

  private

!<constants>

!<constantblock description = "Identifiers for the type of error to be computed.">

  ! $L_2$-error/norm
  integer, parameter, public :: PPERR_L2ERROR = 1

  ! $H_1$-error/norm
  integer, parameter, public :: PPERR_H1ERROR = 2

  ! <tex>$ L_1 $</tex>-error/norm
  integer, parameter, public :: PPERR_L1ERROR = 3

  ! Mean value error/norm
  integer, parameter, public :: PPERR_MEANERROR = 4

!</constantblock>

!<constantblock description = "Identifiers for the type of weighting function to be computed.">

  ! Weighting function for the error
  integer, parameter, public :: PPERR_ERRORWEIGHT = 1

  ! Weighting function for the reference solution
  integer, parameter, public :: PPERR_REFWEIGHT   = 2

  ! Weighting function for the approximate solution
  integer, parameter, public :: PPERR_FEWEIGHT    = 3

!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of elements to handle simultaneously when building vectors
#ifndef PPERR_NELEMSIM
  integer, parameter, public :: PPERR_NELEMSIM = 256
#endif

!</constantblock>

!</constants>

!<types>

!<typeblock>

  type t_errorScVec

    ! IN: An array of scalar coefficient vectors which represents a
    ! FE function or vector field. If p_rdiscr is null, then all vectors in
    ! the array must have the same spatial discretisation!
    type(t_vectorScalar), dimension(:), pointer :: p_RvecCoeff => null()

    ! IN, OPTIONAL: A spatial discretisation structure that is to be used
    ! for the scalar coefficient vectors. If given, the total number of DOFs
    ! of the spatial discretisation must be equal to the number of entries of
    ! each scalar vector given in p_RvecCoeff. If not given, the spatial
    ! discretisation of p_RvecCoeff is used.
    type(t_spatialDiscretisation), pointer :: p_rdiscr => null()

    ! OUT, OPTIONAL: If given, recieves the calculated L2-errors for each
    ! component of p_RvecCoeff.
    real(DP), dimension(:), pointer :: p_DerrorL2 => null()

    ! OUT, OPTIONAL: If given, recieves the calculated H1-errors for each
    ! component of p_RvecCoeff. If given, but the spatial discretisation does
    ! not offer first derivatives, then all components of p_DerrorH1 are set
    ! to SYS_INFINITY_DP to indicate that the H1-error cannot be calculated.
    real(DP), dimension(:), pointer :: p_DerrorH1 => null()

    ! OUT, OPTIONAL: If given, recieves the calculated L1-errors for each
    ! component of p_RvecCoeff.
    real(DP), dimension(:), pointer :: p_DerrorL1 => null()

    ! OUT, OPTIONAL: If given, recieves the calculated MEAN-errors for each
    ! component of p_RvecCoeff.
    real(DP), dimension(:), pointer :: p_DerrorMean => null()

    ! OUT, OPTIONAL: An array of scalar vectors which recieve the L2-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorL2 => null()

    ! OUT, OPTIONAL: An array of scalar vectors which recieve the H1-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorH1 => null()

    ! OUT, OPTIONAL: An array of scalar vectors which recieve the L1-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorL1 => null()

    ! OUT, OPTIONAL: An array of scalar vectors which recieve the MEAN-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorMean => null()

  end type

  public :: t_errorScVec

!</typeblock>

!</types>

  interface pperr_scalar
    module procedure pperr_scalar_conf
    module procedure pperr_scalarObsolete
  end interface

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: pperr_perfconfig

  !************************************************************************

  public :: pperr_initPerfConfig
  public :: pperr_scalar
  public :: pperr_scalarVec
  public :: pperr_scalarBoundary2d
  public :: pperr_scalarErrorEstimate
  public :: pperr_blockErrorEstimate
  public :: pperr_scalarStandardDeviation
  public :: pperr_blockStandardDeviation

contains

  !****************************************************************************

!<subroutine>

  subroutine pperr_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      pperr_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(pperr_perfconfig)
      pperr_perfconfig%NELEMSIM = PPERR_NELEMSIM
    end if

  end subroutine pperr_initPerfConfig

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarVec(rerror, frefFunction, rcollection, rcubatureInfo, &
      ntempArrays, rperfconfig)

!<description>
  ! This routine calculates the errors of a set of given FE functions given
  ! analytical callback function frefFunction.
  ! In contrast to pperr_scalar, this routine is able to compute multiple
  ! errors in multiple components of a vector field at once.
!</description>

!<input>
  ! OPTIONAL: A callback function that provides the analytical reference
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionScVec.inc'
  optional :: frefFunction

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: Number of temp arrays.
  ! If this is specified (and larger than zero), a number of temporary
  ! arrays is allocated and provided to the callback routine.
  ! This temporary memory is user defined and can be used during the
  ! assembly, e.g., for the temporary evaluation of FEM functions.
  !
  ! The temporary memory is available in the callback function as
  !    rdomainIntSubset%p_DtempArrays !
  integer, intent(in), optional :: ntempArrays

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! A structure which defines what errors are to be calculated.
  type(t_errorScVec), intent(inout), target :: rerror
!</inputoutput>

!</subroutine>

  ! A pointer to the discretisation that is to be used
  type(t_spatialDiscretisation), pointer :: p_rdiscr

  ! A pointer to an element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! An array holding the element list
  integer, dimension(:), pointer :: p_IelementList

  ! A pointer to the triangulation
  type(t_triangulation), pointer :: p_rtria

  ! Which errors do we have to calculate?
  logical :: bcalcL2, bcalcH1, bcalcL1, bcalcMean

  ! The total number of components
  integer :: ncomp,icomp

  ! Indices of first and last derivative
  integer :: ifirstDer, ilastDer

  ! Arrays concerning the cubature formula
  real(DP), dimension(:,:), pointer :: p_DcubPts
  real(DP), dimension(:), pointer :: p_Domega

  ! An array needed for the DOF-mapping
  integer, dimension(:,:), pointer :: p_Idofs

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: revalElementSet

  ! Two arrays for the function values and derivatives
  real(DP), dimension(:,:), allocatable :: DvalFunc
  real(DP), dimension(:,:,:), allocatable :: DvalDer

  ! Two arrays for the element evaluation
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), pointer :: p_Dbas

  ! A pointer to the data array of the currently active coefficient vector
  real(DP), dimension(:), pointer :: p_Dcoeff

  ! Pointers to the arrays that recieve th element-wise errors
  real(DP), dimension(:), pointer :: p_DerrL2, p_DerrH1, p_DerrL1, p_DerrMean

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock,ielementdistr

  ! Some other local variables
  integer :: i,j,k,ndofs,ncubp,iel,ider
  integer :: IELset,IELmax,NEL
  integer(I32) :: ctrafoType,ccubature
  integer(I32) :: cevalTag, celement
  real(DP) :: derrL1, derrL2, derrH1, derrMean, dom, daux, daux2
  real(DP), dimension(:,:), pointer :: p_Ddetj
  integer :: icubatureBlock

  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
  type(t_stdCubatureData) :: rcubatureData
  type(t_stdFEBasisEvalData) :: rfeBasisEvalData

  ! Pointer to temporary memory for callback functions.
  real(DP), dimension(:,:,:), pointer :: p_DtempArrays

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => pperr_perfconfig
    end if

    ! Make sure we have the coefficient vectors
    if(.not. associated(rerror%p_RvecCoeff)) then
      call output_line("Coefficient vectors missing!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
      call sys_halt()
    end if

    ! Get the number of components
    ncomp = ubound(rerror%p_RvecCoeff,1)

    ! Do we have a separate discretisation?
    if(associated(rerror%p_rdiscr)) then
      ! Yes, so check whether the discretisation is compatible to the
      ! coefficient vectors.
      p_rdiscr => rerror%p_rdiscr
      k = dof_igetNDofGlob(p_rdiscr)
      do i = 1, ncomp
        if(k .ne. rerror%p_RvecCoeff(i)%NEQ) then
          call output_line("Discretisation and coefficient vectors incompatible!",&
              OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
          call sys_halt()
        end if
      end do

    else
      ! No, so grab the discretisation of the first coefficient vector.
      p_rdiscr => rerror%p_RvecCoeff(1)%p_rspatialDiscr
      if(.not. associated(p_rdiscr)) then
        call output_line("No discretisation assigned!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
    end if

    ! Get the triangulation and the number of elements
    p_rtria => p_rdiscr%p_rtriangulation
    NEL = p_rtria%NEL

    ! Okay, now that we have the discretisation, determine the dimension
    ! to figure out which is the first and the last derivative we need to
    ! evaluate in the case that we want to compute H1-errors.
    select case(p_rdiscr%ndimension)
    case (NDIM1D)
      ifirstDer = DER_DERIV1D_X
      ilastDer  = DER_DERIV1D_X
    case (NDIM2D)
      ifirstDer = DER_DERIV2D_X
      ilastDer  = DER_DERIV2D_Y
    case (NDIM3D)
      ifirstDer = DER_DERIV3D_X
      ilastDer  = DER_DERIV3D_Z
    case default
      call output_line("Invalid discretisation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
      call sys_halt()
    end select

    ! Now determine which errors we are going to calculate
    bcalcL2 = .false.
    bcalcH1 = .false.
    bcalcL1 = .false.
    bcalcMean = .false.

    if(associated(rerror%p_DerrorL2)) then
      bcalcL2 = .true.
      if(ubound(rerror%p_DerrorL2,1) .lt. ncomp) then
        call output_line("Dimension of p_DerrorL2 array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      rerror%p_DerrorL2 = 0.0_DP
    end if

    if(associated(rerror%p_DerrorH1)) then
      bcalcH1 = .true.
      if(ubound(rerror%p_DerrorH1,1) .lt. ncomp) then
        call output_line("Dimension of p_DerrorH1 array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      rerror%p_DerrorH1 = 0.0_DP
    end if

    if(associated(rerror%p_DerrorL1)) then
      bcalcL1 = .true.
      if(ubound(rerror%p_DerrorL1,1) .lt. ncomp) then
        call output_line("Dimension of p_DerrorL1 array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      rerror%p_DerrorL1 = 0.0_DP
    end if

    if(associated(rerror%p_DerrorMean)) then
      bcalcMean = .true.
      if(ubound(rerror%p_DerrorMean,1) .lt. ncomp) then
        call output_line("Dimension of p_DerrorMean array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      rerror%p_DerrorMean = 0.0_DP
    end if

    if(associated(rerror%p_RvecErrorL2)) then
      bcalcL2 = .true.
      if(ubound(rerror%p_RvecErrorL2,1) .lt. ncomp) then
        call output_line("Dimension of p_RvecErrorL2 array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorL2(i)%NEQ .lt. NEL) then
          call output_line("Length of p_RvecErrorL2("//trim(sys_siL(i,4))//&
              ") is too small!",OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorL2(i))
      end do
    end if

    if(associated(rerror%p_RvecErrorH1)) then
      bcalcH1 = .true.
      if(ubound(rerror%p_RvecErrorH1,1) .lt. ncomp) then
        call output_line("Dimension of p_RvecErrorH1 array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorH1(i)%NEQ .lt. NEL) then
          call output_line("Length of p_RvecErrorH1("//trim(sys_siL(i,4))//&
              ") is too small!",OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorH1(i))
      end do
    end if

    if(associated(rerror%p_RvecErrorL1)) then
      bcalcL1 = .true.
      if(ubound(rerror%p_RvecErrorL1,1) .lt. ncomp) then
        call output_line("Dimension of p_RvecErrorL1 array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorL1(i)%NEQ .lt. NEL) then
          call output_line("Length of p_RvecErrorL1("//trim(sys_siL(i,4))//&
              ") is too small!",OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorL1(i))
      end do
    end if

    if(associated(rerror%p_RvecErrorMean)) then
      bcalcMean = .true.
      if(ubound(rerror%p_RvecErrorMean,1) .lt. ncomp) then
        call output_line("Dimension of p_RvecErrorMean array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorMean(i)%NEQ .lt. NEL) then
          call output_line("Length of p_RvecErrorMean("//trim(sys_siL(i,4))//&
              ") is too small!",OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarVec")
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorMean(i))
      end do
    end if

    ! Do not we have anything to do?
    if(.not. (bcalcL2 .or. bcalcH1 .or. bcalcL1 .or. bcalcMean)) return

    ! Do we have to calculate H1 errors?
    if(bcalcH1) then

      ! Okay, in this case all element distributions of the spatial
      ! discretisation must support first derivatives.
      do i = 1, p_rdiscr%inumFEspaces

        ! Get a pointer to the element distribution
        p_relementDistribution => p_rdiscr%RelementDistr(i)

        ! If the maximum supported derivative is 1, then the element does not
        ! support first derivatives!
        if(elem_getMaxDerivative(p_relementDistribution%celement) .le. 1) then
          bcalcH1 = .false.
          exit
        end if
      end do

      ! Now if bcalcH1 is .false. now, then at least one element distribution
      ! does not support first derivatives...
      if(.not. bcalcH1) then
        ! If p_DerrorH1 is given, set its entries to SYS_INFTY to indicate that
        ! the H1-errors are not available
        if(associated(rerror%p_DerrorH1)) &
          rerror%p_DerrorH1 = SYS_INFINITY_DP

      end if

    end if

    ! Set up the Bder array
    Bder = .false.
    Bder(DER_FUNC) = bcalcL2 .or. bcalcL1 .or. bcalcMean
    if(bcalcH1) Bder(ifirstDer:ilastDer) = .true.

    ! Do we have an assembly structure?
    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscr,rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Loop over the cubature blocks. Each defines a separate cubature formula.
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get typical information: Number of elements, element list,...
      call spdiscr_getStdDiscrInfo (icubatureBlock,p_rcubatureInfo,p_rdiscr,&
          ielementDistr,celement,ccubature,NEL,p_IelementList)

      ! Cancel if this element list is empty.
      if (NEL .le. 0) cycle

      ! For saving some memory in smaller discretisations, we calculate
      ! the number of elements per block. For smaller triangulations,
      ! this is NEL. If there are too many elements, it is at most
      ! NELEMSIM. This is only used for allocating some arrays.
      nelementsPerBlock = min(p_rperfconfig%NELEMSIM, NEL)

      ! Initialise cubature and evaluation of the FE basis
      call easminfo_initStdCubature(ccubature,rcubatureData)

      ! Get the pointers for faster access
      p_DcubPts => rcubatureData%p_DcubPtsRef
      p_Domega => rcubatureData%p_Domega
      ncubp = rcubatureData%ncubp

      ! Get the trafo
      ctrafoType = elem_igetTrafoType(celement)

      ! Get the evaluation tag
      cevalTag = elem_getEvaluationTag(celement)

      ! If a reference function is given, it will surely need real points.
      if(present(frefFunction)) &
        cevalTag = ior(cevalTag, EL_EVLTAG_REALPOINTS)

      ! And we definately need jacobian determinants for integration.
      cevalTag = ior(cevalTag, EL_EVLTAG_DETJ)

      ! OpenMP-Extension: Open threads here.
      ! Each thread will allocate its own local memory...
      !
      !$omp parallel default(shared)&
      !$omp private(DvalDer,DvalFunc,IELmax,daux,daux2,dom,&
      !$omp         i,ider,iel,j,k,ndofs,p_Dbas,p_Ddetj,p_Idofs,revalElementSet,&
      !$omp         rfeBasisEvalData,rintSubset,p_DtempArrays,&
      !$omp         icomp,p_Dcoeff,p_DerrL2,p_DerrH1,p_DerrL1,p_DerrMean)&
      !$omp         firstprivate(cevalTag)&
      !$omp         reduction(+:derrL1,derrL2,derrH1,derrMean)&
      !$omp if (NEL > p_rperfconfig%NELEMMIN_OMP)

      ! Allocate memory for user defined data.
      nullify(p_DtempArrays)
      if (present(ntempArrays)) then
        allocate(p_DtempArrays(ncubp,nelementsPerBlock,ntempArrays))
      end if

      ! Initialise the evaluation structure for the FE basis
      call easminfo_initStdFEBasisEval(celement,&
          elem_getMaxDerivative(celement),ncubp,nelementsPerBlock,rfeBasisEvalData)

      ! Quick reference to the values of the basis functions
      p_Dbas => rfeBasisEvalData%p_Dbas
      p_Idofs => rfeBasisEvalData%p_Idofs
      ndofs = rfeBasisEvalData%ndofLocal

      ! Allocate two arrays for the evaluation
      if(bcalcL2 .or. bcalcL1 .or. bcalcMean) &
        allocate(DvalFunc(ncubp,nelementsPerBlock))
        
      if(bcalcH1) &
        allocate(DvalDer(ncubp,nelementsPerBlock,ifirstDer:ilastDer))

      ! Initialise the element evaluation set
      call elprep_init(revalElementSet)

      ! Loop over the elements - blockwise.
      !$omp do schedule(static,1)
      do IELset = 1, NEL, nelementsPerBlock

        ! We always handle nelementsPerBlock elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most
        ! nelementsPerBlock elements simultaneously.

        IELmax = min(NEL,IELset-1+nelementsPerBlock)

        ! First, let us perform the DOF-mapping
        call dof_locGlobMapping_mult(p_rdiscr, p_IelementList(IELset:IELmax), p_Idofs)

        ! Prepare the element for evaluation
        call elprep_prepareSetForEvaluation (revalElementSet, cevalTag, p_rtria, &
            p_IelementList(IELset:IELmax), ctrafoType, p_DcubPts, &
            rperfconfig=p_rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj(:,1:IELmax-IELset+1)

        ! Remove the ref-points eval tag for the next loop iteration
        cevalTag = iand(cevalTag,not(EL_EVLTAG_REFPOINTS))

        ! Prepare the domain integration structure
        call domint_initIntegrationByEvalSet(revalElementSet, rintSubset)

        ! Associate temp memory for the callback.
        if (associated(p_DtempArrays)) then
          call domint_setTempMemory (rintSubset,p_DtempArrays)
        end if

        !rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs(:,1:IELmax-IELset+1)
        rintSubset%celement = celement

        ! Evaluate the element
        call elem_generic_sim2(celement, revalElementSet, Bder, p_Dbas)

        ! Now loop over all vector components
        do icomp = 1, ncomp

          ! Get the coefficient vector`s data array
          call lsyssc_getbase_double(rerror%p_RvecCoeff(icomp), p_Dcoeff)

          ! Get the element-wise error arrays, if given
          if(associated(rerror%p_RvecErrorL2)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorL2(icomp), p_DerrL2)
          else
            nullify(p_DerrL2)
          end if
          if(associated(rerror%p_RvecErrorH1)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorH1(icomp), p_DerrH1)
          else
            nullify(p_DerrH1)
          end if
          if(associated(rerror%p_RvecErrorL1)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorL1(icomp), p_DerrL1)
          else
            nullify(p_DerrL1)
          end if
          if(associated(rerror%p_RvecErrorMean)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorMean(icomp), p_DerrMean)
          else
            nullify(p_DerrMean)
          end if

          ! Reset errors for this component
          derrL1 = 0.0_DP
          derrL2 = 0.0_DP
          derrH1 = 0.0_DP
          derrMean = 0.0_DP

          ! Evaluate function values?
          if(allocated(DvalFunc)) then

            ! Do we have a reference function? If yes, then evaluate it,
            ! otherwise simply format DvalFunc to zero.
            if(present(frefFunction)) then
              call frefFunction(icomp,DER_FUNC,p_rdiscr,IELmax-IELset+1,ncubp,&
                  revalElementSet%p_DpointsReal,rintSubset,&
                  DvalFunc(:,1:IELmax-IELset+1),rcollection)
            else
              DvalFunc = 0.0_DP
            end if

            ! Now subtract the function values of the FE function.
            do j = 1,IELmax-IELset+1
              do i = 1, ncubp
                daux = 0.0_DP
                do k = 1, ndofs
                  daux = daux + p_Dbas(k,DER_FUNC,i,j)*p_Dcoeff(p_Idofs(k,j))
                end do ! k
                DvalFunc(i,j) = DvalFunc(i,j) - daux
              end do ! i
            end do ! j

          end if ! function values evaluation

          ! Evaluate derivatives?
          if(allocated(DvalDer)) then

            ! Do we have a reference function? If yes, then evaluate its
            ! derivatives.
            if(present(frefFunction)) then
              do ider = ifirstDer, ilastDer
                call frefFunction(icomp,ider,p_rdiscr,IELmax-IELset+1,ncubp, &
                    revalElementSet%p_DpointsReal,rintSubset,&
                    DvalDer(:,1:IELmax-IELset+1,ider),rcollection)
              end do
            else
              DvalDer = 0.0_DP
            end if

            ! Now subtract the derivatives of the FE function.
            do j = 1,IELmax-IELset+1
              do i = 1, ncubp
                do ider = ifirstDer, ilastDer
                  daux = 0.0_DP
                  do k = 1, ndofs
                    daux = daux + p_Dbas(k,ider,i,j)*p_Dcoeff(p_Idofs(k,j))
                  end do ! k
                  DvalDer(i,j,ider) = DvalDer(i,j,ider) - daux
                end do ! ider
              end do ! i
            end do ! j

          end if ! derivatives evaluation

          ! Do we calculate L2-errors?
          if(bcalcL2 .and. associated(p_DerrL2)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)**2
              end do ! i
              p_DerrL2(iel) = sqrt(daux)
              derrL2 = derrL2 + daux
            end do ! j
          else if(bcalcL2) then
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)**2
              end do ! i
              derrL2 = derrL2 + daux
            end do ! j
          end if

          ! Do we calculate H1-errors?
          if(bcalcH1 .and. associated(p_DerrH1)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux2 = 0.0_DP
                do ider = ifirstDer, ilastDer
                  daux2 = daux2 + DvalDer(i,j,ider)**2
                end do ! ider
                daux = daux + dom*daux2
              end do ! i
              p_DerrH1(iel) = sqrt(daux)
              derrH1 = derrH1 + daux
            end do ! j
          else if(bcalcH1) then
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux2 = 0.0_DP
                do ider = ifirstDer, ilastDer
                  daux2 = daux2 + DvalDer(i,j,ider)**2
                end do ! ider
                daux = daux + dom*daux2
              end do ! i
              derrH1 = derrH1 + daux
            end do ! j
          end if

          ! Do we calculate L1-errors?
          if(bcalcL1 .and. associated(p_DerrL1)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*abs(DvalFunc(i,j))
              end do ! i
              p_DerrL1(iel) = daux
              derrL1 = derrL1 + daux
            end do ! j
          else if(bcalcL1) then
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*abs(DvalFunc(i,j))
              end do ! i
              derrL1 = derrL1 + daux
            end do ! j
          end if

          ! Do we calculate Mean-errors?
          if(bcalcMean .and. associated(p_DerrMean)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)
              end do ! i
              p_DerrMean(iel) = daux
              derrMean = derrMean + daux
            end do ! j
          else if(bcalcMean) then
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)
              end do ! i
              derrMean = derrMean + daux
            end do ! j
          end if

          !$omp critical
          if(bcalcL2 .and. associated(rerror%p_DerrorL2)) &
              rerror%p_DerrorL2(icomp) = rerror%p_DerrorL2(icomp) + derrL2
          if(bcalcH1 .and. associated(rerror%p_DerrorH1)) &
              rerror%p_DerrorH1(icomp) = rerror%p_DerrorH1(icomp) + derrH1
          if(bcalcL1 .and. associated(rerror%p_DerrorL1)) &
              rerror%p_DerrorL1(icomp) = rerror%p_DerrorL1(icomp) + derrL1
          if(bcalcMean .and. associated(rerror%p_DerrorMean)) &
              rerror%p_DerrorMean(icomp) = rerror%p_DerrorMean(icomp) + derrMean
          !$omp end critical

        end do ! icomp

        ! Release the domain integration structure
        call domint_doneIntegration (rintSubset)

      end do ! IELset
      !$omp end do

      ! Release the element evaluation set
      call elprep_releaseElementSet(revalElementSet)

      ! Release FE evaluation
      call easminfo_doneStdFEBasisEval(rfeBasisEvalData)

      ! Release the temp memory
      if (associated(p_DtempArrays)) then
        deallocate(p_DtempArrays)
      end if

      ! Deallocate all arrays
      if(allocated(DvalDer)) deallocate(DvalDer)
      if(allocated(DvalFunc)) deallocate(DvalFunc)
      !$omp end parallel

      ! Release cubature information
      call easminfo_doneStdCubature(rcubatureData)

    end do ! icurrentElementDistr

    ! Do not forget to take the square roots of the L2- and H1-errors
    if(associated(rerror%p_DerrorL2) .and. bcalcL2) then
      do icomp = 1, ncomp
        rerror%p_DerrorL2(icomp) = sqrt(rerror%p_DerrorL2(icomp))
      end do
    end if
    if(associated(rerror%p_DerrorH1) .and. bcalcH1) then
      do icomp = 1, ncomp
        rerror%p_DerrorH1(icomp) = sqrt(rerror%p_DerrorH1(icomp))
      end do
    end if

    ! That is it

  end subroutine pperr_scalarVec

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarObsolete (rvector, cerrortype, derror,&
                                   ffunctionReference, rcollection,&
                                   rdiscretisation, relementError,&
                                   rperfconfig)

!<description>
  ! This routine calculates the error or the norm, respectively, of a
  ! given finite element function in rvector to a given analytical
  ! callback function ffunctionReference.
  !
  ! If ffunctionReference is specified, the routine calculates
  !   <tex> $$ ||y-z||_{L_1}, ||y-z||_{L_2}  \textrm{ or }  ||y-z||_{H_1} $$ </tex>
  ! with $y$=rvectorScalar and $z$=ffunctionReference.
  !
  ! If ffunctionReference is not specified, the routine calculates
  !   <tex> $$ ||y||_{L_1}, ||y||_{L_2}  \textrm{ or }  ||y||_{H_1} $$ </tex>
  !
  ! If ffunctionWeight is specified, the routine calculates the
  ! desired norm over the selected subdomain and/or scales the error
  ! by the local weighting coefficients.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
  !
  ! If the H1-error is desired and the element does not provide first
  ! derivatives, then this routine sets derror to -1 to indicate that the
  ! calculation of the H1-error is not available.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(in), target :: rvector

  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(in) :: cerrortype

  ! OPTIONAL: A callback function that provides the analytical reference
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionSc.inc'
  optional :: ffunctionReference

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! If not specified, the discretisation structure in the vector is used.
  ! If specified, the discretisation structure must be 'compatible' to the
  ! vector (concerning NEQ,...). pperr_scalar uses the cubature formula
  ! specifier of the linear form in rdiscretisation to compute the integrals
  ! for the error.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! OPTIONAL: A scalar vector which holds the calculated error per element
  type(t_vectorScalar), intent(inout), optional :: relementError
!</inputoutput>

!<output>
  ! The calculated error.
  real(DP), intent(out) :: derror
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_scalarCubatureInfo), target :: rcubatureInfo

#if WARN_DEPREC
    call output_line ("Using deprecated feature. Please update your code.", &
        OU_CLASS_WARNING,OU_MODE_STD,"pperr_scalarObsolete")
#endif

    ! Create an assembly structure and take the associated cubature rule.
    if (present(rdiscretisation)) then
      call spdiscr_createDefCubStructure(rdiscretisation,&
          rcubatureInfo,CUB_GEN_DEPR_EVAL)
    else
      call spdiscr_createDefCubStructure(rvector%p_rspatialDiscr,&
          rcubatureInfo,CUB_GEN_DEPR_EVAL)   
    end if

    if (present(relementError)) then

      ! Get the data array for the error
      call lsyssc_getbase_double (relementError,p_Ddata)

      ! Call the new routine with diffenret ordering of parameters
      call pperr_scalar_conf(cerrortype, derror, rvector,&
                      ffunctionReference, rcollection,&
                      DelementError=p_Ddata,rcubatureInfo=rcubatureInfo,&
                      rperfconfig=rperfconfig)
    else

      ! Call the new routine with diffenret ordering of parameters
      call pperr_scalar_conf(cerrortype, derror, rvector,&
                      ffunctionReference, rcollection,&
                      rcubatureInfo=rcubatureInfo,&
                      rperfconfig=rperfconfig)

    end if

    ! Release the temporary assembly structure.
    call spdiscr_releaseCubStructure(rcubatureInfo)

  end subroutine pperr_scalarObsolete

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalar_conf (cerrortype, derror, rvectorScalar, ffunctionReference,&
      rcollection, ffunctionWeight, DelementError, rcubatureInfo, ntempArrays, rperfconfig)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(in) :: cerrortype

  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(in), target :: rvectorScalar

  ! OPTIONAL: A callback function that provides the analytical reference
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionSc.inc'
  optional :: ffunctionReference

  ! OPTIONAL: A callback function that provides the weighting function
  ! by which the computed error is multipled.
  ! If not specified, the reference function is assumed to be =1!
  optional :: ffunctionWeight

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: Number of temp arrays.
  ! If this is specified (and larger than zero), a number of temporary
  ! arrays is allocated and provided to the callback routine.
  ! This temporary memory is user defined and can be used during the
  ! assembly, e.g., for the temporary evaluation of FEM functions.
  !
  ! The temporary memory is available in the callback function as
  !    rdomainIntSubset%p_DtempArrays !
  integer, intent(in), optional :: ntempArrays

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information for callback routines.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: A scalar array which holds the calculated error per element.
  ! The length of this array must be as large as there are cells in
  ! the underlying triangulation; each entry corresponds to one cell
  ! in the discretisation.
  real(DP), dimension(:), intent(inout), optional :: DelementError
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(out) :: derror
!</output>

!</subroutine>

  ! A pointer to the discretisation that is to be used
  type(t_spatialDiscretisation), pointer :: p_rdiscr

  ! An array holding the element list
  integer, dimension(:), pointer :: p_IelementList

  ! A pointer to the triangulation
  type(t_triangulation), pointer :: p_rtria

  ! Indices of first and last derivative
  integer :: ifirstDer, ilastDer

  ! Arrays concerning the cubature formula
  real(DP), dimension(:,:), pointer :: p_DcubPts
  real(DP), dimension(:), pointer :: p_Domega

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: revalElementSet

  ! Two arrays for the function values and derivatives
  real(DP), dimension(:,:), allocatable :: DvalFunc
  real(DP), dimension(:,:), allocatable :: DvalWeight
  real(DP), dimension(:,:,:), allocatable :: DvalDer

  ! Array for the element evaluation
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), pointer :: p_Dbas
  integer, dimension(:,:), pointer :: p_Idofs

  ! A pointer to the data array of the currently active coefficient vector
  real(DP), dimension(:), pointer :: p_Dcoeff

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! Current assembly block, cubature formula, element type,...
  integer :: ielementDistr,icubatureBlock
  integer(I32) :: cevalTag, celement, ccubature, ctrafoType
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
  type(t_stdCubatureData) :: rcubatureData
  type(t_stdFEBasisEvalData) :: rfeBasisEvalData

  ! Pointer to temporary memory for callback functions.
  real(DP), dimension(:,:,:), pointer :: p_DtempArrays

  ! Some other local variables
  integer :: i,j,k,ndofs,ncubp,iel,ider
  integer :: IELset,IELmax,NEL
  real(DP) :: dom,daux,daux2
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => pperr_perfconfig
    end if

    if (rvectorScalar%bisSorted) then
      call output_line("Vector must be unsorted!",&
          OU_CLASS_ERROR, OU_MODE_STD, "pperr_scalar_conf")
      call sys_halt()
    end if

    ! Get the discretisation
    p_rdiscr => rvectorScalar%p_rspatialDiscr

    if(.not. associated(p_rdiscr)) then
      call output_line("No discretisation assigned!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalar_conf")
      call sys_halt()
    end if

    if ((p_rdiscr%ccomplexity .ne. SPDISC_UNIFORM) .and.&
        (p_rdiscr%ccomplexity .ne. SPDISC_CONFORMAL)) then
      call output_line("Unsupported discretisation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalar_conf")
      call sys_halt()
    end if

    ! Get the triangulation
    p_rtria => p_rdiscr%p_rtriangulation

    ! Okay, now that we have the discretisation, determine the dimension
    ! to figure out which is the first and the last derivative we need to
    ! evaluate in the case that we want to compute H1-errors.
    select case(p_rdiscr%ndimension)
    case (NDIM1D)
      ifirstDer = DER_DERIV1D_X
      ilastDer  = DER_DERIV1D_X
    case (NDIM2D)
      ifirstDer = DER_DERIV2D_X
      ilastDer  = DER_DERIV2D_Y
    case (NDIM3D)
      ifirstDer = DER_DERIV3D_X
      ilastDer  = DER_DERIV3D_Z
    case default
      call output_line("Invalid discretisation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalar_conf")
      call sys_halt()
    end select

    ! Reset the error
    derror = 0.0_DP

    ! Assure that the element error array is large enough.
    if(present(DelementError)) then
      if(ubound(DelementError,1) .lt. p_rtria%NEL) then
        call output_line("Dimension of Derror array is too small!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalar_conf")
        call sys_halt()
      end if

      ! Reset the error
      DelementError(:) = 0.0_DP
    end if

    ! Do we have to calculate H1 errors?
    if(cerrortype .eq. PPERR_H1ERROR) then

      ! Okay, in this case all element distributions of the spatial
      ! discretisation must support first derivatives.
      do i = 1, p_rdiscr%inumFEspaces

        ! If the maximum supported derivative is 1, then the element does not
        ! support first derivatives!
        if(elem_getMaxDerivative(p_rdiscr%RelementDistr(i)%celement) .le. 1) then
          call output_line(&
              "Finite element does not support the calculation of the H1 error!",&
              OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalar_conf")
          if (present(DelementError)) then
            DelementError(:) = SYS_INFINITY_DP
          end if
          return
        end if
      end do

    end if

    ! Set up the Bder array
    Bder = .false.
    Bder(DER_FUNC) = (cerrortype .eq. PPERR_L2ERROR) .or. &
                     (cerrortype .eq. PPERR_L1ERROR) .or. &
                     (cerrortype .eq. PPERR_MEANERROR)
    if (cerrortype .eq. PPERR_H1ERROR) then
      Bder(ifirstDer:ilastDer) = .true.
    end if

    ! Do we have an assembly structure?
    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscr,rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Get the coefficient vector`s data array
    call lsyssc_getbase_double(rvectorScalar, p_Dcoeff)

    ! Loop over the element blocks. Each defines a separate cubature formula.
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get typical information: Number of elements, element list,...
      call spdiscr_getStdDiscrInfo (icubatureBlock,p_rcubatureInfo,p_rdiscr,&
          ielementDistr,celement,ccubature,NEL,p_IelementList)

      ! Cancel if this element list is empty.
      if (NEL .le. 0) cycle

      ! For saving some memory in smaller discretisations, we calculate
      ! the number of elements per block. For smaller triangulations,
      ! this is NEL. If there are too many elements, it is at most
      ! NELEMSIM. This is only used for allocating some arrays.
      nelementsPerBlock = min(p_rperfconfig%NELEMSIM, NEL)

      ! Initialise cubature and evaluation of the FE basis
      call easminfo_initStdCubature(ccubature,rcubatureData)

      ! Get the pointers for faster access
      p_DcubPts => rcubatureData%p_DcubPtsRef
      p_Domega => rcubatureData%p_Domega
      ncubp = rcubatureData%ncubp

      ! Get the trafo
      ctrafoType = elem_igetTrafoType(celement)

      ! Get the evaluation tag
      cevalTag = elem_getEvaluationTag(celement)

      ! If a reference function or a weighting function is given, 
      ! it will surely need real points.
      if(present(ffunctionReference)) &
        cevalTag = ior(cevalTag, EL_EVLTAG_REALPOINTS)

      if(present(ffunctionWeight)) &
        cevalTag = ior(cevalTag, EL_EVLTAG_REALPOINTS)

      ! And we definately need jacobian determinants for integration.
      cevalTag = ior(cevalTag, EL_EVLTAG_DETJ)

      ! OpenMP-Extension: Open threads here.
      ! Each thread will allocate its own local memory...
      !
      !$omp parallel default(shared)&
      !$omp private(DvalDer,DvalFunc,DvalWeight,IELmax,daux,daux2,dom,&
      !$omp         i,ider,iel,j,k,ndofs,p_Dbas,p_Ddetj,p_Idofs,revalElementSet,&
      !$omp         rfeBasisEvalData,rintSubset,p_DtempArrays)&
      !$omp         firstprivate(cevalTag)&
      !$omp         reduction(+:derror)&
      !$omp if (NEL > p_rperfconfig%NELEMMIN_OMP)

      ! Allocate memory for user defined data.
      nullify(p_DtempArrays)
      if (present(ntempArrays)) then
        allocate(p_DtempArrays(ncubp,nelementsPerBlock,ntempArrays))
      end if

      ! Initialise the evaluation structure for the FE basis
      call easminfo_initStdFEBasisEval(celement,&
          elem_getMaxDerivative(celement),ncubp,nelementsPerBlock,rfeBasisEvalData)

      ! Quick reference to the values of the basis functions
      p_Dbas => rfeBasisEvalData%p_Dbas
      p_Idofs => rfeBasisEvalData%p_Idofs
      ndofs = rfeBasisEvalData%ndofLocal

      ! Allocate two arrays for the evaluation
      if ((cerrortype .eq. PPERR_L2ERROR) .or. &
          (cerrortype .eq. PPERR_L1ERROR) .or. &
          (cerrortype .eq. PPERR_MEANERROR)) then
        allocate(DvalFunc(ncubp,nelementsPerBlock))
      end if

      if (cerrortype .eq. PPERR_H1ERROR) then
        allocate(DvalDer(ncubp,nelementsPerBlock,ifirstDer:ilastDer))
      end if

      ! If a weighting function is specified, we need an additional
      ! array for its computed values
      if(present(ffunctionWeight)) then
        allocate(DvalWeight(ncubp,nelementsPerBlock))
      end if

      ! Initialise the element evaluation set
      call elprep_init(revalElementSet)

      ! Loop over the elements - blockwise.
      !$omp do schedule(static,1)
      do IELset = 1, NEL, nelementsPerBlock

        ! We always handle nelementsPerBlock elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most
        ! nelementsPerBlock elements simultaneously.

        IELmax = min(NEL,IELset-1+nelementsPerBlock)

        ! First, let us perform the DOF-mapping
        call dof_locGlobMapping_mult(p_rdiscr, p_IelementList(IELset:IELmax), p_Idofs)

        ! Prepare the element for evaluation
        call elprep_prepareSetForEvaluation (revalElementSet, cevalTag, p_rtria, &
            p_IelementList(IELset:IELmax), ctrafoType, p_DcubPts, &
            rperfconfig=p_rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj(:,1:IELmax-IELset+1)

        ! Remove the ref-points eval tag for the next loop iteration
        cevalTag = iand(cevalTag,not(EL_EVLTAG_REFPOINTS))

        ! Prepare the domain integration structure
        call domint_initIntegrationByEvalSet(revalElementSet, rintSubset)

        ! Associate temp memory for the callback.
        if (associated(p_DtempArrays)) then
          call domint_setTempMemory (rintSubset,p_DtempArrays)
        end if

        !rintSubset%ielementDistribution = ielementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs(:,1:IELmax-IELset+1)
        rintSubset%celement = celement

        ! Evaluate the element
        call elem_generic_sim2(celement, revalElementSet, Bder, p_Dbas)

        ! Evaluate function values?
        if(allocated(DvalFunc)) then

          ! Do we have a reference function? If yes, then evaluate it,
          ! otherwise simply format DvalFunc to zero.
          if(present(ffunctionReference)) then
            call ffunctionReference(DER_FUNC,p_rdiscr,IELmax-IELset+1,ncubp,&
                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
                p_Idofs(:,1:IELmax-IELset+1),rintSubset,&
                DvalFunc(:,1:IELmax-IELset+1),rcollection)
          else
            DvalFunc = 0.0_DP
          end if

          ! Evaluate the FEM function.
          ! Subtract the function values of the FE function.
          do j = 1,IELmax-IELset+1
            do i = 1, ncubp
              daux = 0.0_DP
              do k = 1, ndofs
                daux = daux + p_Dbas(k,DER_FUNC,i,j)*p_Dcoeff(p_Idofs(k,j))
              end do ! k
              DvalFunc(i,j) = DvalFunc(i,j) - daux
            end do ! i
          end do ! j

          ! If a weighting function is specified, evaluate it and multiply
          if(present(ffunctionWeight)) then
            call ffunctionWeight(p_rdiscr,IELmax-IELset+1,ncubp,&
                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
                p_Idofs(:,1:IELmax-IELset+1),rintSubset,&
                DvalWeight(:,1:IELmax-IELset+1),rcollection)

            do j = 1,IELmax-IELset+1
              do i = 1, ncubp
                DvalFunc(i,j) = DvalFunc(i,j) * DvalWeight(i,j)
              end do ! i
            end do ! j

          end if

        end if ! function values evaluation

        ! Evaluate derivatives?
        if(allocated(DvalDer)) then

          ! Do we have a reference function? If yes, then evaluate its
          ! derivatives.
          if(present(ffunctionReference)) then
            do ider = ifirstDer, ilastDer
              call ffunctionReference(ider,p_rdiscr,IELmax-IELset+1,ncubp, &
                  revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
                  p_Idofs(:,1:IELmax-IELset+1),rintSubset,&
                  DvalDer(:,1:IELmax-IELset+1,ider),rcollection)
            end do
          else
            DvalDer = 0.0_DP
          end if

          ! Evaluate the FEM function.
          ! Subtract the function values of the FE function.
          do j = 1,IELmax-IELset+1
            do i = 1, ncubp
              do ider = ifirstDer, ilastDer
                daux = 0.0_DP
                do k = 1, ndofs
                  daux = daux + p_Dbas(k,ider,i,j)*p_Dcoeff(p_Idofs(k,j))
                end do ! k
                DvalDer(i,j,ider) = DvalDer(i,j,ider) - daux
              end do ! ider
            end do ! i
          end do ! j

          ! If a weighting function is specified, evaluate it and multiply
          if(present(ffunctionWeight)) then
            call ffunctionWeight(p_rdiscr,IELmax-IELset+1,ncubp,&
                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
                p_Idofs(:,1:IELmax-IELset+1),rintSubset,&
                DvalWeight(:,1:IELmax-IELset+1),rcollection)

            do ider = ifirstDer, ilastDer
              do j = 1,IELmax-IELset+1
                do i = 1, ncubp
                  DvalDer(i,j,ider) = DvalDer(i,j,ider) * DvalWeight(i,j)
                end do ! i
              end do ! j
            end do ! ider
          end if

        end if ! derivatives evaluation

        ! Which error to calculate?
        select case (cerrorType)
        case (PPERR_MEANERROR)

          ! Calculate the MEAN error
          if(present(DelementError)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)
              end do ! i
              DelementError(iel) = daux
              derror = derror + daux
            end do ! j
          else
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)
              end do ! i
              derror = derror + daux
            end do ! j
          end if

        case (PPERR_L1ERROR)

          ! Calculate the L1 error
          if(present(DelementError)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*abs(DvalFunc(i,j))
              end do ! i
              DelementError(iel) = daux
              derror = derror + daux
            end do ! j
          else 
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*abs(DvalFunc(i,j))
              end do ! i
              derror = derror + daux
            end do ! j
          end if

        case (PPERR_L2ERROR)

          ! Calculate the L2 error
          if(present(DelementError)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)**2
              end do ! i
              DelementError(iel) = sqrt(daux)
              derror = derror + daux
            end do ! j
          else
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)**2
              end do ! i
              derror = derror + daux
            end do ! j
          end if

        case (PPERR_H1ERROR)
          ! Do we calculate H1-errors?
          if(present(DelementError)) then
            do j = 1,IELmax-IELset+1
              iel = p_IelementList(IELset+j-1)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))

                ! Sum up the squared derivatives.
                daux2 = 0.0_DP
                do ider = ifirstDer, ilastDer
                  daux2 = daux2 + DvalDer(i,j,ider)**2
                end do ! ider
                daux = daux + dom*daux2
              end do ! i
              DelementError(iel) = sqrt(daux)
              derror = derror + daux
            end do ! j
          else 
            do j = 1,IELmax-IELset+1
              daux = 0.0_DP
              do i = 1, ncubp
                dom = p_Domega(i) * abs(p_Ddetj(i,j))

                ! Sum up the squared derivatives.
                daux2 = 0.0_DP
                do ider = ifirstDer, ilastDer
                  daux2 = daux2 + DvalDer(i,j,ider)**2
                end do ! ider
                daux = daux + dom*daux2
              end do ! i
              derror = derror + daux
            end do ! j
          end if

        end select

        ! Release the domain integration structure
        call domint_doneIntegration (rintSubset)

      end do ! IELset
      !$omp end do

      ! Release the element evaluation set
      call elprep_releaseElementSet(revalElementSet)

      ! Release FE evaluation
      call easminfo_doneStdFEBasisEval(rfeBasisEvalData)

      ! Release the temp memory
      if (associated(p_DtempArrays)) then
        deallocate(p_DtempArrays)
      end if

      ! Deallocate all arrays
      if(allocated(DvalDer)) deallocate(DvalDer)
      if(allocated(DvalFunc)) deallocate(DvalFunc)
      if(allocated(DvalWeight)) deallocate(DvalWeight)
      !$omp end parallel

      ! Release cubature information
      call easminfo_doneStdCubature(rcubatureData)

    end do ! icubatureBlock

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

    ! Do not forget to take the square roots of the L2- and H1-errors.
    ! derror is ||error||^2, so take the square root at last.
    if ((cerrortype .eq. PPERR_L2ERROR) .or. (cerrortype .eq. PPERR_H1ERROR)) then
      derror = sqrt(derror)
    end if

    ! That is it

  end subroutine pperr_scalar_conf

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarBoundary2D (cerrortype, ccubType, derror,&
                                     rboundaryRegion, rvectorScalar,&
                                     ffunctionReference, rcollection,&
                                     rdiscretisation, ffunctionWeight)

!<description>
  ! This routine calculates the error or the norm, respectively, of a given
  ! finite element function in rvector to a given analytical
  ! callback function ffunctionReference.
  !
  ! If ffunctionReference is specified, the routine calculates
  !   <tex> $$ ||y-z||_{L_2}  \textrm{ or }  ||y-z||_{L_1}  \textrm{ or }  ||y-z||_{H_1} $$ </tex>
  ! with $y$=rvectorScalar and $z$=ffunctionReference.
  !
  ! If ffunctionReference is not specified, the routine calculates
  !   <tex> $$ ||y||_{L_2}  \textrm{ or }  ||y||_{L_1}  \textrm{ or }  ||y||_{H_1}. $$ </tex>
  !
  ! If the vector rvectorScalar is not specified, it is assumed to be =0.
  !
  ! If ffunctionWeight is specified, the routine calculates the
  ! desired norm over the selected boundary and/or scales the error
  ! by the local weighting coefficients.
  !
  ! rboundaryRegion is a t_boundaryRegion object that allows to
  ! specify the boundary region where the error should be computed.
  ! If not specified, the error is computed over the whole boundary.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(in) :: cerrortype

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType

  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion

  ! OPTIONAL: The FE solution vector. Represents a scalar FE function.
  ! If omitted, the function is assumed to be constantly =1.
  type(t_vectorScalar), intent(in), optional, target :: rvectorScalar

  ! OPTIONAL: A callback function that provides the analytical reference
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include "intf_refFunctionScBdr2D.inc"
  optional :: ffunctionReference

  ! OPTIONAL: A callback function that provides the weighting function
  ! by which the computed error is multipled.
  ! If not specified, the reference function is assumed to be =1!
  optional :: ffunctionWeight

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! Must be specified if rvectorScalar is not specified as this
  ! describes the domain/triangulation/...
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation
!</input>

!<output>
  ! The calculated error.
  real(DP), intent(out) :: derror
!</output>

!</subroutine>

    ! local variables
    type(t_boundaryRegion) :: rboundaryReg
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(DP) :: dlocalError
    integer :: ibdc,cdataType

    ! Get the correct discretisation structure and check if we can use it.
    if (present(rdiscretisation)) then
      p_rdiscretisation => rdiscretisation
      if (present(rvectorScalar)) then
        call lsyssc_checkDiscretisation (rvectorScalar, p_rdiscretisation)
      end if
    elseif (present(rvectorScalar)) then
      p_rdiscretisation => rvectorScalar%p_rspatialdiscr
    else
      nullify(p_rdiscretisation)
    end if

    if (.not. associated(p_rdiscretisation)) then
      call output_line("No discretisation structure!",&
                       OU_CLASS_ERROR, OU_MODE_STD, "pperr_scalarBoundary2D")
      call sys_halt()
    end if

    if (p_rdiscretisation%ndimension .ne. NDIM2D) then
      call output_line("Only 2D discretisations allowed.",&
                       OU_CLASS_ERROR, OU_MODE_STD, "pperr_scalarBoundary2D")
      call sys_halt()
    end if

    ! The vector must be unsorted, otherwise we can not set up the vector.
    if (present(rvectorScalar)) then
      if (rvectorScalar%bisSorted) then
        call output_line("Vector must be unsorted!",&
                         OU_CLASS_ERROR, OU_MODE_STD, "pperr_scalarBoundary2D")
        call sys_halt()
      end if
      cdataType = rvectorScalar%cdataType
    else
      cdataType = ST_DOUBLE
    end if

    ! If the boundary region is specified, call pperr_scalarBoundary2d_conf
    ! for that boundary region. Otherwise, call pperr_scalarBoundary2d_conf
    ! for all possible boundary regions and sum up the errors.
    if (present(rboundaryRegion)) then

      select case(cdataType)

      case (ST_DOUBLE)
        call pperr_scalarBoundary2d_conf (cerrortype, ccubType, derror,&
                                          rboundaryRegion, p_rdiscretisation,&
                                          rvectorScalar, ffunctionReference,&
                                          rcollection, ffunctionWeight)
      case DEFAULT
        call output_line("Single precision vectors currently not supported!",&
                         OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarBoundary2D")
        call sys_halt()
      end select

    else

      select case(cdataType)

      case (ST_DOUBLE)

        derror = 0.0_DP
        ! Create a boundary region for each boundary component and call
        ! the calculation routine for that.
        do ibdc = 1,boundary_igetNBoundComp(p_rdiscretisation%p_rboundary)
          call boundary_createRegion (p_rdiscretisation%p_rboundary, &
                                      ibdc, 0, rboundaryReg)
          call pperr_scalarBoundary2d_conf (cerrortype, ccubType, dlocalError,&
                                            rboundaryReg, p_rdiscretisation,&
                                            rvectorScalar, ffunctionReference,&
                                            rcollection, ffunctionWeight)
          derror = derror + dlocalError
        end do

      case DEFAULT
        call output_line("Single precision vectors currently not supported!",&
                         OU_CLASS_ERROR,OU_MODE_STD,"pperr_scalarBoundary2D")
        call sys_halt()
      end select

    end if

  end subroutine pperr_scalarBoundary2D

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarBoundary2d_conf (cerrortype, ccubType, derror,&
                                          rboundaryRegion, rdiscretisation,&
                                          rvectorScalar, ffunctionReference,&
                                          rcollection, ffunctionWeight)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(in) :: cerrortype

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType

  ! A t_boundaryRegion specifying the boundary region where
  ! to calculate.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation

  ! OPTIONAL: The FE solution vector. Represents a scalar FE function.
  ! If omitted, the function is assumed to be constantly =0.
  type(t_vectorScalar), intent(in), optional, target :: rvectorScalar

  ! OPTIONAL: A callback function that provides a coefficient in front
  ! of the FE function. If not specified, a value of 1 is assumed.
  include 'intf_refFunctionScBdr2D.inc'
  optional :: ffunctionReference

  ! OPTIONAL: A callback function that provides the weighting function
  ! by which the computed error is multipled.
  ! If not specified, the reference function is assumed to be =1!
  optional :: ffunctionWeight
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional
  ! information for callback routines.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(out) :: derror
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), allocatable :: IelementOrientation
    integer, dimension(:), allocatable :: IelementList
    real(DP), dimension(:,:), allocatable :: DedgePosition

    integer :: NELbdc,NVE,iel,i,k
    integer(I32) :: ctrafoType

    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Arrays for cubature points
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:), allocatable :: Dxi2D,Dpoints,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointPar
    real(DP), dimension(CUB_MAXCUBP) :: Domega1D
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    real(DP), dimension(NDIM2D,TRIA_MAXNVE) :: Dcoord
    integer :: ncubp,ipoint,ielementDistr
    integer(I32) :: celement
    integer(i32) :: icoordSystem
    real(DP) :: dlen


    ! Get some pointers and arrays for quicker access
    p_rtriangulation => rdiscretisation%p_rtriangulation

    call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoordinates)

    ! Number of elements on that boundary component?
    NELbdc = bdraux_getNELAtRegion(rboundaryRegion, p_rtriangulation)

    ! In a first step, we figure out the elements on the boundary and their
    ! orientation. Allocate arrays that are large enough to hold
    ! even all elements on the boundary if necessary.
    allocate(IelementList(NELbdc), IelementOrientation(NELbdc))

    ! Allocate an array saving the start- and end-parameter values
    ! of the edges on the boundary.
    allocate(DedgePosition(2,NELbdc))

    ! Get the parameter values of the 1D cubature formula
    ! as well as the number of cubature points ncubp
    call cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega1D)


    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do ielementDistr = 1, rdiscretisation%inumFESpaces

      ! Set type of elements for this distribution
      celement = rdiscretisation%RelementDistr(ielementDistr)%celement

      ! Calculate the list of elements adjacent to the boundary
      call bdraux_getElementsAtRegion(rboundaryRegion,&
          rdiscretisation, NELbdc, IelementList, IelementOrientation,&
          DedgePosition, celement, BDR_PAR_LENGTH)


      ! Check if element distribution is empty
      if (NELbdc .le. 0) cycle


      ! Map the 1D cubature points to the edges in 2D.
      allocate(Dxi2D(ncubp,NDIM2D+1,NELbdc))

      ! Get the type of the coordinate system
      icoordSystem = elem_igetCoordSystem(celement)

      ! Transform the coordinates
      do iel = 1,NELbdc
        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      end do

      ! Transpose the coordinate array such that we get coordinates
      ! we can work with.
      allocate(DpointsRef(NDIM2D+1,ncubp,NELbdc))
      do iel = 1,NELbdc
        do i = 1,ncubp
          do k = 1,ubound(DpointsRef,1)
            DpointsRef(k,i,iel) = Dxi2D(i,k,iel)
          end do
        end do
      end do

      ! Dxi2D is not needed anymore.
      deallocate(Dxi2D)

      ! If the reference function or the weighting function exist,
      ! calculate the coordinates of the points on world coordinates
      if (present(ffunctionReference) .or.&
          present(ffunctionWeight)) then

        ! We need the real coordinates of the points.
        allocate(Dpoints(NDIM2D,ncubp,NELbdc))

        ! We need the parameter values of the points.
        allocate (DpointPar(ncubp,NELbdc))

        ! Get the number of corner vertices of the element
        NVE = elem_igetNVE(celement)

        ! All elements with the same transformation
        ctrafoType = elem_igetTrafoType(celement)

        do iel = 1,NELbdc

          ! Get the points forming the element
          do ipoint = 1, NVE
            Dcoord(1,ipoint) = &
                p_DvertexCoordinates(1, p_IverticesAtElement(ipoint,IelementList(iel)))
            Dcoord(2,ipoint) = &
                p_DvertexCoordinates(2, p_IverticesAtElement(ipoint,IelementList(iel)))
          end do

          do ipoint = 1,ncubp

            ! Transform the cubature points
            call trafo_calcRealCoords (ctrafoType, Dcoord(:,1:NVE),&
                DpointsRef(:,ipoint,iel), Dpoints(:,ipoint,iel))

            ! Calculate the parameter values of the points on the boundary
            ! Dxi1D is in [-1,1] while the current edge has parmeter values
            ! [DedgePosition(1),DedgePosition(2)]. So do a linear
            ! transformation to transform Dxi1D into that interval, this
            ! gives the parameter values in length parametrisation
            call mprim_linearRescale(Dxi1D(ipoint,1), -1.0_DP, 1.0_DP,&
                DedgePosition(1,iel), DedgePosition(2,iel), DpointPar(ipoint,iel))

          end do
        end do

      end if

      ! So Dxi2 defines the coordinates on the reference element for all
      ! elements. Generally said, we have to evaluate the elements in these
      ! points now. That can be done by using fevl_evaluate_mult.
      !
      ! Which type of integral is to calculate? H1 or L2 or L1?
      select case (cerrortype)
      case (PPERR_L2ERROR)

        allocate (Dcoefficients(ncubp,NELbdc,3))

        ! If the FE function exists, evaluate it.
        if (present(rvectorScalar)) then
          do iel = 1,NELbdc
            ! Evaluate on the element, write results to Dcoefficients
            call fevl_evaluate_mult (DER_FUNC, Dcoefficients(:,iel,1), rvectorScalar, &
                                     IelementList(iel), DpointsRef(:,:,iel))
          end do
        else
          Dcoefficients(:,:,1) = 0.0_DP
        end if

        ! If the reference function exists, evaluate it.
        if (present(ffunctionReference)) then

          ! Evaluate the reference function on the boundary
          call ffunctionReference (DER_FUNC, rdiscretisation, DpointsRef,&
                                   Dpoints, rboundaryRegion%iboundCompIdx,&
                                   DpointPar, IelementList(1:NELbdc),&
                                   Dcoefficients(:,:,2), rcollection)
        else
          Dcoefficients(:,:,2) = 0.0_DP
        end if

        ! If the weighting function exists, evaluate it.
        if (present(ffunctionWeight)) then

          ! Evaluate the reference function on the boundary
          call ffunctionWeight (rdiscretisation, DpointsRef, Dpoints,&
                                rboundaryRegion%iboundCompIdx, DpointPar,&
                                IelementList(1:NELbdc), Dcoefficients(:,:,3), rcollection)
        else
          Dcoefficients(:,:,3) = 1.0_DP
        end if

        ! Linear combination to get the actual values in the cubature points.
        do iel = 1,NELbdc
          do ipoint = 1,ncubp
            Dcoefficients(ipoint,iel,1) = Dcoefficients(ipoint,iel,2)-Dcoefficients(ipoint,iel,1)
          end do
        end do

        ! Now, Dcoefficients contains in Dcoefficients(:,:,1) the term
        ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
        ! element. We are finally able to calculate the integral!
        ! That means, run over all the edges and sum up...
        ! (ok, if rvectorScalar is not specified, we have
        !  -u_h(x,y) in Dcoefficients(:,:,1), but as we take the square,
        !  it does not matter if we have u_h or -u_h there!)

        derror = 0.0_DP
        do iel = 1, NELbdc

          ! Get the length of the edge. Let us use the parameter values
          ! on the boundary for that purpose; this is a more general
          ! implementation than using simple lines as it will later
          ! support isoparametric elements.
          !
          ! The length of the current edge serves as a "determinant"
          ! in the cubature, so we have to divide it by 2 as an edge on
          ! the unit interval [-1,1] has length 2.
          dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))

          do ipoint = 1, ncubp
            derror = derror + dlen * Domega1D(ipoint) *&
                Dcoefficients(ipoint,iel,3) * (Dcoefficients(ipoint,iel,1)**2)
          end do
        end do

        deallocate(Dcoefficients)

      case (PPERR_L1ERROR)

        allocate (Dcoefficients(ncubp,NELbdc,3))

        ! If the FE function exists, evaluate it.
        if (present(rvectorScalar)) then
          do iel = 1,NELbdc
            ! Evaluate on the element, write results to Dcoefficients
            call fevl_evaluate_mult (DER_FUNC, Dcoefficients(:,iel,1), rvectorScalar, &
                                     IelementList(iel), DpointsRef(:,:,iel))
          end do
        else
          Dcoefficients(:,:,1) = 0.0_DP
        end if

        ! If the reference function exists, evaluate it.
        if (present(ffunctionReference)) then

          ! Evaluate the reference function on the boundary
          call ffunctionReference (DER_FUNC, rdiscretisation, DpointsRef,&
                                   Dpoints, rboundaryRegion%iboundCompIdx,&
                                   DpointPar, IelementList(1:NELbdc),&
                                   Dcoefficients(:,:,2), rcollection)
        else
          Dcoefficients(:,:,2) = 0.0_DP
        end if

        ! If the weighting function exists, evaluate it.
        if (present(ffunctionWeight)) then

          ! Evaluate the reference function on the boundary
          call ffunctionWeight (rdiscretisation, DpointsRef, Dpoints,&
                                rboundaryRegion%iboundCompIdx, DpointPar,&
                                IelementList(1:NELbdc), Dcoefficients(:,:,3), rcollection)
        else
          Dcoefficients(:,:,3) = 1.0_DP
        end if

        ! Linear combination to get the actual values in the cubature points.
        do iel = 1,NELbdc
          do ipoint = 1,ncubp
            Dcoefficients(ipoint,iel,1) = Dcoefficients(ipoint,iel,2)-Dcoefficients(ipoint,iel,1)
          end do
        end do

        ! Now, Dcoefficients contains in Dcoefficients(:,:,1) the term
        ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
        ! element. We are finally able to calculate the integral!
        ! That means, run over all the edges and sum up...
        ! (ok, if rvectorScalar is not specified, we have
        !  -u_h(x,y) in Dcoefficients(:,:,1), but as we take the square,
        !  it does not matter if we have u_h or -u_h there!)

        derror = 0.0_DP
        do iel = 1,NELbdc

          ! Get the length of the edge. Let us use the parameter values
          ! on the boundary for that purpose; this is a more general
          ! implementation than using simple lines as it will later
          ! support isoparametric elements.
          !
          ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on
          ! the unit interval [-1,1] has length 2.
          dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))

          do ipoint = 1,ncubp
            derror = derror + dlen * Domega1D(ipoint) *&
                Dcoefficients(ipoint,iel,3) * abs(Dcoefficients(ipoint,iel,1))
          end do
        end do

        deallocate(Dcoefficients)

      case (PPERR_H1ERROR)

        allocate (Dcoefficients(ncubp,NELbdc,5))
        Dcoefficients = 0.0_DP

        ! If the FE function exists, evaluate it.
        if (present(rvectorScalar)) then
          do iel = 1,NELbdc
            ! Evaluate on the element, write results to Dcoefficients.
            !
            ! X-derivative
            call fevl_evaluate_mult (DER_DERIV_X, Dcoefficients(:,iel,1), rvectorScalar, &
                                     IelementList(iel), DpointsRef(:,:,iel))

            ! Y-derivative
            call fevl_evaluate_mult (DER_DERIV_Y, Dcoefficients(:,iel,2), rvectorScalar, &
                                     IelementList(iel), DpointsRef(:,:,iel))
          end do
        end if

        ! If the reference function exists, evaluate it.
        if (present(ffunctionReference)) then

          ! Evaluate the reference function on the boundary
          !
          ! X-derivative
          call ffunctionReference (DER_DERIV_X, rdiscretisation, DpointsRef,&
                                   Dpoints, rboundaryRegion%iboundCompIdx,&
                                   DpointPar, IelementList(1:NELbdc),&
                                   Dcoefficients(:,:,3), rcollection)

          ! Y-derivative
          call ffunctionReference (DER_DERIV_Y, rdiscretisation, DpointsRef,&
                                   Dpoints, rboundaryRegion%iboundCompIdx,&
                                   DpointPar, IelementList(1:NELbdc),&
                                   Dcoefficients(:,:,4), rcollection)
        end if

        ! If the weighting function exists, evaluate it.
        if (present(ffunctionWeight)) then

          ! Evaluate the reference function on the boundary
          call ffunctionWeight (rdiscretisation, DpointsRef, Dpoints,&
                                rboundaryRegion%iboundCompIdx, DpointPar,&
                                IelementList(1:NELbdc), Dcoefficients(:,:,5), rcollection)
        else
          Dcoefficients(:,:,5) = 1.0_DP
        end if

        ! Linear combination to get the actual values in the cubature points.
        ! ||u-u_h||_H1 = int ( grad(u-u_h) * grad(u-u_h) )
        !              ~ sum grad_x(u-u_h)**2 + grad_y(u-u_h)
        do iel = 1,NELbdc
          do ipoint = 1,ncubp
            Dcoefficients(ipoint,iel,1) =&
                (Dcoefficients(ipoint,iel,1)-Dcoefficients(ipoint,iel,3))**2 + &
                (Dcoefficients(ipoint,iel,2)-Dcoefficients(ipoint,iel,4))**2
          end do
        end do

        ! Now, Dcoefficients contains in Dcoefficients(:,:,1) the term
        ! "grad(u(x,y))-grad(u_h(x,y))" -- in every cubature point on every
        ! element. We are finally able to calculate the integral!
        ! That means, run over all the edges and sum up...

        derror = 0.0_DP
        do iel = 1,NELbdc

          ! Get the length of the edge. Let us use the parameter values
          ! on the boundary for that purpose; this is a more general
          ! implementation than using simple lines as it will later
          ! support isoparametric elements.
          !
          ! The length of the current edge serves as a "determinant"
          ! in the cubature, so we have to divide it by 2 as an edge on
          ! the unit interval [-1,1] has length 2.
          dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))

          do ipoint = 1,ncubp
            derror = derror + dlen * Domega1D(ipoint) *&
                Dcoefficients(ipoint,iel,5) * (Dcoefficients(ipoint,iel,1)**2)
          end do
        end do

        deallocate(Dcoefficients)

      case (PPERR_MEANERROR)

        allocate (Dcoefficients(ncubp,NELbdc,3))

        ! If the FE function exists, evaluate it.
        if (present(rvectorScalar)) then
          do iel = 1,NELbdc
            ! Evaluate on the element, write results to Dcoefficients
            call fevl_evaluate_mult (DER_FUNC, Dcoefficients(:,iel,1), rvectorScalar, &
                                     IelementList(iel), DpointsRef(:,:,iel))
          end do
        else
          Dcoefficients(:,:,1) = 0.0_DP
        end if

        ! If the reference function exists, evaluate it.
        if (present(ffunctionReference)) then

          ! Evaluate the reference function on the boundary
          call ffunctionReference (DER_FUNC, rdiscretisation, DpointsRef,&
                                   Dpoints, rboundaryRegion%iboundCompIdx,&
                                   DpointPar, IelementList(1:NELbdc),&
                                   Dcoefficients(:,:,2), rcollection)
        else
          Dcoefficients(:,:,2) = 0.0_DP
        end if

        ! If the weighting function exists, evaluate it.
        if (present(ffunctionWeight)) then

          ! Evaluate the reference function on the boundary
          call ffunctionWeight (rdiscretisation, DpointsRef, Dpoints,&
                                rboundaryRegion%iboundCompIdx, DpointPar,&
                                IelementList(1:NELbdc), Dcoefficients(:,:,3), rcollection)
        else
          Dcoefficients(:,:,3) = 1.0_DP
        end if

        ! Linear combination to get the actual values in the cubature points.
        do iel = 1,NELbdc
          do ipoint = 1,ncubp
            Dcoefficients(ipoint,iel,1) = Dcoefficients(ipoint,iel,2)-Dcoefficients(ipoint,iel,1)
          end do
        end do

        ! Now, Dcoefficients contains in Dcoefficients(:,:,1) the term
        ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
        ! element. We are finally able to calculate the integral!
        ! That means, run over all the edges and sum up...
        ! (ok, if rvectorScalar is not specified, we have
        !  -u_h(x,y) in Dcoefficients(:,:,1), but as we take the square,
        !  it does not matter if we have u_h or -u_h there!)

        derror = 0.0_DP
        do iel = 1,NELbdc

          ! Get the length of the edge. Let us use the parameter values
          ! on the boundary for that purpose; this is a more general
          ! implementation than using simple lines as it will later
          ! support isoparametric elements.
          !
          ! The length of the current edge serves as a "determinant"
          ! in the cubature, so we have to divide it by 2 as an edge on
          ! the unit interval [-1,1] has length 2.
          dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))

          do ipoint = 1,ncubp
            derror = derror + dlen * Domega1D(ipoint) *&
                Dcoefficients(ipoint,iel,3) * Dcoefficients(ipoint,iel,1)
          end do
        end do

        deallocate(Dcoefficients)

      case default

        ! This case realises
        !
        !  int (w u) dx
        !
        ! with w being the weight and u the analytical function given by
        ! ffunctionReference. This function must be present such that it can
        ! be evaluated -- otherwise the user made a mistake in calling
        ! this routine.

        allocate (Dcoefficients(ncubp,NELbdc,2))

        ! If the reference function exists, evaluate it.
        if (present(ffunctionReference)) then

          ! Evaluate the reference function on the boundary
          call ffunctionReference (DER_FUNC, rdiscretisation, DpointsRef,&
                                   Dpoints, rboundaryRegion%iboundCompIdx,&
                                   DpointPar, IelementList(1:NELbdc),&
                                   Dcoefficients(:,:,1), rcollection)
        else
          call output_line('Reference function missing in user-defined error type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2d_conf')
          call sys_halt()
        end if

        ! If the weighting function exists, evaluate it.
        if (present(ffunctionWeight)) then

          ! Evaluate the reference function on the boundary
          call ffunctionWeight (rdiscretisation, DpointsRef, Dpoints,&
                                rboundaryRegion%iboundCompIdx, DpointPar,&
                                IelementList(1:NELbdc), Dcoefficients(:,:,2), rcollection)
        else
          Dcoefficients(:,:,2) = 1.0_DP
        end if

        ! Now, Dcoefficients contains in Dcoefficients(:,:,1) the term
        ! "f(u(x,y)-u_h(x,y)" -- in every cubature point on every
        ! element. We are finally able to calculate the integral!
        ! That means, run over all the edges and sum up...

        derror = 0.0_DP
        do iel = 1,NELbdc

          ! Get the length of the edge. Let us use the parameter values
          ! on the boundary for that purpose; this is a more general
          ! implementation than using simple lines as it will later
          ! support isoparametric elements.
          !
          ! The length of the current edge serves as a "determinant"
          ! in the cubature, so we have to divide it by 2 as an edge on
          ! the unit interval [-1,1] has length 2.
          dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))

          do ipoint = 1,ncubp
            derror = derror + dlen * Domega1D(ipoint) *&
                Dcoefficients(ipoint,iel,2) * Dcoefficients(ipoint,iel,1)
          end do
        end do

        deallocate(Dcoefficients)

      end select

      ! Release memory
      deallocate(DpointsRef)

      if (present(ffunctionReference) .or.&
          present(ffunctionWeight)) deallocate(Dpoints, DpointPar)

    end do

    ! Release memory
    deallocate(IelementList, IelementOrientation, DedgePosition)

  end subroutine pperr_scalarBoundary2d_conf

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarErrorEstimate (rvector, rvectorRef, ctype, derror,&
                                        rcubatureInfo, relementError,&
                                        rperfconfig)

!<description>
  ! This routine calculates the error of a given FE function in
  ! rvector and a reference vector given in rvectorRef. As an example,
  ! one can think of the consistent FE gradient and some recovered
  ! reference gradient, c.f. ZZ-technique.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution vector
    type(t_vectorScalar), intent(in), target :: rvector

    ! FE reference solution vector
    type(t_vectorScalar), intent(in), target :: rvectorRef

    ! Type of error to compute. A PPERR_xxERROR constant.
    ! PPERR_L2ERROR computes the L2-error, PPERR_L1ERROR the L1-error.
    integer, intent(in) :: ctype

    ! OPTIONAL: A scalar cubature information structure that gives additional information
    ! about how to set up the matrix (e.g. cubature formula). If not specified,
    ! default settings are used.
    type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated error on each element.
    type(t_vectorScalar), intent(inout), optional :: relementError
!</inputoutput>

!<output>
    ! The calculated error.
    real(DP), intent(out) :: derror
!</output>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorBlock,rvectorBlockRef
    type(t_blockDiscretisation) :: rDiscr,rDiscrRef

    ! Create block discretisations with one component
    if (associated(rvector%p_rspatialdiscr)) then
      call spdiscr_createBlockDiscrInd(rvector%p_rspatialdiscr, rDiscr)
    else
      call output_line('Vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarL2ErrorEstimate')
      call sys_halt()
    end if

    if (associated(rvectorRef%p_rspatialdiscr)) then
      call spdiscr_createBlockDiscrInd(rvectorRef%p_rspatialdiscr, rDiscrRef)
    else
      call output_line('Reference vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarL2ErrorEstimate')
      call sys_halt()
    end if

    ! Create block vectors with one block
    call lsysbl_createVecFromScalar(rvector, rvectorBlock, rDiscr)
    call lsysbl_createVecFromScalar(rvectorRef, rvectorBlockRef, rDiscrRef)

    ! Call block version
    call pperr_blockErrorEstimate(rvectorBlock, rvectorBlockRef, ctype,&
        derror, rcubatureInfo, relementError, rperfconfig)

    ! Release auxiliary block discretisations
    call spdiscr_releaseBlockDiscr(rDiscr)
    call spdiscr_releaseBlockDiscr(rDiscrRef)

    ! Release auxiliary block vectors
    call lsysbl_releaseVector(rvectorBlock)
    call lsysbl_releaseVector(rvectorBlockRef)

  end subroutine pperr_scalarErrorEstimate

  !****************************************************************************

!<subroutine>

  subroutine pperr_blockErrorEstimate (rvector, rvectorRef, ctype, derror,&
                                       rcubatureInfo, relementError,&
                                       rperfconfig)

!<description>
  ! This routine calculates the error of a given FE function in rvector
  ! and a reference vector given in rvectorRef. Both vectors must have the
  ! same number of blocks. As an example, one can think of the consistent
  ! FE gradient and some recovered reference gradient, c.f. ZZ-technique.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution block vector
    type(t_vectorBlock), intent(in), target :: rvector

    ! FE reference solution block vector
    type(t_vectorBlock), intent(in), target :: rvectorRef

    ! Type of error to compute. A PPERR_xxERROR constant.
    ! PPERR_L2ERROR computes the L2-error.
    ! PPERR_L1ERROR computes the L1-error.
    ! PPERR_H1ERROR computes the H1-error.
    integer, intent(in) :: ctype

    ! OPTIONAL: A scalar cubature information structure that gives additional information
    ! about how to set up the matrix (e.g. cubature formula). If not specified,
    ! default settings are used.
    type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated error on each element.
    type(t_vectorScalar), intent(inout), optional :: relementError
!</inputoutput>

!<output>
    ! The calculated error.
    real(DP), intent(out) :: derror
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationRef
    integer :: i,k,ielementDistr,ielementDistrRef,iblock,ICUBP,NVE
    integer :: IEL, IELmax, IELset,IELGlobal
    real(DP) :: OM,delementError

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
    integer :: indofTrial,indofTrialRef

    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation

    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! Jacobian determinant in all points
    real(DP), dimension(:,:), pointer :: p_Ddetj

    ! Pointer to the element error
    real(DP), dimension(:), pointer :: p_DelementError

    ! Number of elements in the current element distribution
    integer :: NEL,NELref

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients

    ! Type of transformation from the reference to the real element
    integer(I32) :: ctrafoType

    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Number of elements in a block. Normally =NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock

    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    type(t_evalElementSet) :: revalElementSet

    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial,IdofsTrialRef

    ! Current assembly block, cubature formula, element type,...
    integer :: iinfoBlock
    integer(I32) :: celement, celementRef, ccubature, ccubatureRef
    type(t_scalarCubatureInfo), target :: rlocalCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => pperr_perfconfig
    end if

    ! Get the correct discretisation structure for the solution vector
    p_rdiscretisation => rvector%p_rblockDiscr%RspatialDiscr(1)
    do iblock=2,rvector%nblocks
      call lsyssc_checkDiscretisation (rvector%RvectorBlock(iblock), p_rdiscretisation)
    end do

    ! The vectors must have the same number of blocks
    if (rvector%nblocks .ne. rvectorRef%nblocks) then
      call output_line("Vectors have different number of blocks!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockErrorEstimate")
      call sys_halt()
    end if

    ! The vector must be unsorted.
    do iblock=1,rvector%nblocks
      if (rvector%RvectorBlock(iblock)%bisSorted .or.&
          rvectorRef%RvectorBlock(iblock)%bisSorted) then
        call output_line("Vectors must be unsorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockErrorEstimate")
        call sys_halt()
      end if
    end do

    ! Do we have a cubature information structure?
    ! If we do not have it, create temporary one that
    ! defines how to do the cubature.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscretisation,&
          rlocalCubatureInfo,CUB_GEN_DEPR_EVAL)
      p_rcubatureInfo => rlocalCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! We only need the function values of basis functions
    Bder = .false.
    Bder(DER_FUNC) = .true.

    ! Get the discretisation from the reference vector
    p_rdiscretisationRef => rvectorRef%p_rblockDiscr%RspatialDiscr(1)

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscretisationRef%p_rtriangulation

    ! Set the current error to 0 and add the error contributions of each element to that.
    derror = 0.0_DP

    ! Get a pointer to the element error (if required)
    if (present(relementError)) then
      call lsyssc_getbase_double(relementError,p_DelementError)
      call lalg_clearVector (p_DelementError)
    end if

    ! Check that both discretisations have the same number of element distributions
    if (p_rdiscretisation%inumFESpaces .ne. &
        p_rdiscretisationRef%inumFESpaces) then
      call output_line("Number of element distributions mismatch!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockErrorEstimate")
      call sys_halt()
    end if

    ! Loop over the element blocks. Each defines a separate cubature formula
    ! and is connected to a combination of test and trial functions
    ! via the corresponding element distribution.
    do iinfoBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get typical information: Number of elements, element list,...
      call spdiscr_getStdDiscrInfo (iinfoBlock,p_rcubatureInfo,p_rdiscretisation,&
          ielementDistr,celement,ccubature,NEL,p_IelementList)
      call spdiscr_getStdDiscrInfo (iinfoBlock,p_rcubatureInfo,p_rdiscretisationRef,&
          ielementDistrRef,celementRef,ccubatureRef,NELref)

      ! Check if element distributions have different number of elements
      if (p_rdiscretisation%RelementDistr(ielementDistr)%NEL .ne. &
          p_rdiscretisation%RelementDistr(ielementDistrRef)%NEL) then
        call output_line("Number of elements in distributions mismatch!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockErrorEstimate")
        call sys_halt()
      end if

      ! Cancel if this element list is empty.
      if (NEL .le. 0) cycle

      ! For saving some memory in smaller discretisations, we calculate
      ! the number of elements per block. For smaller triangulations,
      ! this is NEL. If there are too many elements, it is at most
      ! NELEMSIM. This is only used for allocating some arrays.
      nelementsPerBlock = min(p_rperfconfig%NELEMSIM,NEL)

      ! Get the number of local DOF`s for trial functions
      indofTrial    = elem_igetNDofLoc(celement)
      indofTrialRef = elem_igetNDofLoc(celementRef)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(celementRef)

      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(ccubatureRef, ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(celementRef)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do

      ! Allocate memory for the DOF`s of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))
      allocate(IdofsTrialRef(indofTrialRef,nelementsPerBlock))

      ! Allocate memory for the coefficients, that is, two times the spatial dimension
      allocate(Dcoefficients(ncubp,nelementsPerBlock,2*rvector%nblocks))

      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(celement)
      cevaluationTag = ior(cevaluationTag,&
                       elem_getEvaluationTag(celementRef))

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, nelementsPerBlock

        ! We always handle NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most NELEMSIM
        ! elements simultaneously.

        IELmax = min(NEL,IELset-1+nelementsPerBlock)

        ! Calculate the global DOF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
        call dof_locGlobMapping_mult(p_rdiscretisationRef, p_IelementList(IELset:IELmax), &
                                     IdofsTrialRef)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we do not have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Depending on ctype, choose now the error to compute.
        select case (ctype)
        case (PPERR_L1ERROR)

          ! L1-error uses only the values of the function.

          ! Calculate the values of the FE solution vector and the reference solution
          ! vector in the cubature points: u_h(x,y) and u_ref(x,y)
          ! Save the result to Dcoefficients(:,:,2*iblock-1) and
          ! Dcoefficients(:,:,2*iblock)

          do iblock=1,rvector%nblocks

            ! solution vector
            call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                    celement, IdofsTrial, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock))

            ! solution reference vector
            call fevl_evaluate_sim3 (rvectorRef%RvectorBlock(iblock), revalElementSet,&
                    celementRef, IdofsTrialRef, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock-1))

          end do

          ! Subtraction of Dcoefficients(:,:,2*iblock-1) from Dcoefficients(:,:,2*iblock)
          ! and summing over all iblock=1,..,nblocks gives the error
          ! $u_h(cubature pt.) - u_ref(cubature pt.)$

          ! Loop through elements in the set and for each element,
          ! loop through the DOF`s and cubature points to calculate the
          ! integral: int_Omega (u_h-u_ref,u_h-u_ref) dx

          do IEL=1,IELmax-IELset+1

            ! Initialise element error by 0
            delementError = 0.0_DP

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!

              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

              ! L1-error is:   int_... abs(u_h-u_ref) dx

              do iblock=1,rvector%nblocks
                delementError = delementError + &
                    OM * abs(Dcoefficients(icubp,IEL,2*iblock-1)-&
                            Dcoefficients(icubp,IEL,2*iblock))
              end do

            end do ! ICUBP

            ! Apply to global error
            derror = derror + delementError

            ! Store in element error (if required)
            if (present(relementError)) then
              IELGlobal = p_IelementList(IELset+IEL-1)
              p_DelementError(IELGlobal) = delementError
            end if

          end do ! IEL

        case (PPERR_L2ERROR)

          ! L2-error uses only the values of the function.

          ! Calculate the values of the FE solution vector and the reference solution
          ! vector in the cubature points: u_h(x,y) and u_ref(x,y)
          ! Save the result to Dcoefficients(:,:,2*iblock-1) and
          ! Dcoefficients(:,:,2*iblock)

          do iblock=1,rvector%nblocks

            ! solution vector
            call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                    celement, IdofsTrial, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock))

            ! solution reference vector
            call fevl_evaluate_sim3 (rvectorRef%RvectorBlock(iblock), revalElementSet,&
                    celementRef, IdofsTrialRef, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock-1))

          end do

          ! Subtraction of Dcoefficients(:,:,2*iblock-1) from Dcoefficients(:,:,2*iblock)
          ! and summing over all iblock=1,..,nblocks gives the error
          ! $u_h(cubature pt.) - u_ref(cubature pt.)$

          ! Loop through elements in the set and for each element,
          ! loop through the DOF`s and cubature points to calculate the
          ! integral: int_Omega (u_h-u_ref,u_h-u_ref) dx

          do IEL=1,IELmax-IELset+1

            ! Initialise element error by 0
            delementError = 0.0_DP

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!

              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

              ! L2-error is:   int_... (u_h-u_ref)*(u_h-u_ref) dx

              do iblock=1,rvector%nblocks
                delementError = delementError + &
                    OM * (Dcoefficients(icubp,IEL,2*iblock-1)-&
                          Dcoefficients(icubp,IEL,2*iblock))**2
              end do

            end do ! ICUBP

            ! Apply to global error
            derror = derror + delementError

            ! Store in element error (if required)
            if (present(relementError)) then
              IELGlobal = p_IelementList(IELset+IEL-1)
              p_DelementError(IELGlobal) = sqrt(delementError)
            end if

          end do ! IEL

        case DEFAULT

          call output_line("Requested error estimate not implemented!",&
              OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockErrorEstimate")
          call sys_halt()

        end select

      end do ! IELset

      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial,IdofsTrialRef)

    end do ! ielementDistr

    if (ctype .ne. PPERR_L1ERROR) then
      ! derror is ||error||^2, so take the square root at last.
      derror = sqrt(derror)
    end if

    ! Release the cubature structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rlocalCubatureInfo)
    end if

  end subroutine pperr_blockErrorEstimate

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarStandardDeviation (rvector, ddeviation,&
                                            relementDeviation, rperfconfig)

!<description>
  ! This routine calculates the standard deviation
  !
  ! <tex> $$ \sigma=\sqrt{\int_\Omega r^2 u dx} $$ </tex>
  !
  ! of a given FE function $u$ in rvector, whereby
  !
  ! <tex> $$ r^2=(x-\hat x)^2 + (y-\hat y)^2 + (z-\hat z)^2 $$ </tex>
  !
  ! and each component is computed from the following relation
  !
  ! <tex> $$ \hat x=\int_\Omega x u dx $$ </tex>
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution vector
    type(t_vectorScalar), intent(in), target :: rvector

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated deviation on each element.
    type(t_vectorScalar), intent(inout), optional :: relementDeviation
!</inputoutput>

!<output>
    ! The calculated standard deviation.
    real(DP), intent(out) :: ddeviation
!</output>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorBlock
    type(t_blockDiscretisation) :: rDiscr

    ! Create block discretisations with one component
    if (associated(rvector%p_rspatialdiscr)) then
      call spdiscr_createBlockDiscrInd(rvector%p_rspatialdiscr, rDiscr)
    else
      call output_line('Vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarStandardDeviation')
      call sys_halt()
    end if

    ! Create block vectors with one block
    call lsysbl_createVecFromScalar(rvector, rvectorBlock, rDiscr)

    ! Call block version
    call pperr_blockStandardDeviation(rvectorBlock, ddeviation,&
        relementDeviation, rperfconfig)

    ! Release auxiliary block discretisations
    call spdiscr_releaseBlockDiscr(rDiscr)

    ! Release auxiliary block vectors
    call lsysbl_releaseVector(rvectorBlock)

  end subroutine pperr_scalarStandardDeviation

  !****************************************************************************

!<subroutine>

  subroutine pperr_blockStandardDeviation (rvector, ddeviation,&
                                           relementDeviation, rperfconfig)

!<description>
  ! This routine calculates the standard deviation
  !
  ! <tex> $$ \sigma=\sqrt{\int_\Omega r^2 u dx} $$ </tex>
  !
  ! of a given FE function $u$ in rvector, whereby
  !
  ! <tex> $$ r^2=(x-\hat x)^2 + (y-\hat y)^2 + (z-\hat z)^2 $$ </tex>
  !
  ! and each component is computed from the following relation
  !
  ! <tex> $$ \hat x=\int_\Omega x u dx $$ </tex>
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution block vector
    type(t_vectorBlock), intent(in), target :: rvector

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated deviation on each element.
    type(t_vectorScalar), intent(inout), optional :: relementDeviation
!</inputoutput>

!<output>
    ! The calculated deviation.
    real(DP), intent(out) :: ddeviation
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    integer :: i,k,ielementDistr,iblock,ICUBP,NVE,idim
    integer :: IEL, IELmax, IELset,IELGlobal
    real(DP) :: OM,delementDeviation

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
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! An array receiving the coordinates of cubature points on
    ! the real element for all elements in a set.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsReal

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj

    ! Pointer to the element deviation
    real(DP), dimension(:), pointer :: p_DelementDeviation

    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution

    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients

    ! Mathematical expectation of the center of mass
    real(DP), dimension(NDIM3D) :: DmassCenter

    ! Type of transformation from the reference to the real element
    integer(I32) :: ctrafoType

    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Number of elements in a block. Normally =NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock

    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    type(t_evalElementSet) :: revalElementSet

    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => pperr_perfconfig
    end if

    ! Get the correct discretisation structure for the solution vector
    p_rdiscretisation => rvector%p_rblockDiscr%RspatialDiscr(1)
    do iblock=2,rvector%nblocks
      call lsyssc_checkDiscretisation (rvector%RvectorBlock(iblock), p_rdiscretisation)
    end do

    if (.not. associated(p_rdiscretisation)) then
      call output_line("No discretisation structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockStandardDeviation")
      call sys_halt()
    end if

    ! The vector must be unsorted.
    do iblock=1,rvector%nblocks
      if (rvector%RvectorBlock(iblock)%bisSorted) then
        call output_line("Vectors must be unsorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"pperr_blockStandardDeviation")
        call sys_halt()
      end if
    end do

    ! We only need the function values of basis functions
    Bder = .false.
    Bder(DER_FUNC) = .true.

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscretisation%p_rtriangulation

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! Set the mathematical expectation of the center of mass to 0
    DmassCenter = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do ielementDistr = 1,p_rdiscretisation%inumFESpaces

      ! Activate the current element distribution
      p_relementDistribution => p_rdiscretisation%RelementDistr(ielementDistr)

      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF`s for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)

      ! Initialise the cubature formula.
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistribution%ccubTypeEval,&
          ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do

      ! Allocate memory for the DOF`s of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients.
      allocate(Dcoefficients(ncubp,nelementsPerBlock,rvector%nblocks))

      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)

      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, p_rperfconfig%NELEMSIM

        ! We always handle NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most NELEMSIM
        ! elements simultaneously.

        IELmax = min(NEL,IELset-1+p_rperfconfig%NELEMSIM)

        ! Calculate the global DOF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj
        p_DcubPtsReal => revalElementSet%p_DpointsReal

        ! In the next loop, we do not have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Standard deviation uses only the values of the function.

        ! Calculate the values of the FE solution vector in the cubature
        ! points u_h(x,y) and save the result to Dcoefficients(:,:,iblock)

        do iblock=1,rvector%nblocks

          ! solution vector
          call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1,iblock))

        end do

        ! Calculate the mathematical expectation of the center of mass
        ! $\hat x_h=\int_\Omega x u_h dx$

        ! Loop through elements in the set and for each element,
        ! loop through the DOF`s and cubature points to calculate the
        ! integral: int_Omega x*u_h dx, for x,y and z

        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!

            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Mathematical expectation of the center of mass is:
            ! int_... x*u_h dx

            do iblock=1,rvector%nblocks
              do idim=1,p_rdiscretisation%ndimension
                DmassCenter(idim) = DmassCenter(idim) + &
                    OM * p_DcubPtsReal(idim,icubp,IEL) * &
                         Dcoefficients(icubp,IEL,iblock)
              end do
            end do

          end do ! ICUBP

        end do ! IEL

      end do ! IELset

      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)

    end do ! ielementDistr

    ! Ok, we have the mathematical expectation of the center of mass.
    ! Let us compute the standard deviation.

    Ddeviation = 0.0_DP

    ! Get a pointer to the element deviation (if required)
    if (present(relementDeviation)) then
      call lsyssc_getbase_double(relementDeviation,p_DelementDeviation)
      call lalg_clearVector (p_DelementDeviation)
    end if

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Make sure that we have determinants.
    cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do ielementDistr = 1,p_rdiscretisation%inumFESpaces

      ! Activate the current element distribution
      p_relementDistribution => p_rdiscretisation%RelementDistr(ielementDistr)

      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF`s for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)

      ! Initialise the cubature formula.
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistribution%ccubTypeEval,&
          ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do

      ! Allocate memory for the DOF`s of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients.
      allocate(Dcoefficients(ncubp,nelementsPerBlock,rvector%nblocks))

      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)

      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, p_rperfconfig%NELEMSIM

        ! We always handle NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most NELEMSIM
        ! elements simultaneously.

        IELmax = min(NEL,IELset-1+p_rperfconfig%NELEMSIM)

        ! Calculate the global DOF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj
        p_DcubPtsReal => revalElementSet%p_DpointsReal

        ! In the next loop, we do not have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Standard deviation uses only the values of the function.

        ! Calculate the values of the FE solution vector in the cubature
        ! points u_h(x,y) and save the result to Dcoefficients(:,:,iblock)

        do iblock=1,rvector%nblocks

          ! solution vector
          call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1,iblock))

        end do

        ! Calculate the standard deviation
        ! $\int_\Omega ((x-\hat x_h)^2 + (y-\hat y_h)^2 + (z-\hat z_h)^2) u_h dx$

        ! Loop through elements in the set and for each element,
        ! loop through the DOF`s and cubature points to calculate the
        ! integral: int_Omega (x-\hat x_h)*(x-\hat x_h)*u_h dx
        ! and sum up for x,y and z

        do IEL=1,IELmax-IELset+1

          ! Initialise element deviation by 0
          delementDeviation = 0.0_DP

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!

            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Standard deviation is: int_... (x-\hat x_h)*(x-\hat x_h)*u_h dx
            ! summed up for all x,y and z

            do iblock=1,rvector%nblocks
              do idim=1,p_rdiscretisation%ndimension
                delementDeviation = delementDeviation + &
                    abs(OM * Dcoefficients(icubp,IEL,iblock)) * &
                       (p_DcubPtsReal(idim,icubp,IEL)-DmassCenter(idim))**2
              end do
            end do

          end do ! ICUBP

          ! Apply to global deviation
          ddeviation = ddeviation + delementDeviation

          ! Store in element deviation (if required)
          if (present(relementDeviation)) then
            IELGlobal = p_IelementList(IELset+IEL-1)
            p_DelementDeviation(IELGlobal) = sqrt(delementDeviation)
          end if

        end do ! IEL

      end do ! IELset

      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)

    end do ! ielementDistr

    ! ddeviation is ||deviation||^2, so take the square root at last.
    if (ddeviation .ge. huge(ddeviation)) then
      ddeviation = 0._DP
    else
      ddeviation = sqrt(ddeviation)
    end if

  end subroutine pperr_blockStandardDeviation

end module pprocerror

!##############################################################################
!# ****************************************************************************
!# <name> errorestimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures and subroutines
!# for error estimation.
!#
!# 1.) errest_initErrorEstimator
!#     -> Initialize the error estimator from parameter file
!#
!# 2.) errest_releaseErrorEstimator
!#     -> Release the error estimator and all of its structures
!#
!# 3.) errest_clearErrorEstimator
!#     -> Clear all temporal vectors from the error estimator
!#
!# 4.) errest_calcH1Error
!#     -> Estimate the solution error for a scalar indicator variable
!#        in the H1-semi norm, that is, the error of the solution gradient
!#
!# 5.) errest_calcFirstDiffIndicator
!#     -> Compute the first-difference indicator by Ill
!#
!# 6.) errest_calcSecondDiffIndicator
!#     -> Compute the second-difference indicator by Loehner
!#
!# 7.) errest_calcProtectionLayers
!#     -> Adjust the error indicator to include some protection layers
!#
!# 8.) errest_calcGridIndicator
!#     -> Compute the grid indicator based on various strategies
!# </purpose>
!##############################################################################

module errorestimation

  use fsystem
  use genoutput
  use geometryaux
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocerror
  use pprocgradients
  use spatialdiscretisation
  use stdoperators
  use storage
  use triangulation

  implicit none

  private
  public :: t_errorEstimator
  public :: t_gradientRecovery
  public :: t_differenceIndicator
  public :: errest_initErrorEstimator
  public :: errest_releaseErrorEstimator
  public :: errest_clearErrorEstimator
  public :: errest_calcH1Error
  public :: errest_calcFirstDiffIndicator
  public :: errest_calcSecondDiffIndicator
  public :: errest_calcProtectionLayers
  public :: errest_calcGridIndicator

!<constants>

!<constantblock description="Global constants for error estimator">

  ! Corrected L2-projection
  integer, parameter, public :: ERREST_CL2PROJECTION = -1

  ! Corrected superconvergent patch recovery (vertex-based)
  integer, parameter, public :: ERREST_CSPR_VERTEX   = -2

  ! Corrected superconvergent patch recovery (element-based)
  integer, parameter, public :: ERREST_CSPR_ELEMENT  = -3

  ! Corrected superconvergent patch recovery (face-based)
  integer, parameter, public :: ERREST_CSPR_FACE     = -4

  ! L2-projection
  integer, parameter, public :: ERREST_L2PROJECTION  = 1

  ! Superconvergent patch recovery (vertex-based)
  integer, parameter, public :: ERREST_SPR_VERTEX    = 2

  ! Superconvergent patch recovery (element-based)
  integer, parameter, public :: ERREST_SPR_ELEMENT   = 3

  ! Superconvergent patch recovery (face-based)
  integer, parameter, public :: ERREST_SPR_FACE      = 4
 
  ! Limited averaging gradient recovery
  integer, parameter, public :: ERREST_LIMAVR        = 5

  ! First-difference indicator (by Ill)
  integer, parameter, public :: ERREST_FIRSTDIFF     = 6

  ! Second-difference indicator (by Löhner)
  integer, parameter, public :: ERREST_SECONDDIFF    = 7

!</constantblock>

!<constantblock description="Global constants for error redistribution">

  ! Use error 'as is'
  integer, parameter, public :: ERREST_ASIS          = 0

  ! Equidistribution of error
  integer, parameter, public :: ERREST_EQUIDIST      = 1

  ! Logarithmic equidistribution of error
  integer, parameter, public :: ERREST_LOGEQUIDIST   = 2

  ! Fixed-rate strategy
  integer, parameter, public :: ERREST_FIXEDRATE     = 3

  ! Automatic treshold based on RMS
  integer, parameter, public :: ERREST_AUTORMS       = 4

!</constantblock>

!</constants>

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  
!<types>
  
!<typeblock>
  
  ! This data structure is used to compute the H1-error
  ! by means of gradient recovery techniques
  type t_gradientRecovery

    ! Type of gradient recovery
    integer :: ierrorestimator       = 0

    ! Type of element
    integer :: ieltype               = 0

    ! Number of error variable
    integer :: ierrorvariable        = 0

    ! Scalar error variable
    type(t_vectorScalar) :: rerrorVariable

    ! Consistent finite element gradient
    type(t_vectorBlock) :: rgradient

    ! Recovered reference gradient
    type(t_vectorBlock) :: rgradientRef

    ! Discretization structure for consistent gradient
    type(t_blockDiscretisation) :: rdiscrBlock

    ! Discretizaiton structure for recovered gradient
    type(t_blockDiscretisation) :: rdiscrBlockRef

    ! Pointer to the triangulation structure
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! Pointer to the boundary structure
    type(t_boundary), pointer :: p_rboundary => null()

  end type t_gradientRecovery
!</typeblock>

!<typeblock>
  
  ! This data structure is used to compute first/second difference indicator
  type t_differenceIndicator

    ! Type of difference indicator
    integer :: ierrorestimator = 0

    ! Number of error variable
    integer :: ierrorvariable  = 0

    ! Tolerance for noise filter
    real(DP) :: dnoiseFilter   = 0._DP

    ! Absolute tolerance for filter
    real(DP) :: dabsFilter     = 0._DP

    ! Scalar indicator variable
    type(t_vectorScalar) :: rerrorVariable

    ! Pointer to the triangulation structure
    type(t_triangulation), pointer :: p_rtriangulation => null()
  end type t_differenceIndicator
!</typeblock>

!<typeblock>

  ! This data structure contains the complete error estimator
  type t_errorEstimator

    ! Type of error estimator
    integer :: ierrorestimator = 0

    ! Type of grid-indicator
    integer :: igridindicator  = 0
    
    ! Number of protection layers
    integer :: nprotectlayers  = 0
    
    ! Pointer to gradient recovery
    type(t_gradientRecovery), pointer :: p_rgradientRecovery => null()

    ! Pointer to difference indicator
    type(t_differenceIndicator), pointer :: p_rdifferenceIndicator => null()

  end type t_errorEstimator
!</typeblock>

!</types>
  
contains

  !*****************************************************************************

!<subroutine>

  subroutine errest_initErrorEstimator(rerrorEstimator, rparlist, ssection)

!<description>
    ! This subroutine initializes the error estimator
    ! with the values supplied by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! name of the section
    character(LEN=*), intent(IN) :: ssection
!</input>

!<output>
    ! error estimator
    type(t_errorEstimator), intent(OUT) :: rerrorEstimator
!</output>
!</subroutine>

    ! Set the type of error estimator
    call parlst_getvalue_int(rparlist, ssection, "ierrorestimator",&
        rerrorEstimator%ierrorestimator, ERREST_SECONDDIFF)
    
    ! Set the type of grid indicator
    call parlst_getvalue_int(rparlist, ssection, "igridindicator",&
        rerrorEstimator%igridindicator, ERREST_ASIS)
    
    ! Set the number of protection layers
    call parlst_getvalue_int(rparlist, ssection, "nprotectlayers",&
        rerrorEstimator%nprotectlayers, 0)
    

    ! What kind of error estimator are we?
    select case(rerrorEstimator%ierrorestimator)

    case (ERREST_CSPR_FACE: ERREST_LIMAVR)

      ! Recovery-based error estimator
      allocate(rerrorEstimator%p_rgradientRecovery)

      ! Set the type of gradient recovery
      rerrorEstimator%p_rgradientRecovery%ierrorestimator =&
          rerrorEstimator%ierrorestimator

      ! Set the error variable
      call parlst_getvalue_int(rparlist, ssection, "ierrorvariable",&
          rerrorEstimator%p_rgradientRecovery%ierrorvariable, 1)

      ! Set the type of element
      call parlst_getvalue_int(rparlist, ssection, "ieltype",&
          rerrorEstimator%p_rgradientRecovery%ieltype)


    case (ERREST_FIRSTDIFF:ERREST_SECONDDIFF)

      ! First-/Second-difference indicator
      allocate(rerrorEstimator%p_rdifferenceIndicator)
      
      ! Set the type of difference indicator
      rerrorEstimator%p_rdifferenceIndicator%ierrorestimator =&
          rerrorEstimator%ierrorestimator

      ! Set the error variable
      call parlst_getvalue_int(rparlist, ssection, "ierrorvariable",&
          rerrorEstimator%p_rdifferenceIndicator%ierrorvariable, 1)

      ! Set the value for the noise filter
      call parlst_getvalue_double(rparlist, ssection, "dnoiseFilter",&
          rerrorEstimator%p_rdifferenceIndicator%dnoiseFilter, 5e-3_DP)

      ! Set the value for the absolute filter
      call parlst_getvalue_double(rparlist, ssection, "dabsFilter",&
          rerrorEstimator%p_rdifferenceIndicator%dabsFilter, 1e-6_DP)
      
    case DEFAULT
      call output_line('Invalid type of error estimator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'errest_initErrorEstimator')
      call sys_halt()
    end select

  end subroutine errest_initErrorEstimator

  !*****************************************************************************

!<subroutine>

  subroutine errest_releaseErrorEstimator(rerrorEstimator)

!<description>
    ! This subroutine releases the error estimator
!</description>

!<inputoutput>
    ! error estimator
    type(t_errorEstimator), intent(INOUT) :: rerrorEstimator
!</inputoutput>
!</subroutine>
    
    ! Release temporal data from error estimator
    call errest_clearErrorEstimator(rerrorEstimator)
    
    ! Reset data
    rerrorEstimator%ierrorestimator = 0
    rerrorEstimator%igridindicator  = 0
    rerrorEstimator%nprotectlayers  = 0

    ! Release substructures
    if (associated(rerrorEstimator%p_rgradientRecovery)) then

      ! Reset data
      rerrorEstimator%p_rgradientRecovery%ierrorestimator = 0
      rerrorEstimator%p_rgradientRecovery%ieltype         = 0
      rerrorEstimator%p_rgradientRecovery%ierrorvariable  = 0

      ! Deallocate memory
      deallocate(rerrorEstimator%p_rgradientRecovery)
      nullify(rerrorEstimator%p_rgradientRecovery)
    end if

    if (associated(rerrorEstimator%p_rdifferenceIndicator)) then

      ! Reset data
      rerrorEstimator%p_rdifferenceIndicator%dnoiseFilter    = 0._DP
      rerrorEstimator%p_rdifferenceIndicator%dabsFilter      = 0._DP
      rerrorEstimator%p_rdifferenceIndicator%ierrorestimator = 0
      rerrorEstimator%p_rdifferenceIndicator%ierrorvariable  = 0

      ! Deallocate memory
      deallocate(rerrorEstimator%p_rdifferenceIndicator)
      nullify(rerrorEstimator%p_rdifferenceIndicator)
    end if

  end subroutine errest_releaseErrorEstimator

  !*****************************************************************************

!<subroutine>

  subroutine errest_clearErrorEstimator(rerrorEstimator)

!<description>
    ! This subroutine clears temporal data of the error estimator
!</description>

!<inputoutput>
    ! error estimator
    type(t_errorEstimator), intent(INOUT) :: rerrorEstimator
!</inputoutput>
!</subroutine>

    ! Clear temporal data for gradient recovery
    if (associated(rerrorEstimator%p_rgradientRecovery)) then
      
      ! Release discretisation structure
      call spdiscr_releaseBlockDiscr(rerrorEstimator%p_rgradientRecovery%rdiscrBlock)
      call spdiscr_releaseBlockDiscr(rerrorEstimator%p_rgradientRecovery%rdiscrBlockRef)
      
      ! Release vectors
      call lsyssc_releaseVector(rerrorEstimator%p_rgradientRecovery%rerrorVariable)
      call lsysbl_releaseVector(rerrorEstimator%p_rgradientRecovery%rgradient)
      call lsysbl_releaseVector(rerrorEstimator%p_rgradientRecovery%rgradientRef)

      ! Nullify pointer
      nullify(rerrorEstimator%p_rgradientRecovery%p_rtriangulation)
      nullify(rerrorEstimator%p_rgradientRecovery%p_rboundary)
    end if

    ! Clear temporal data for difference indicator
    if (associated(rerrorEstimator%p_rdifferenceIndicator)) then

      ! Release vectors
      call lsyssc_releaseVector(rerrorEstimator%p_rdifferenceIndicator%rerrorVariable)

      ! Nullify pointer
      nullify(rerrorEstimator%p_rdifferenceIndicator%p_rtriangulation)
    end if

  end subroutine errest_clearErrorEstimator

  !*****************************************************************************

!<subroutine>

  subroutine errest_calcH1Error(rgradientRecovery, rerrorH1, derrorH1Opt)

!<description>
    ! This subroutine estimates the error in the H1-semi-norm by
    ! means of gradient recovery techniques.
!</description>
    
!<inputoutput>
    ! gradient recovery structure
    type(t_gradientRecovery), intent(INOUT) :: rgradientRecovery
!</inputoutput>

!<output>
    ! local H1-error (per element)
    type(t_vectorScalar), intent(OUT) :: rerrorH1

    ! OPTIONAL: global H1-error
    real(DP), intent(OUT), optional :: derrorH1Opt
!</output>
!</subroutine>
    
    ! How many spatial dimensions are we?
    select case (rgradientRecovery%p_rtriangulation%ndim)
    case (NDIM1D)
      call calcH1Error1D(rgradientRecovery, rerrorH1, derrorH1Opt)
      
    case (NDIM2D)
      call calcH1Error2D(rgradientRecovery, rerrorH1, derrorH1Opt)
      
    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'errest_calcH1Error')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Calculate the H1 error in 1D
    
    subroutine calcH1Error1D(rgradientRecovery, rerrorH1, derrorH1Opt)

      type(t_gradientRecovery), intent(INOUT) :: rgradientRecovery
      type(t_vectorScalar), intent(OUT) :: rerrorH1
      real(DP), intent(OUT), optional :: derrorH1Opt

      ! local variables
      type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
      real(DP) :: derrorH1

      ! Create new vector for local H1-error
      call lsyssc_createVector(rerrorH1, rgradientRecovery%p_rtriangulation%NEL, .true.)

      ! What error estimator should be used?
      select case(rgradientRecovery%ierrorestimator)
      case (ERREST_CSPR_FACE:ERREST_SPR_FACE)
        !------------------------------------------------------------------------
        ! Superconvergent patch recovery
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rgradientRecovery%rerrorVariable%p_rspatialDiscr
        
        ! Initialise block discretisations
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlock, 1 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlockRef, 1 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        
        ! What kind of element type is used for the FE solution?
        select case(rgradientRecovery%ieltype)
        case(-1,1,11)
          ! Initialise spatial discretisations for gradient with P0-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E000_1D, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          
          ! Initialise spatial discretisations for reference gradient with P1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          
        case(2,13)
          ! Initialise spatial discretisations for gradient with P1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          
          ! Initialise spatial discretisations for reference gradient with P2-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
                    
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error1D')
          call sys_halt()
        end select
        
        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlock,&
                                           rgradientRecovery%rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlockRef,&
                                           rgradientRecovery%rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
                                 rgradientRecovery%rgradient)
        
        ! How should the smoothed gradient be recovered?
        select case(rgradientRecovery%ierrorestimator)
        case (ERREST_L2PROJECTION)
          call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
                                   rgradientRecovery%rgradientRef)
          
        case (ERREST_SPR_VERTEX)
          call ppgrd_calcGradSuperPatchRecov (rgradientRecovery%rerrorVariable,&
                                              rgradientRecovery%rgradientRef, PPGRD_NODEPATCH)
          
        case (ERREST_SPR_ELEMENT)
          call ppgrd_calcGradSuperPatchRecov (rgradientRecovery%rerrorVariable,&
                                              rgradientRecovery%rgradientRef, PPGRD_ELEMPATCH)
          
        case (ERREST_SPR_FACE)
          call ppgrd_calcGradSuperPatchRecov (rgradientRecovery%rerrorVariable,&
                                              rgradientRecovery%rgradientRef, PPGRD_FACEPATCH)
        end select
        
        ! Compute estimated gradient error
        call pperr_blockErrorEstimate(rgradientRecovery%rgradient,&
                                      rgradientRecovery%rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
        
      case (ERREST_LIMAVR)
        !------------------------------------------------------------------------
        ! Limited gradient averaging
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rgradientRecovery%rerrorVariable%p_rspatialDiscr
        
        ! Initialise block discretisations
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlock, 1 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlockRef, 1 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        
        ! What kind of element type is used for the FE solution?
        select case(rgradientRecovery%ieltype)
        case(-1,1,11)
          ! Initialise spatial discretisations for gradient with P0-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E000_1D, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          
          ! Initialise spatial discretisations for reference gradient with P1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
                    
        case DEFAULT
          call output_line('Unsupported element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error1D')
          call sys_halt()
        end select
        
        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlock,&
                                           rgradientRecovery%rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlockRef,&
                                           rgradientRecovery%rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
                                 rgradientRecovery%rgradient)
        
        ! Recover limited averaged gradient
        call ppgrd_calcGradLimAvgP1Q1cnf (rgradientRecovery%rerrorVariable,&
                                          rgradientRecovery%rgradientRef)
        
        ! Compute gradient error
        call pperr_blockErrorEstimate(rgradientRecovery%rgradient,&
                                      rgradientRecovery%rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
      case DEFAULT
        call output_line('Invalid type of error estimator!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error1D')
        call sys_halt()
      end select
      
      ! Return global gradient error?
      if (present(derrorH1Opt)) derrorH1Opt = derrorH1
    end subroutine calcH1Error1D


    !**************************************************************
    ! Calculate the H1 error in 2D
    
    subroutine calcH1Error2D(rgradientRecovery, rerrorH1, derrorH1Opt)

      type(t_gradientRecovery), intent(INOUT) :: rgradientRecovery
      type(t_vectorScalar), intent(OUT) :: rerrorH1
      real(DP), intent(OUT), optional :: derrorH1Opt

      ! local variables
      type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
      real(DP) :: derrorH1

      ! Create new vector for local H1-error
      call lsyssc_createVector(rerrorH1, rgradientRecovery%p_rtriangulation%NEL, .true.)

      ! What error estimator should be used?
      select case(rgradientRecovery%ierrorestimator)
      case (ERREST_CSPR_FACE:ERREST_SPR_FACE)
        !------------------------------------------------------------------------
        ! Superconvergent patch recovery
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rgradientRecovery%rerrorVariable%p_rspatialDiscr
        
        ! Initialise block discretisations
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlock, 2 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlockRef, 2 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        
        ! What kind of element type is used for the FE solution?
        select case(rgradientRecovery%ieltype)
        case(1)
          ! Initialise spatial discretisations for gradient with P0-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E000, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E000, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with P1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case(2)
          ! Initialise spatial discretisations for gradient with P1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E001, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with P2-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E002, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E002, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case(11)
          ! Initialise spatial discretisations for gradient with Q0-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E010, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E010, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with Q1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E011, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E011, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case(13)
          ! Initialise spatial discretisations for gradient with Q1-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E011, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E011, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with Q2-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E013, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E013, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case(-1)
          ! Initialise spatial discretisations for gradient with P0/Q0-elements
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with P1/Q1-elements
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case DEFAULT
          call output_line('Unsupproted element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error2D')
          call sys_halt()
        end select
        
        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlock,&
                                           rgradientRecovery%rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlockRef,&
                                           rgradientRecovery%rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
                                 rgradientRecovery%rgradient)
        
        ! How should the smoothed gradient be recovered?
        select case(rgradientRecovery%ierrorestimator)
        case (ERREST_L2PROJECTION)
          call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
                                   rgradientRecovery%rgradientRef)
          
        case (ERREST_SPR_VERTEX)
          call ppgrd_calcGradSuperPatchRecov (rgradientRecovery%rerrorVariable,&
                                              rgradientRecovery%rgradientRef, PPGRD_NODEPATCH)
          
        case (ERREST_SPR_ELEMENT)
          call ppgrd_calcGradSuperPatchRecov (rgradientRecovery%rerrorVariable,&
                                              rgradientRecovery%rgradientRef, PPGRD_ELEMPATCH)
          
        case (ERREST_SPR_FACE)
          call ppgrd_calcGradSuperPatchRecov (rgradientRecovery%rerrorVariable,&
                                              rgradientRecovery%rgradientRef, PPGRD_FACEPATCH)
        end select
        
        ! Compute estimated gradient error
        call pperr_blockErrorEstimate(rgradientRecovery%rgradient,&
                                      rgradientRecovery%rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
        
      case (ERREST_LIMAVR)
        !------------------------------------------------------------------------
        ! Limited gradient averaging
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rgradientRecovery%rerrorVariable%p_rspatialDiscr
        
        ! Initialise block discretisations
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlock, 2 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        call spdiscr_initBlockDiscr2D (rgradientRecovery%rdiscrBlockRef, 2 ,&
                                       rgradientRecovery%p_rtriangulation,&
                                       rgradientRecovery%p_rboundary)
        
        ! What kind of element type is used for the FE solution?
        select case(rgradientRecovery%ieltype)
        case(1)
          ! Initialise spatial discretisations for gradient with P0-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E000, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E000, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with P1~-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E020, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E020, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case(11)
          ! Initialise spatial discretisations for gradient with Q0-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E010, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E010, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with Q1~-elements
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E030, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
              EL_E030, SPDISC_CUB_AUTOMATIC,&
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case(-1)
          ! Initialise spatial discretisations for gradient with P0/Q0-elements
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
          
          ! Initialise spatial discretisations for reference gradient with P1~/Q1~-elements
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E020, EL_E030, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
              EL_E020, EL_E030, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
          
        case DEFAULT
          call output_line('Unsupported element type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error2D')
          call sys_halt()
        end select
        
        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlock,&
                                           rgradientRecovery%rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlockRef,&
                                           rgradientRecovery%rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
                                 rgradientRecovery%rgradient)
        
        ! Recover limited averaged gradient
        call ppgrd_calcGradLimAvgP1Q1cnf (rgradientRecovery%rerrorVariable,&
                                          rgradientRecovery%rgradientRef)
        
        ! Compute gradient error
        call pperr_blockErrorEstimate(rgradientRecovery%rgradient,&
                                      rgradientRecovery%rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
      case DEFAULT
        call output_line('Invalid type of error estimator!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error2D')
        call sys_halt()
      end select
      
      ! Return global gradient error?
      if (present(derrorH1Opt)) derrorH1Opt = derrorH1
    end subroutine calcH1Error2D
  end subroutine errest_calcH1Error

  !*****************************************************************************

!<subroutine>

  subroutine errest_calcFirstDiffIndicator(rdifferenceIndicator, rindicator)

!<description>
    ! This subroutine computes the first-difference indicator by Loehner.
!</description>

!<inputoutput>
    ! difference indicator
    type(t_differenceIndicator), intent(INOUT) :: rdifferenceIndicator
!</inputoutput>

!<output>
    ! second-difference indicator
    type(t_vectorScalar), intent(OUT) :: rindicator
!</output>
!</subroutine>
    
    ! Scalar vector to compute median
    type(t_vectorScalar) :: rvector

    ! Pointer to data vector
    real(DP), dimension(:), pointer :: p_Ddata

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to error variable
    real(DP), dimension(:), pointer :: p_DerrorVariable

    ! Pointer to element indicator
    real(DP), dimension(:), pointer :: p_Dindicator

    ! Pointer to vertices at element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! local variables
    integer :: iel,i1,i2,i3,i4
    real(DP) :: u1,u2,u3,u4
    real(DP) :: stdDev,meanDev,median
    
    ! Set pointer to triangulation
    p_rtriangulation => rdifferenceIndicator%p_rtriangulation

    ! Create new nodal vectors
    call lsyssc_createVector(rindicator,   p_rtriangulation%NEL, .true.)
    call lsyssc_createVector(rvector,      p_rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lsyssc_getbase_double(rvector,    p_Ddata)

    ! Set pointers
    call lsyssc_getbase_double(rdifferenceIndicator%rerrorVariable,&
                               p_DerrorVariable)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)

    ! Loop over all elements in triangulation
    do iel = 1, p_rtriangulation%NEL

      ! What type of element are we
      select case(tria_getNVE(p_rtriangulation, iel))
        
      case(TRIA_NVETRI2D)

        ! Determine global degrees of freedom
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Determine solution at global degrees of freedom
        u1 = p_DerrorVariable(i1)
        u2 = p_DerrorVariable(i2)
        u3 = p_DerrorVariable(i3)

        ! Compute absolute value of first difference
        p_Dindicator(iel) = max( abs(u2-u1), abs(u3-u2), abs(u1-u3))

      case (TRIA_NVEQUAD2D)

        ! Determine global degrees of freedom
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        i4 = p_IverticesAtElement(4, iel)
        
        ! Determine solution at global degrees of freedom
        u1 = p_DerrorVariable(i1)
        u2 = p_DerrorVariable(i2)
        u3 = p_DerrorVariable(i3)
        u4 = p_DerrorVariable(i4)

        ! Compute absolute value of first difference
        p_Dindicator(iel) = max( abs(u2-u1), abs(u3-u2), abs(u4-u3), abs(u1-u4))

      case DEFAULT
        call output_line('Unsupproted element type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'errest_calcFirstDifferenceIndicator')
        call sys_halt()
      end select
    end do

    ! Compute standard and mean deviation
    stdDev  = mprim_stdDeviationDble(p_Dindicator)
    meanDev = mprim_meanDeviationDble(p_Dindicator)

    ! Loop over all element in triangulation
    do iel = 1, p_rtriangulation%NEL

      ! Normalize first difference
      p_Dindicator(iel) = (p_Dindicator(iel)-meanDev)/&
          max(stdDev, rdifferenceIndicator%dnoiseFilter*meanDev,&
              rdifferenceIndicator%dabsFilter)
    end do

    ! Compute median
    call lsyssc_copyVector(rindicator, rvector)
    call sort_dp(p_Ddata, SORT_QUICK)

    if (mod(p_rtriangulation%NEL,2) .eq. 0) then
      median = (p_Ddata(p_rtriangulation%NEL/2)+&
                p_Ddata(1+p_rtriangulation%NEL/2))/2._DP
    else
      median = p_Ddata(1+p_rtriangulation%NEL/2)
    end if

    ! Subtract median from grid indicator
    p_Dindicator = p_Dindicator-median

    ! Release temporal vector
    call lsyssc_releaseVector(rvector)
  end subroutine errest_calcFirstDiffIndicator

  !*****************************************************************************

!<subroutine>

  subroutine errest_calcSecondDiffIndicator(rdifferenceIndicator, rindicator)

!<description>
    ! This subroutine computes the second-difference indicator by Loehner.
!</description>

!<inputoutput>
    ! difference indicator
    type(t_differenceIndicator), intent(INOUT) :: rdifferenceIndicator
!</inputoutput>

!<output>
    ! second-difference indicator
    type(t_vectorScalar), intent(OUT) :: rindicator
!</output>
!</subroutine>

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to element indicator
    real(DP), dimension(:), pointer :: p_Dindicator
    
    ! Pointer to scalar error vector
    real(DP), dimension(:), pointer :: p_DerrorVariable

    ! Pointer to vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    ! Pointer to vertices at element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to neighbours at element
    integer, dimension(:,:), pointer :: p_IneighboursAtElement

    ! Pointer to nodal property
    integer(I32), dimension(:), pointer :: p_InodalProperty

    ! Handle for nodal coefficients
    integer :: h_Dcoefficients
    
    ! Pointer to p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_Dcoefficients

    ! Handle for nodal indicator
    integer :: h_DnodalIndicator
    
    ! Pointer to p_DnodalIndicator
    real(DP), dimension(:), pointer :: p_DnodalIndicator

    ! local variable
    integer(I32), dimension(2) :: Isize

    ! Set pointer to triangulation
    p_rtriangulation => rdifferenceIndicator%p_rtriangulation

    ! Check if triangulation exists
    if (.not.associated(p_rtriangulation)) then
      call output_line('No triangulation structure available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'errest_calcSecondDiffIndicator')
      call sys_halt()
    end if

    ! Create scalar vector for element indicator
    call lsyssc_createVector(rindicator, p_rtriangulation%NEL, .true.)

    ! Create temporal memory
    h_Dcoefficients = ST_NOHANDLE
    Isize = (/12, p_rtriangulation%NVT/)
    call storage_new('errest_calcSecondDiffIndicator',' Dcoefficients',&
        Isize, ST_DOUBLE, h_Dcoefficients, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double2d(h_Dcoefficients, p_Dcoefficients)

    ! Create temporal memory
    h_DnodalIndicator = ST_NOHANDLE
    call storage_new('errest_calcSecondDiffIndicator',' DnodalIndicator',&
        p_rtriangulation%NVT, ST_DOUBLE, h_DnodalIndicator, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double(h_DnodalIndicator, p_DnodalIndicator)
    
    ! Set pointers
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call storage_getbase_int(p_rtriangulation%h_InodalProperty, p_InodalProperty)
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call lsyssc_getbase_double(rdifferenceIndicator%rerrorVariable, p_DerrorVariable)
    call lsyssc_getbase_double(rindicator, p_Dindicator)

    ! How many spatial dimensions are we?
    select case (p_rtriangulation%ndim)
    case (NDIM2D)
      call doSecondDiffIndicator2D(p_DerrorVariable, p_DvertexCoords, &
          p_IverticesAtElement, p_IneighboursAtElement, p_InodalProperty,&
          rdifferenceIndicator%dnoiseFilter, rdifferenceIndicator%dabsFilter,&
          p_rtriangulation%NEL, p_rtriangulation%NVT,&
          p_Dcoefficients, p_DnodalIndicator, p_Dindicator)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                        OU_CLASS_ERROR,OU_MODE_STD,'errest_calcSecondDiffIndicator')
      call sys_halt()
    end select

    ! Release temporal memory
    call storage_free(h_Dcoefficients)
    call storage_free(h_DnodalIndicator)

  contains
    
    ! Here, the real working routines follow.

    !**************************************************************
    ! Second-difference indicator in 2D

    subroutine doSecondDiffIndicator2D(DerrorVariable, DvertexCoords,&
        IverticesAtElement, IneighboursAtElement, InodalProperty, &
        dweight, dfilter, NEL, NVT, Dcoefficients, DnodalIndicator, Dindicator)
      
      real(DP), dimension(:), intent(IN) :: DerrorVariable
      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      integer, dimension(:,:), intent(IN) :: IverticesAtElement
      integer, dimension(:,:), intent(IN) :: IneighboursAtElement
      integer(I32), dimension(:), intent(IN) :: InodalProperty
      real(DP), intent(IN) :: dweight,dfilter
      integer, intent(IN) :: NEL
      integer, intent(IN) :: NVT

      real(DP), dimension(:,:), intent(OUT) :: Dcoefficients
      real(DP), dimension(:), intent(OUT) :: DnodalIndicator
      real(DP), dimension(:), intent(OUT) :: Dindicator

      ! local variables
      real(DP), dimension(2) :: DabsDeriv,Dderiv,Dfunc
      real(DP), dimension(2,3) :: Dbas,dabsBas
      real(DP) :: darea,ddet,dlength,dbdrInt
      real(DP) :: daux1,daux2,daux3
      real(DP) :: u1,u2,u3,u4,x1,x2,x3,x4,y1,y2,y3,y4
      integer :: iel,ivt,i1,i2,i3,i4,ive,ncontributions



      ! Clear array for coefficients
      call lalg_clearVectorDble2D(Dcoefficients)

      ! Loop over all elements in triangulation
      do iel = 1, NEL
        
        ! What type of element are we
        select case(tria_getNVE(IverticesAtElement, iel))
          
        case(TRIA_NVETRI2D)
          
          ! Determine global degrees of freedom
          i1 = IverticesAtElement(1, iel)
          i2 = IverticesAtElement(2, iel)
          i3 = IverticesAtElement(3, iel)

          ! Determine vertex coordinates
          x1 = DvertexCoords(1, i1)
          y1 = DvertexCoords(2, i1)
          x2 = DvertexCoords(1, i2)
          y2 = DvertexCoords(2, i2)
          x3 = DvertexCoords(1, i3)
          y3 = DvertexCoords(2, i3)

          ! Determine solution values at vertices
          u1 = DerrorVariable(i1)
          u2 = DerrorVariable(i2)
          u3 = DerrorVariable(i3)


          ! Calculate determinant and area for triangle
          ddet  = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
          darea = 0.5_DP*abs(ddet)

          ! Calculate averages of solution on triangle
          Dfunc(1)  = (abs(u1*(y2-y3)) + abs(u2*(y3-y1)) + abs(u3*(y1-y2))) / abs(ddet)
          Dfunc(2)  = (abs(u1*(x3-x2)) + abs(u2*(x1-x3)) + abs(u3*(x2-x1))) / abs(ddet)
          
          ! Calculate derivatives of solution on triangle
          Dderiv(1) = (u1*(y2-y3) + u2*(y3-y1) + u3*(y1-y2)) / ddet
          Dderiv(2) = (u1*(x3-x2) + u2*(x1-x3) + u3*(x2-x1)) / ddet

          ! Calculate absolute values of derivatives
          DabsDeriv = abs(Dderiv)

          ! Calculate derivatives of basis functions on triangle
          Dbas(1,1) = (y2-y3) / ddet
          Dbas(2,1) = (x3-x2) / ddet
          Dbas(1,2) = (y3-y1) / ddet
          Dbas(2,2) = (x1-x3) / ddet
          Dbas(1,3) = (y1-y2) / ddet
          Dbas(2,3) = (x2-x1) / ddet

          ! Calculate absolute values of basis functions
          DabsBas = abs(Dbas)
          
          ! Update nodal coefficient vector for node I1
          Dcoefficients(1,i1) = Dcoefficients(1,i1) + Dbas(1,1) * Dderiv(1) * darea
          Dcoefficients(2,i1) = Dcoefficients(2,i1) + Dbas(1,1) * Dderiv(2) * darea
          Dcoefficients(3,i1) = Dcoefficients(3,i1) + Dbas(2,1) * Dderiv(1) * darea
          Dcoefficients(4,i1) = Dcoefficients(4,i1) + Dbas(2,1) * Dderiv(2) * darea

          Dcoefficients(5,i1) = Dcoefficients(5,i1) + DabsBas(1,1) * DabsDeriv(1) * darea
          Dcoefficients(6,i1) = Dcoefficients(6,i1) + DabsBas(1,1) * DabsDeriv(2) * darea
          Dcoefficients(7,i1) = Dcoefficients(7,i1) + DabsBas(2,1) * DabsDeriv(1) * darea
          Dcoefficients(8,i1) = Dcoefficients(8,i1) + DabsBas(2,1) * DabsDeriv(2) * darea

          Dcoefficients(9, i1) = Dcoefficients(9, i1) + DabsBas(1,1) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i1) = Dcoefficients(10,i1) + DabsBas(1,1) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i1) = Dcoefficients(11,i1) + DabsBas(2,1) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i1) = Dcoefficients(12,i1) + DabsBas(2,1) * dweight * Dfunc(2) * darea

          
          ! Update nodal coefficient vector for node I2
          Dcoefficients(1,i2) = Dcoefficients(1,i2) + Dbas(1,2) * Dderiv(1) * darea
          Dcoefficients(2,i2) = Dcoefficients(2,i2) + Dbas(1,2) * Dderiv(2) * darea
          Dcoefficients(3,i2) = Dcoefficients(3,i2) + Dbas(2,2) * Dderiv(1) * darea
          Dcoefficients(4,i2) = Dcoefficients(4,i2) + Dbas(2,2) * Dderiv(2) * darea

          Dcoefficients(5,i2) = Dcoefficients(5,i2) + DabsBas(1,2) * DabsDeriv(1) * darea
          Dcoefficients(6,i2) = Dcoefficients(6,i2) + DabsBas(1,2) * DabsDeriv(2) * darea
          Dcoefficients(7,i2) = Dcoefficients(7,i2) + DabsBas(2,2) * DabsDeriv(1) * darea
          Dcoefficients(8,i2) = Dcoefficients(8,i2) + DabsBas(2,2) * DabsDeriv(2) * darea

          Dcoefficients(9, i2) = Dcoefficients(9, i2) + DabsBas(1,2) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i2) = Dcoefficients(10,i2) + DabsBas(1,2) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i2) = Dcoefficients(11,i2) + DabsBas(2,2) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i2) = Dcoefficients(12,i2) + DabsBas(2,2) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I3
          Dcoefficients(1,i3) = Dcoefficients(1,i3) + Dbas(1,3) * Dderiv(1) * darea
          Dcoefficients(2,i3) = Dcoefficients(2,i3) + Dbas(1,3) * Dderiv(2) * darea
          Dcoefficients(3,i3) = Dcoefficients(3,i3) + Dbas(2,3) * Dderiv(1) * darea
          Dcoefficients(4,i3) = Dcoefficients(4,i3) + Dbas(2,3) * Dderiv(2) * darea

          Dcoefficients(5,i3) = Dcoefficients(5,i3) + DabsBas(1,3) * DabsDeriv(1) * darea
          Dcoefficients(6,i3) = Dcoefficients(6,i3) + DabsBas(1,3) * DabsDeriv(2) * darea
          Dcoefficients(7,i3) = Dcoefficients(7,i3) + DabsBas(2,3) * DabsDeriv(1) * darea
          Dcoefficients(8,i3) = Dcoefficients(8,i3) + DabsBas(2,3) * DabsDeriv(2) * darea

          Dcoefficients(9, i3) = Dcoefficients(9, i3) + DabsBas(1,3) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i3) = Dcoefficients(10,i3) + DabsBas(1,3) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i3) = Dcoefficients(11,i3) + DabsBas(2,3) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i3) = Dcoefficients(12,i3) + DabsBas(2,3) * dweight * Dfunc(2) * darea


        case (TRIA_NVEQUAD2D)
          
          ! Determine global degrees of freedom
          i1 = IverticesAtElement(1, iel)
          i2 = IverticesAtElement(2, iel)
          i3 = IverticesAtElement(3, iel)
          i4 = IverticesAtElement(4, iel)
          
          ! Determine vertex coordinates
          x1 = DvertexCoords(1, i1)
          y1 = DvertexCoords(2, i1)
          x2 = DvertexCoords(1, i2)
          y2 = DvertexCoords(2, i2)
          x3 = DvertexCoords(1, i3)
          y3 = DvertexCoords(2, i3)
          x4 = DvertexCoords(1, i4)
          y4 = DvertexCoords(2, i4)

          ! Determine solution values at vertices
          u1 = DerrorVariable(i1)
          u2 = DerrorVariable(i2)
          u3 = DerrorVariable(i3)
          u4 = DerrorVariable(i4)


          ! Calculate determinant and area for triangle (1-2-3)
          ddet  = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
          darea = 0.5_DP*abs(ddet)

          ! Calculate averages of solution on triangle
          Dfunc(1)  = (abs(u1*(y2-y3)) + abs(u2*(y3-y1)) + abs(u3*(y1-y2))) / abs(ddet)
          Dfunc(2)  = (abs(u1*(x3-x2)) + abs(u2*(x1-x3)) + abs(u3*(x2-x1))) / abs(ddet)

          ! Calculate derivatives of solution on triangle
          Dderiv(1) = (u1*(y2-y3) + u2*(y3-y1) + u3*(y1-y2)) / ddet
          Dderiv(2) = (u1*(x3-x2) + u2*(x1-x3) + u3*(x2-x1)) / ddet
          
          ! Calculate absolute values of derivatives
          DabsDeriv = abs(Dderiv)

          ! Calculate derivatives of basis functions on triangle
          Dbas(1,1) = (y2-y3) / ddet
          Dbas(2,1) = (x3-x2) / ddet
          Dbas(1,2) = (y3-y1) / ddet
          Dbas(2,2) = (x1-x3) / ddet
          Dbas(1,3) = (y1-y2) / ddet
          Dbas(2,3) = (x2-x1) / ddet

          ! Calculate absolute values of basis functions
          DabsBas = abs(Dbas)
          

          ! Update nodal coefficient vector for node I1
          Dcoefficients(1,i1) = Dcoefficients(1,i1) + Dbas(1,1) * Dderiv(1) * darea
          Dcoefficients(2,i1) = Dcoefficients(2,i1) + Dbas(1,1) * Dderiv(2) * darea
          Dcoefficients(3,i1) = Dcoefficients(3,i1) + Dbas(2,1) * Dderiv(1) * darea
          Dcoefficients(4,i1) = Dcoefficients(4,i1) + Dbas(2,1) * Dderiv(2) * darea

          Dcoefficients(5,i1) = Dcoefficients(5,i1) + DabsBas(1,1) * DabsDeriv(1) * darea
          Dcoefficients(6,i1) = Dcoefficients(6,i1) + DabsBas(1,1) * DabsDeriv(2) * darea
          Dcoefficients(7,i1) = Dcoefficients(7,i1) + DabsBas(2,1) * DabsDeriv(1) * darea
          Dcoefficients(8,i1) = Dcoefficients(8,i1) + DabsBas(2,1) * DabsDeriv(2) * darea

          Dcoefficients(9, i1) = Dcoefficients(9, i1) + DabsBas(1,1) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i1) = Dcoefficients(10,i1) + DabsBas(1,1) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i1) = Dcoefficients(11,i1) + DabsBas(2,1) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i1) = Dcoefficients(12,i1) + DabsBas(2,1) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I2
          Dcoefficients(1,i2) = Dcoefficients(1,i2) + Dbas(1,2) * Dderiv(1) * darea
          Dcoefficients(2,i2) = Dcoefficients(2,i2) + Dbas(1,2) * Dderiv(2) * darea
          Dcoefficients(3,i2) = Dcoefficients(3,i2) + Dbas(2,2) * Dderiv(1) * darea
          Dcoefficients(4,i2) = Dcoefficients(4,i2) + Dbas(2,2) * Dderiv(2) * darea

          Dcoefficients(5,i2) = Dcoefficients(5,i2) + DabsBas(1,2) * DabsDeriv(1) * darea
          Dcoefficients(6,i2) = Dcoefficients(6,i2) + DabsBas(1,2) * DabsDeriv(2) * darea
          Dcoefficients(7,i2) = Dcoefficients(7,i2) + DabsBas(2,2) * DabsDeriv(1) * darea
          Dcoefficients(8,i2) = Dcoefficients(8,i2) + DabsBas(2,2) * DabsDeriv(2) * darea

          Dcoefficients(9, i2) = Dcoefficients(9, i2) + DabsBas(1,2) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i2) = Dcoefficients(10,i2) + DabsBas(1,2) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i2) = Dcoefficients(11,i2) + DabsBas(2,2) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i2) = Dcoefficients(12,i2) + DabsBas(2,2) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I3
          Dcoefficients(1,i3) = Dcoefficients(1,i3) + Dbas(1,3) * Dderiv(1) * darea
          Dcoefficients(2,i3) = Dcoefficients(2,i3) + Dbas(1,3) * Dderiv(2) * darea
          Dcoefficients(3,i3) = Dcoefficients(3,i3) + Dbas(2,3) * Dderiv(1) * darea
          Dcoefficients(4,i3) = Dcoefficients(4,i3) + Dbas(2,3) * Dderiv(2) * darea
          
          Dcoefficients(5,i3) = Dcoefficients(5,i3) + DabsBas(1,3) * DabsDeriv(1) * darea
          Dcoefficients(6,i3) = Dcoefficients(6,i3) + DabsBas(1,3) * DabsDeriv(2) * darea
          Dcoefficients(7,i3) = Dcoefficients(7,i3) + DabsBas(2,3) * DabsDeriv(1) * darea
          Dcoefficients(8,i3) = Dcoefficients(8,i3) + DabsBas(2,3) * DabsDeriv(2) * darea

          Dcoefficients(9, i3) = Dcoefficients(9, i3) + DabsBas(1,3) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i3) = Dcoefficients(10,i3) + DabsBas(1,3) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i3) = Dcoefficients(11,i3) + DabsBas(2,3) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i3) = Dcoefficients(12,i3) + DabsBas(2,3) * dweight * Dfunc(2) * darea

          
          ! Calculate determinant and area for triangle (1-3-4)
          ddet  = (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)
          darea = 0.5_DP*abs(ddet)

          ! Calculate averages of solution on triangle
          Dfunc(1)  = (abs(u1*(y3-y4)) + abs(u3*(y4-y1)) + abs(u4*(y1-y3))) / abs(ddet)
          Dfunc(2)  = (abs(u1*(x4-x3)) + abs(u3*(x1-x4)) + abs(u4*(x3-x1))) / abs(ddet)

          ! Calculate derivatives of solution on triangle
          Dderiv(1) = (u1*(y3-y4) + u3*(y4-y1) + u4*(y1-y3)) / ddet
          Dderiv(2) = (u1*(x4-x3) + u3*(x1-x4) + u4*(x3-x1)) / ddet

          ! Calculate absolute values of derivatives
          DabsDeriv = abs(Dderiv)

          ! Calculate derivatives of basis functions on triangle
          Dbas(1,1) = (y3-y4) / ddet
          Dbas(2,1) = (x4-x3) / ddet
          Dbas(1,2) = (y4-y1) / ddet
          Dbas(2,2) = (x1-x4) / ddet
          Dbas(1,3) = (y1-y3) / ddet
          Dbas(2,3) = (x3-x1) / ddet
          
          ! Calculate absolute values of basis functions
          DabsBas = abs(Dbas)

          ! Update nodal coefficient vector for node I1
          Dcoefficients(1,i1) = Dcoefficients(1,i1) + Dbas(1,1) * Dderiv(1) * darea
          Dcoefficients(2,i1) = Dcoefficients(2,i1) + Dbas(1,1) * Dderiv(2) * darea
          Dcoefficients(3,i1) = Dcoefficients(3,i1) + Dbas(2,1) * Dderiv(1) * darea
          Dcoefficients(4,i1) = Dcoefficients(4,i1) + Dbas(2,1) * Dderiv(2) * darea

          Dcoefficients(5,i1) = Dcoefficients(5,i1) + DabsBas(1,1) * DabsDeriv(1) * darea
          Dcoefficients(6,i1) = Dcoefficients(6,i1) + DabsBas(1,1) * DabsDeriv(2) * darea
          Dcoefficients(7,i1) = Dcoefficients(7,i1) + DabsBas(2,1) * DabsDeriv(1) * darea
          Dcoefficients(8,i1) = Dcoefficients(8,i1) + DabsBas(2,1) * DabsDeriv(2) * darea

          Dcoefficients(9, i1) = Dcoefficients(9, i1) + DabsBas(1,1) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i1) = Dcoefficients(10,i1) + DabsBas(1,1) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i1) = Dcoefficients(11,i1) + DabsBas(2,1) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i1) = Dcoefficients(12,i1) + DabsBas(2,1) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I3
          Dcoefficients(1,i3) = Dcoefficients(1,i3) + Dbas(1,2) * Dderiv(1) * darea
          Dcoefficients(2,i3) = Dcoefficients(2,i3) + Dbas(1,2) * Dderiv(2) * darea
          Dcoefficients(3,i3) = Dcoefficients(3,i3) + Dbas(2,2) * Dderiv(1) * darea
          Dcoefficients(4,i3) = Dcoefficients(4,i3) + Dbas(2,2) * Dderiv(2) * darea

          Dcoefficients(5,i3) = Dcoefficients(5,i3) + DabsBas(1,2) * DabsDeriv(1) * darea
          Dcoefficients(6,i3) = Dcoefficients(6,i3) + DabsBas(1,2) * DabsDeriv(2) * darea
          Dcoefficients(7,i3) = Dcoefficients(7,i3) + DabsBas(2,2) * DabsDeriv(1) * darea
          Dcoefficients(8,i3) = Dcoefficients(8,i3) + DabsBas(2,2) * DabsDeriv(2) * darea

          Dcoefficients(9, i3) = Dcoefficients(9, i3) + DabsBas(1,2) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i3) = Dcoefficients(10,i3) + DabsBas(1,2) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i3) = Dcoefficients(11,i3) + DabsBas(2,2) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i3) = Dcoefficients(12,i3) + DabsBas(2,2) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I4
          Dcoefficients(1,i4) = Dcoefficients(1,i4) + Dbas(1,3) * Dderiv(1) * darea
          Dcoefficients(2,i4) = Dcoefficients(2,i4) + Dbas(1,3) * Dderiv(2) * darea
          Dcoefficients(3,i4) = Dcoefficients(3,i4) + Dbas(2,3) * Dderiv(1) * darea
          Dcoefficients(4,i4) = Dcoefficients(4,i4) + Dbas(2,3) * Dderiv(2) * darea

          Dcoefficients(5,i4) = Dcoefficients(5,i4) + DabsBas(1,3) * DabsDeriv(1) * darea
          Dcoefficients(6,i4) = Dcoefficients(6,i4) + DabsBas(1,3) * DabsDeriv(2) * darea
          Dcoefficients(7,i4) = Dcoefficients(7,i4) + DabsBas(2,3) * DabsDeriv(1) * darea
          Dcoefficients(8,i4) = Dcoefficients(8,i4) + DabsBas(2,3) * DabsDeriv(2) * darea

          Dcoefficients(9, i4) = Dcoefficients(9, i4) + DabsBas(1,3) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i4) = Dcoefficients(10,i4) + DabsBas(1,3) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i4) = Dcoefficients(11,i4) + DabsBas(2,3) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i4) = Dcoefficients(12,i4) + DabsBas(2,3) * dweight * Dfunc(2) * darea


        case DEFAULT
          call output_line('Unsupproted element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'doSecondDiffIndicator2D')
          call sys_halt()
        end select
      end do

      ! Loop over all vertices in triangulation
      do ivt = 1, NVT
        ! Compute numerator
        daux1 = Dcoefficients(1, ivt)*Dcoefficients(1, ivt) +&
                Dcoefficients(2, ivt)*Dcoefficients(2, ivt) +&
                Dcoefficients(3, ivt)*Dcoefficients(3, ivt) +&
                Dcoefficients(4, ivt)*Dcoefficients(4, ivt)

        ! Compute derivative part of denominator
        daux2 = Dcoefficients(5, ivt)*Dcoefficients(5, ivt) +&
                Dcoefficients(6, ivt)*Dcoefficients(6, ivt) +&
                Dcoefficients(7, ivt)*Dcoefficients(7, ivt) +&
                Dcoefficients(8, ivt)*Dcoefficients(8, ivt)

        ! Compute average part of denominator
        daux3 = Dcoefficients( 9, ivt)*Dcoefficients( 9, ivt) +&
                Dcoefficients(10, ivt)*Dcoefficients(10, ivt) +&
                Dcoefficients(11, ivt)*Dcoefficients(11, ivt) +&
                Dcoefficients(12, ivt)*Dcoefficients(12, ivt)

        ! Compute nodal indicator
        DnodalIndicator(ivt) = sqrt(daux1/(daux2+max(dfilter, daux3)))
      end do
      
      
      ! Clear array
      call lalg_clearVectorDble(Dindicator)

      ! Loop over all elements in triangulation
      do iel = 1, NEL
        
        ! Initialize number of contributions
        ncontributions = 0

        ! Loop over all vertices of the element
        do ive = 1, tria_getNVE(IverticesAtElement, iel)

          ! Determine global degree of freedom
          ivt = IverticesAtElement(ive, iel)

          ! Skip boundary nodes
          if (InodalProperty(ivt) .ne. 0) cycle

          ! Apply nodal contribution to element indicator
          Dindicator(iel) = Dindicator(iel) + DnodalIndicator(ivt)

          ! Increase number of contributions
          ncontributions = ncontributions+1
        end do
        
        ! Scale element indicator by number of contributions
        if (ncontributions .ne. 0) then
          Dindicator(iel) = Dindicator(iel) / real(ncontributions, DP)
        else
          Dindicator(iel) = 0._DP
        end if
      end do
      
    end subroutine doSecondDiffIndicator2D
  end subroutine errest_calcSecondDiffIndicator
  
  !*****************************************************************************

!<subroutine>

  subroutine errest_calcProtectionLayers(rerrorEstimator, rindicator, dthreshold)

!<description>
    ! This subroutine adjusts the grid indicator to include a prescribed
    ! number of protection layers based on the given threshold level
!</description>

!<input>
    ! error estimator
    type(t_errorEstimator), intent(IN) :: rerrorEstimator

    ! threshold value
    real(DP), intent(IN) :: dthreshold
!</input>

!<inputoutput>
    ! elementwise grid indicator
    type(t_vectorScalar), intent(INOUT) :: rindicator
!</inputoutput>
!</subroutine>

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to element indicator
    real(DP), dimension(:), pointer :: p_Dindicator

    ! Pointer to vertices at element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to neighbours at element
    integer, dimension(:,:), pointer :: p_IneighboursAtElement

    ! Pointer to BisactiveElement
    logical, dimension(:), pointer :: p_BisactiveElement

    ! Handle for h_Bisactiveelement
    integer :: h_BisactiveElement
    

    ! local variables
    integer :: iprotectlayer


    ! What kind of error estimator are we?
    select case(rerrorEstimator%ierrorestimator)

    case (ERREST_CSPR_FACE: ERREST_LIMAVR)
      ! Set pointer to triangulation
      p_rtriangulation => rerrorEstimator%p_rgradientRecovery%p_rtriangulation

    case (ERREST_FIRSTDIFF:ERREST_SECONDDIFF)
      ! Set pointer to triangulation
      p_rtriangulation => rerrorEstimator%p_rdifferenceIndicator%p_rtriangulation

    case DEFAULT
      call output_line('Invalid type of error estimator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'errest_calcProtectionLayers')
      call sys_halt()
    end select

    ! Check if triangulation exists
    if (.not.associated(p_rtriangulation)) then
      call output_line('No triangulation structure available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'errest_calcProtectionLayers')
      call sys_halt()
    end if

    ! Create memory
    h_BisactiveElement = ST_NOHANDLE
    call storage_new('errest_calcProtectionLayers',' BisactiveElement',&
        p_rtriangulation%NEL, ST_LOGICAL, h_BisactiveElement, ST_NEWBLOCK_NOINIT)
    call storage_getbase_logical(h_BisactiveElement, p_BisactiveElement)

    ! Set pointers
    call storage_getbase_int2D(p_rtriangulation%h_IneighboursAtElement,&
                               p_IneighboursAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
    call lsyssc_getbase_double(rindicator, p_Dindicator)

    ! Compute protection layers
    do iprotectlayer = 1, rerrorEstimator%nprotectlayers

      p_BisActiveElement = .false.     
      call doProtectionLayer(p_IverticesAtElement, p_IneighboursAtElement,&
          p_rtriangulation%NEL, dthreshold, p_Dindicator, p_BisActiveElement)
    end do

    ! Release memory
    call storage_free(h_BisactiveElement)

  contains
    
    ! Here, the real working routines follow.

    !**************************************************************
    ! Compute one protection layer

    subroutine doProtectionLayer(IverticesAtElement, IneighboursAtElement, NEL,&
                                 dthreshold, Dindicator, BisactiveElement)

      integer, dimension(:,:), intent(IN) :: IverticesAtElement
      integer, dimension(:,:), intent(IN) :: IneighboursAtElement     
      real(DP), intent(IN) :: dthreshold
      integer, intent(IN) :: NEL
      real(DP), dimension(:), intent(INOUT) :: Dindicator
      logical, dimension(:), intent(INOUT) :: BisactiveElement
      

      ! local variables
      integer :: iel,jel,ive

      ! Loop over all elements in triangulation
      do iel = 1, NEL
        
        ! Do nothing if element belongs to active layer
        if (BisactiveElement(iel)) cycle

        ! Do nothing if element indicator does not exceed threshold
        if (Dindicator(iel) .lt. dthreshold) cycle

        ! Loop over neighbouring elements
        do ive = 1, tria_getNVE(IverticesAtElement, iel)
          
          ! Get number of neighbouring element
          jel = IneighboursAtElement(ive, iel)

          ! Do nothing at the boundary
          if (jel .eq. 0) cycle

          ! Check if element belongs to active layer
          if (BisactiveElement(jel)) then
            ! If yes, then just update the element indicator
            Dindicator(jel) = max(Dindicator(jel), Dindicator(iel))
          else
            ! Otherwise, we have to check if the neighbouring element
            ! exceeds the prescribed threshold level. If this is the case
            ! it will be processed later or has already been processed
            if (Dindicator(jel) .lt. dthreshold) then
              Dindicator(jel) = max(Dindicator(jel), Dindicator(iel))
              BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do
    end subroutine doProtectionLayer
  end subroutine errest_calcProtectionLayers

  !*****************************************************************************

!<subroutine>

  subroutine errest_calcGridIndicator(rerrorEstimator, rindicator, dglobalError)
    
!<description>
    ! This subroutine computes the elementwise grid indicator
    ! based on the adaptation strategy given in the error estimator
!</description>

!<input>
    ! error estimator
    type(t_errorEstimator), intent(IN) :: rerrorEstimator

    ! OPTIONAL: norm of global error
    real(DP), intent(IN), optional :: dglobalError
!</input>

!<inputoutput>
    ! elementwise grid indicator
    type(t_vectorScalar), intent(INOUT) :: rindicator
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP) :: derror,dsolution,dfactor,dmeanval,dmaxval
    integer :: i

    ! What kind of strategy should be used?
    select case(rerrorEstimator%igridindicator)

    case (ERREST_ASIS)   
      ! That's simple, do nothing.


    case (ERREST_EQUIDIST)
      ! Try to equidistribute the relative percentage error
      if (present(dglobalError)) then
        derror = dglobalError
      else
        ! Now we compute the global error
      end if
      
      ! We also need the global norm of the scalar error variable
      select case(rerrorEstimator%ierrorestimator)
        
      case (ERREST_CSPR_FACE:ERREST_LIMAVR)
        call pperr_scalar(rerrorEstimator%p_rgradientRecovery%rerrorVariable,&
                          PPERR_L2ERROR, dsolution)

      case (ERREST_FIRSTDIFF:ERREST_SECONDDIFF)
        call pperr_scalar(rerrorEstimator%p_rdifferenceIndicator%rerrorVariable,&
                          PPERR_L2ERROR, dsolution)

      case DEFAULT
        call output_line('Invalid type of error estimator!',&
            OU_CLASS_ERROR,OU_MODE_STD,'errest_calcGridIndicator')
        call sys_halt()
      end select
      
      ! Compute permissible element error
      dfactor = sqrt(dsolution*dsolution + derror*derror)/sqrt(real(rindicator%NEQ, DP))

      ! Scale element error by permissible error
      call lsyssc_scaleVector(rindicator, 1/dfactor)


    case (ERREST_LOGEQUIDIST)

      ! Set pointer
      call lsyssc_getbase_double(rindicator, p_Ddata)
      
      ! Determine largest error value
      dmaxval = - SYS_MAXREAL
      do i = 1, size(p_Ddata)
        dmaxval = max(dmaxval, p_Ddata(i))
      end do

      ! Normalize error by largest value
      do i = 1, size(p_Ddata)
        p_Ddata(i) = p_Ddata(i)/dmaxval
      end do
      
      ! Initialize mean value
      dmeanval = 0.0_DP

      ! Loop over all contributions
      do i = 1, size(p_Ddata)
        p_Ddata(i) = log(max(exp(-20.0_DP), p_Ddata(i)))
        dmeanval   = dmeanval + p_Ddata(i)
      end do
      
      ! Calculate mean
      dmeanval = dmeanval/real(size(p_Ddata), DP)
      
      ! Subtract mean value from grid indicator
      do i = 1, size(p_Ddata)
        p_Ddata(i) = p_Ddata(i)-dmeanval
      end do


    case (ERREST_AUTORMS)

      ! Set pointer
      call lsyssc_getbase_double(rindicator, p_Ddata)

      ! Initialize mean value
      dmeanval = 0.0_DP

      ! Loop over all  contributions
      do i = 1, size(p_Ddata)
        dmeanval = dmeanval + p_Ddata(i)**2
      end do

      ! Calculate root mean value
      dmeanval = sqrt(dmeanval/real(size(p_Ddata), DP))

      ! Divide grid indicator by RMS
      do i = 1, size(p_Ddata)
        p_Ddata(i) = p_Ddata(i)/dmeanval
      end do
      

    case DEFAULT
      call output_line('Invalid type of grid indicator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'errest_calcGridIndicator')
      call sys_halt()
    end select

  end subroutine errest_calcGridIndicator
end module errorestimation

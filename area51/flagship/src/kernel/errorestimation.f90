!##############################################################################
!# ****************************************************************************
!# <name> errorestimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures and subroutines
!# for error estimation.
!#
!# 1.) errest_initFromParameterlist
!#     -> Initializes the error estimator from parameter file
!#
!# 2.) errest_releaseErrorEstimator
!#     -> Releases the error estimator and all of its structures
!#
!# 3.) errest_addVariable
!#     -> Adds a scalar variable to the error estimator
!#
!# 4.) errest_clearErrorEstimator
!#     -> Clear all temporal vectors from the error estimator
!#
!# 5.) errest_calcH1Error
!#     -> Estimate the solution error for a scalar indicator variable
!#        in the H1-semi norm, that is, the error of the solution gradient
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

  public :: errest_initFromParameterlist
  public :: errest_releaseErrorEstimator
  public :: errest_addVariable
  public :: errest_clearErrorEstimator
  public :: errest_calcH1Error
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

  ! Goal-oriented error estimation
  integer, parameter, public :: ERREST_GOALORIENTED  = 8

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
    integer :: ierrorestimator = 0

    ! Type of element
    integer :: ieltype = 0

    ! Number of error variable
    integer :: ierrorvariable = 0

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
    real(DP) :: dnoiseFilter = 0._DP

    ! Absolute tolerance for filter
    real(DP) :: dabsFilter = 0._DP

    ! Scalar indicator variable
    type(t_vectorScalar) :: rerrorVariable

    ! Pointer to the triangulation structure
    type(t_triangulation), pointer :: p_rtriangulation => null()
  end type t_differenceIndicator
!</typeblock>

!<typeblock>

  ! This data structure contains the complete error estimator
  type t_errorEstimator

    ! Number of currecntly attached variables
    integer :: nvariables = 0

    ! Type of error estimator
    integer :: ierrorestimator = 0

    ! Type of grid-indicator
    integer :: igridindicator = 0
    
    ! Number of protection layers
    integer :: nprotectlayers = 0
    
    ! Pointer to the underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! A pointer to a list of handles of double precision pointers. 
    ! p_Hvariables(I) points to the data of variable I.
    integer, dimension(:), pointer :: p_Hvariables => null()
    

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

  subroutine errest_initFromParameterlist(rparlist, ssectionName, rerrorEstimator)

!<description>
    ! This subroutine initializes the error estimator
    ! with the values supplied by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! name of the section
    character(LEN=*), intent(IN) :: ssectionName
!</input>

!<output>
    ! error estimator
    type(t_errorEstimator), intent(OUT) :: rerrorEstimator
!</output>
!</subroutine>

    integer :: nerrorVariables

    ! Set the type of error estimator
    call parlst_getvalue_int(rparlist, ssectionName, 'ierrorestimator',&
                             rerrorEstimator%ierrorestimator)
    
    ! Set the type of grid indicator
    call parlst_getvalue_int(rparlist, ssectionName, 'igridindicator',&
                             rerrorEstimator%igridindicator, ERREST_ASIS)
    
    ! Set the number of protection layers
    call parlst_getvalue_int(rparlist, ssectionName, 'nprotectlayers',&
                             rerrorEstimator%nprotectlayers, 0)
    
!!$   
!!$
!!$    ! What kind of error estimator are we?
!!$    select case(rerrorEstimator%ierrorestimator)
!!$
!!$    case (ERREST_CSPR_FACE:ERREST_LIMAVR)
!!$
!!$      ! Recovery-based error estimator
!!$      allocate(rerrorEstimator%p_rgradientRecovery)
!!$      
!!$      ! Set the type of gradient recovery
!!$      rerrorEstimator%p_rgradientRecovery%ierrorestimator = rerrorEstimator%ierrorestimator
!!$
!!$      ! Set the error variable
!!$      call parlst_getvalue_int(rparlist, ssectionName, "ierrorvariable",&
!!$                               rerrorEstimator%p_rgradientRecovery%ierrorvariable, 1)
!!$
!!$      ! Set the type of element
!!$      call parlst_getvalue_int(rparlist, ssectionName, "ieltype",&
!!$                               rerrorEstimator%p_rgradientRecovery%ieltype)
!!$
!!$      !-------------------------------------------------------------------------
!!$
!!$    case (ERREST_FIRSTDIFF:ERREST_SECONDDIFF)
!!$
!!$      ! First-/Second-difference indicator
!!$      allocate(rerrorEstimator%p_rdifferenceIndicator)
!!$      
!!$      ! Set the type of difference indicator
!!$      rerrorEstimator%p_rdifferenceIndicator%ierrorestimator = rerrorEstimator%ierrorestimator
!!$      
!!$      ! Set the error variable
!!$      call parlst_getvalue_int(rparlist, ssectionName, "ierrorvariable",&
!!$                               rerrorEstimator%p_rdifferenceIndicator%ierrorvariable, 1)
!!$
!!$      ! Set the value for the noise filter
!!$      call parlst_getvalue_double(rparlist, ssectionName, "dnoiseFilter",&
!!$                                  rerrorEstimator%p_rdifferenceIndicator%dnoiseFilter, 5e-3_DP)
!!$
!!$      ! Set the value for the absolute filter
!!$      call parlst_getvalue_double(rparlist, ssectionName, "dabsFilter",&
!!$                                  rerrorEstimator%p_rdifferenceIndicator%dabsFilter, 1e-6_DP)
!!$      
!!$    case DEFAULT
!!$      call output_line('Invalid type of error estimator!',&
!!$          OU_CLASS_ERROR,OU_MODE_STD,'errest_initErrorEstimator')
!!$      call sys_halt()
!!$    end select

  end subroutine errest_initFromParameterlist

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
    
!!$    ! Release temporal data from error estimator
!!$    call errest_clearErrorEstimator(rerrorEstimator)
!!$    
!!$    ! Reset data
!!$    rerrorEstimator%ierrorestimator = 0
!!$    rerrorEstimator%igridindicator  = 0
!!$    rerrorEstimator%nprotectlayers  = 0
!!$
!!$    ! Release substructures
!!$    if (associated(rerrorEstimator%p_rgradientRecovery)) then
!!$
!!$      ! Reset data
!!$      rerrorEstimator%p_rgradientRecovery%ierrorestimator = 0
!!$      rerrorEstimator%p_rgradientRecovery%ieltype         = 0
!!$      rerrorEstimator%p_rgradientRecovery%ierrorvariable  = 0
!!$
!!$      ! Deallocate memory
!!$      deallocate(rerrorEstimator%p_rgradientRecovery)
!!$      nullify(rerrorEstimator%p_rgradientRecovery)
!!$    end if
!!$
!!$    if (associated(rerrorEstimator%p_rdifferenceIndicator)) then
!!$
!!$      ! Reset data
!!$      rerrorEstimator%p_rdifferenceIndicator%dnoiseFilter    = 0._DP
!!$      rerrorEstimator%p_rdifferenceIndicator%dabsFilter      = 0._DP
!!$      rerrorEstimator%p_rdifferenceIndicator%ierrorestimator = 0
!!$      rerrorEstimator%p_rdifferenceIndicator%ierrorvariable  = 0
!!$
!!$      ! Deallocate memory
!!$      deallocate(rerrorEstimator%p_rdifferenceIndicator)
!!$      nullify(rerrorEstimator%p_rdifferenceIndicator)
!!$    end if

  end subroutine errest_releaseErrorEstimator

  !*****************************************************************************

!<subroutine>

  subroutine errest_addVariable(rerrorEstimator, Ddata)

!<description>
    ! This subroutine adds variable data to the error estimator
!</description>

!<input>
    ! Ddata(I) is the value of the variable in vertex I of the triangulation)
    real(DP), dimension(:), intent(IN) :: Ddata
!</input>

!<inputoutput>
    ! error estimator
    type(t_errorEstimator), intent(INOUT) :: rerrorEstimator
!</inputoutput>
!</subroutine>

  end subroutine errest_addVariable
    
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

  subroutine errest_calcH1Error(rtriangulation, rvectorScalar, ierrorEstimator,&
                                rerrorH1, derrorH1Opt)

!<description>
    ! This subroutine estimates the error in the H1-semi-norm by
    ! means of various gradient recovery techniques.
!</description>
    
!<input>
    ! triangulation
    type(t_triangulation), intent(IN) :: rtriangulation
    
    ! FE solution vector
    type(t_vectorScalar), intent(IN) :: rvectorScalar

    ! type of error estimator
    integer, intent(IN) :: ierrorEstimator
!</input>

!<output>
    ! local H1-error (per element)
    type(t_vectorScalar), intent(OUT) :: rerrorH1

    ! OPTIONAL: global H1-error
    real(DP), intent(OUT), optional :: derrorH1Opt
!</output>
!</subroutine>
    
    ! How many spatial dimensions are we?
    select case (rtriangulation%ndim)
    case (NDIM1D)
      call calcH1Error1D(rtriangulation, rvectorScalar, ierrorEstimator,&
                         rerrorH1, derrorH1Opt)
      
    case (NDIM2D)
      call calcH1Error2D(rtriangulation, rvectorScalar, ierrorEstimator,&
                         rerrorH1, derrorH1Opt)
      
    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'errest_calcH1Error')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Calculate the H1 error in 1D
    
    subroutine calcH1Error1D(rtriangulation, rvectorScalar, ierrorEstimator,&
                             rerrorH1, derrorH1Opt)

      type(t_triangulation), intent(IN) :: rtriangulation
      type(t_vectorScalar), intent(IN) :: rvectorScalar
      integer, intent(IN) :: ierrorEstimator
      
      type(t_vectorScalar), intent(OUT) :: rerrorH1
      real(DP), intent(OUT), optional :: derrorH1Opt

      ! Pointer to the spatial discretization
      type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

      ! local variables
      type(t_blockDiscretisation) :: rdiscrBlock,rdiscrBlockRef
      type(t_vectorBlock) :: rgradient,rgradientRef
      real(DP) :: derrorH1
      

      ! Create new vector for local H1-error
      call lsyssc_createVector(rerrorH1, rtriangulation%NEL, .true.)

      ! What error estimator should be used?
      select case(ierrorEstimator)
      case (ERREST_CSPR_FACE:ERREST_SPR_FACE)
        !------------------------------------------------------------------------
        ! Superconvergent patch recovery
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
        
        ! Initialise block discretisations
        call spdiscr_initBlockDiscr(rdiscrBlock,    1, rtriangulation)
        call spdiscr_initBlockDiscr(rdiscrBlockRef, 1, rtriangulation)
        
!!$        ! What kind of element type is used for the FE solution?
!!$        select case(rgradientRecovery%ieltype)
!!$        case(-1,1,11)
!!$          ! Initialise spatial discretisations for gradient with P0-elements
!!$          call spdiscr_deriveSimpleDiscrSc(p_rspatialDiscr,&
!!$              EL_E000_1D, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
!!$          
!!$          ! Initialise spatial discretisations for reference gradient with P1-elements
!!$          call spdiscr_deriveSimpleDiscrSc(p_rspatialDiscr,&
!!$              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
!!$          
!!$        case(2,13)
!!$          ! Initialise spatial discretisations for gradient with P1-elements
!!$          call spdiscr_deriveSimpleDiscrSc(p_rspatialDiscr,&
!!$              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
!!$          
!!$          ! Initialise spatial discretisations for reference gradient with P2-elements
!!$          call spdiscr_deriveSimpleDiscrSc(p_rspatialDiscr,&
!!$              EL_E002_1D, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
!!$                    
!!$        case DEFAULT
!!$          call output_line('Unsupproted element type!',&
!!$                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error1D')
!!$          call sys_halt()
!!$        end select
        
        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr(rdiscrBlock,    rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr(rdiscrBlockRef, rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient (rvectorScalar, rgradient)
        
        ! How should the smoothed gradient be recovered?
        select case(ierrorEstimator)
        case (ERREST_L2PROJECTION)
          call ppgrd_calcGradient (rvectorScalar, rgradientRef)
          
        case (ERREST_SPR_VERTEX)
          call ppgrd_calcGradSuperPatchRecov (rvectorScalar, rgradientRef, PPGRD_NODEPATCH)
          
        case (ERREST_SPR_ELEMENT)
          call ppgrd_calcGradSuperPatchRecov (rvectorScalar, rgradientRef, PPGRD_ELEMPATCH)
          
        case (ERREST_SPR_FACE)
          call ppgrd_calcGradSuperPatchRecov (rvectorScalar, rgradientRef, PPGRD_FACEPATCH)
        end select
        
        ! Compute estimated gradient error
        call pperr_blockErrorEstimate(rgradient, rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
        ! Release temporal discretizations
        call spdiscr_releaseBlockDiscr(rdiscrBlock)
        call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
        
        ! Release temporal vectors
        call lsysbl_releaseVector(rgradient)
        call lsysbl_releaseVector(rgradientRef)

        
      case (ERREST_LIMAVR)
        !------------------------------------------------------------------------
        ! Limited gradient averaging
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
        
        ! Initialise block discretisations
        call spdiscr_initBlockDiscr(rdiscrBlock,    1, rtriangulation)
        call spdiscr_initBlockDiscr(rdiscrBlockRef, 1, rtriangulation)
        
!!$        ! What kind of element type is used for the FE solution?
!!$        select case(rgradientRecovery%ieltype)
!!$        case(-1,1,11)
!!$          ! Initialise spatial discretisations for gradient with P0-elements
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E000_1D, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
!!$          
!!$          ! Initialise spatial discretisations for reference gradient with P1-elements
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E001_1D, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
!!$                    
!!$        case DEFAULT
!!$          call output_line('Unsupported element type!',&
!!$                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error1D')
!!$          call sys_halt()
!!$        end select
        
        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr(rdiscrBlock,    rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr(rdiscrBlockRef, rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient(rvectorScalar, rgradient)
        
        ! Recover limited averaged gradient
        call ppgrd_calcGradLimAvgP1Q1cnf(rvectorScalar, rgradientRef)
        
        ! Compute gradient error
        call pperr_blockErrorEstimate(rgradient, rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
        ! Release temporal discretizations
        call spdiscr_releaseBlockDiscr(rdiscrBlock)
        call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
        
        ! Release temporal vectors
        call lsysbl_releaseVector(rgradient)
        call lsysbl_releaseVector(rgradientRef)

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
    
    subroutine calcH1Error2D(rtriangulation, rvectorScalar, ierrorEstimator,&
                             rerrorH1, derrorH1Opt)

      type(t_triangulation), intent(IN) :: rtriangulation
      type(t_vectorScalar), intent(IN) :: rvectorScalar
      integer, intent(IN) :: ierrorEstimator

      type(t_vectorScalar), intent(OUT) :: rerrorH1
      real(DP), intent(OUT), optional :: derrorH1Opt

      ! local variables
      type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
      type(t_blockDiscretisation) :: rdiscrBlock,rdiscrBlockRef
      type(t_vectorBlock) :: rgradient,rgradientRef
      real(DP) :: derrorH1

      integer :: i,idim,ccub,ieltype

      ! Create new vector for local H1-error
      call lsyssc_createVector(rerrorH1, rtriangulation%NEL, .true.)

      ! What error estimator should be used?
      select case(ierrorestimator)
      case (ERREST_CSPR_FACE:ERREST_SPR_FACE)
        !------------------------------------------------------------------------
        ! Superconvergent patch recovery
        !------------------------------------------------------------------------
        
        ! Set pointer
        p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
        
        ! Initialise block discretisations for the consistent/reconstructed gradient vectors
        call spdiscr_initBlockDiscr(rdiscrBlock, p_rspatialDiscr%ndimension, rtriangulation)
        call spdiscr_initBlockDiscr(rdiscrBlockRef, p_rspatialDiscr%ndimension, rtriangulation)
        
        ! Duplicate the discretisation from the scalar vector and adjust 
        ! the FE spaces for the consistent finite element gradient
        do idim = 1, p_rspatialDiscr%ndimension

          call spdiscr_duplicateDiscrSc(p_rspatialDiscr, rdiscrBlock%RspatialDiscr(idim))
          call spdiscr_duplicateDiscrSc(p_rspatialDiscr, rdiscrBlockRef%RspatialDiscr(idim))
          
          ! Adjust the FE space for the consistent gradient values
          do i = 1, rdiscrBlock%RspatialDiscr(idim)%inumFESpaces
            
            select case(rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%celement)
            case (EL_P1_1D)
              ieltype = EL_P0_1D
            case (EL_P2_1D)
              ieltype = EL_P1_1D
              
            case (EL_P1_2D)
              ieltype = EL_P0_2D
            case (EL_P2_2D)
              ieltype = EL_P1_2D
            case (EL_P3_2D)
              ieltype = EL_P2_2D

            case (EL_Q1_2D)
              ieltype = EL_Q0_2D
            case (EL_Q2_2D)
              ieltype = EL_Q1_2D
            case (EL_Q3_2D)
              ieltype = EL_Q2_2D

            case (EL_P1_3D)
              ieltype = EL_P0_3D
            case (EL_P2_3D)
              ieltype = EL_P1_3D

            case (EL_Q1_3D)
              ieltype = EL_Q0_3D
            case (EL_Q2_3D)
              ieltype = EL_Q1_3D

            case DEFAULT
              call output_line('Unsupproted element type!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error2D')
              call sys_halt()
            end select
            
            ! Compute natural cubature rule
            ccub = spdiscr_getStdCubature(ieltype)

            ! Adjust element distribution
            rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%celement        = ieltype
            rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeBilForm = ccub
            rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeLinForm = ccub
            rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeEval    = ccub
            rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ctrafoType      = elem_igetTrafoType(ieltype)
          end do
        end do

        ! Create block vector for gradient values
        call lsysbl_createVecBlockByDiscr (rdiscrBlock,    rgradient,    .true.)
        call lsysbl_createVecBlockByDiscr (rdiscrBlockRef, rgradientRef, .true.)
        
        ! Recover consistent gradient
        call ppgrd_calcGradient (rvectorScalar, rgradient)
        
        ! How should the smoothed gradient be recovered?
        select case(ierrorestimator)
        case (ERREST_L2PROJECTION)
          call ppgrd_calcGradient (rvectorScalar, rgradientRef)
          
        case (ERREST_SPR_VERTEX)
          call ppgrd_calcGradSuperPatchRecov (rvectorScalar, rgradientRef, PPGRD_NODEPATCH)
          
        case (ERREST_SPR_ELEMENT)
          call ppgrd_calcGradSuperPatchRecov (rvectorScalar, rgradientRef, PPGRD_ELEMPATCH)
          
        case (ERREST_SPR_FACE)
          call ppgrd_calcGradSuperPatchRecov (rvectorScalar, rgradientRef, PPGRD_FACEPATCH)
        end select
        
        ! Compute estimated gradient error
        call pperr_blockErrorEstimate(rgradient, rgradientRef, PPERR_L2ERROR,&
                                      derrorH1, relementError=rerrorH1)
        
        
!!$      case (ERREST_LIMAVR)
!!$        !------------------------------------------------------------------------
!!$        ! Limited gradient averaging
!!$        !------------------------------------------------------------------------
!!$        
!!$        ! Set pointer
!!$        p_rspatialDiscr => rgradientRecovery%rerrorVariable%p_rspatialDiscr
!!$        
!!$        ! Initialise block discretisations
!!$        call spdiscr_initBlockDiscr (rgradientRecovery%rdiscrBlock, 2 ,&
!!$                                     rgradientRecovery%p_rtriangulation,&
!!$                                     rgradientRecovery%p_rboundary)
!!$        call spdiscr_initBlockDiscr (rgradientRecovery%rdiscrBlockRef, 2 ,&
!!$                                     rgradientRecovery%p_rtriangulation,&
!!$                                     rgradientRecovery%p_rboundary)
!!$        
!!$        ! What kind of element type is used for the FE solution?
!!$        select case(rgradientRecovery%ieltype)
!!$        case(1)
!!$          ! Initialise spatial discretisations for gradient with P0-elements
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E000, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E000, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
!!$          
!!$          ! Initialise spatial discretisations for reference gradient with P1~-elements
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E020, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E020, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
!!$          
!!$        case(11)
!!$          ! Initialise spatial discretisations for gradient with Q0-elements
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E010, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E010, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
!!$          
!!$          ! Initialise spatial discretisations for reference gradient with Q1~-elements
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E030, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
!!$          call spdiscr_deriveSimpleDiscrSc (p_rspatialDiscr,&
!!$              EL_E030, SPDISC_CUB_AUTOMATIC,&
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
!!$          
!!$        case(-1)
!!$          ! Initialise spatial discretisations for gradient with P0/Q0-elements
!!$          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
!!$              EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(1))
!!$          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
!!$              EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
!!$              rgradientRecovery%rdiscrBlock%RspatialDiscr(2))
!!$          
!!$          ! Initialise spatial discretisations for reference gradient with P1~/Q1~-elements
!!$          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
!!$              EL_E020, EL_E030, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(1))
!!$          call spdiscr_deriveDiscr_triquad (p_rspatialDiscr,&
!!$              EL_E020, EL_E030, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
!!$              rgradientRecovery%rdiscrBlockRef%RspatialDiscr(2))
!!$          
!!$        case DEFAULT
!!$          call output_line('Unsupported element type!',&
!!$                           OU_CLASS_ERROR,OU_MODE_STD,'calcH1Error2D')
!!$          call sys_halt()
!!$        end select
!!$        
!!$        ! Create block vector for gradient values
!!$        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlock,&
!!$                                           rgradientRecovery%rgradient,    .true.)
!!$        call lsysbl_createVecBlockByDiscr (rgradientRecovery%rdiscrBlockRef,&
!!$                                           rgradientRecovery%rgradientRef, .true.)
!!$        
!!$        ! Recover consistent gradient
!!$        call ppgrd_calcGradient (rgradientRecovery%rerrorVariable,&
!!$                                 rgradientRecovery%rgradient)
!!$        
!!$        ! Recover limited averaged gradient
!!$        call ppgrd_calcGradLimAvgP1Q1cnf (rgradientRecovery%rerrorVariable,&
!!$                                          rgradientRecovery%rgradientRef)
!!$        
!!$        ! Compute gradient error
!!$        call pperr_blockErrorEstimate(rgradientRecovery%rgradient,&
!!$                                      rgradientRecovery%rgradientRef, PPERR_L2ERROR,&
!!$                                      derrorH1, relementError=rerrorH1)
        
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

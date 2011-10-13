!#########################################################################
!# ***********************************************************************
!# <name> pprocgradients </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating gradients
!# of finite element functions.
!#
!# The following routines can be found in this module:
!#
!# 1.) ppgrd_calcGradient
!#     -> Standard implementation for the calculation of the recovered 
!#        gradient of a scalar finite element function.
!#        A parameter admits to choose a method which is used with
!#        default parameters.
!#
!# 2.) ppgrd_calcGradientError
!#     -> Calculate the gradient error by comparing the consistent FE
!#        smoothed gradient obtained by gradient reconstruction
!#
!# .) ppgrd_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# Auxiliary routines, called internally.
!#
!# 1.) ppgrd_calcGradInterpP12Q12cnf
!#     -> Calculate the reconstructed gradient as P1, Q1, P2 or Q2
!#        vector for an arbitrary conformal discretisation.
!#        Uses 1st order interpolation to reconstruct the gradient.
!#
!# 2.) ppgrd_calcGradSuperPatchRecov
!#     -> Calculate the reconstructed gradient as P1, Q1, P2 or Q2
!#        vector for an arbitrary conformal discretisation.
!#        Uses the superconvergent patch recovery technique suggested
!#        by Zienkiewicz and Zhu. 1D, 2D, 3D.
!#
!# 3.) ppgrd_calcGradLimAvgP1Q1cnf
!#     -> Calculate the reconstructed gradient as ~P1 or ~Q1 vector
!#        for an arbitrary conformal discretisation.
!#        Uses the limited gradient averaging technique by M. Möller
!#        and D. Kuzmin.
!# </purpose>
!#########################################################################

module pprocgradients

  use basicgeometry
  use boundary
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
  use pprocerror
  use spatialdiscretisation
  use storage
  use storage
  use transformation
  use triangulation

  implicit none

  private

!<constants>

!<constantblock description = "Identifiers for the method how to calculate a gradient vector.">

  ! Use standard interpolation to calculate a gradient vector. 1st order. 
  integer, parameter, public :: PPGRD_INTERPOL = 0
  
  ! ZZ-technique for recovering a gradient. 2nd order on regular meshes.
  integer, parameter, public :: PPGRD_ZZTECHNIQUE = 1
  
  ! Limited gradient averaging technique. Only usable for special-type
  ! meshes consisting of only P1/Q1 finite elements.
  integer, parameter, public :: PPGRD_LATECHNIQUE = 2
  
!</constantblock>

!<constantblock description = "Identifiers for the type of patch used to recover the gradient vector.">

  ! Node-based patch: Use elements surrounding a particular node
  integer, parameter, public :: PPGRD_NODEPATCH = 0

  ! Element-based patch: Use elements surrounding a particular element
  integer, parameter, public :: PPGRD_ELEMPATCH = 1

  ! Face-based patch: Use subset of element-based patch which has common face
  integer, parameter, public :: PPGRD_FACEPATCH = 2

!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of elements to handle simultaneously when building vectors
#ifndef PPGRD_NELEMSIM
  integer, parameter, public :: PPGRD_NELEMSIM   = 1000
#endif

  ! Number of patches to handle simultaneously when performing gradient recovery
#ifndef PPGRD_NPATCHSIM
  integer, parameter, public :: PPGRD_NPATCHSIM  = 100
#endif
  
!</constantblock>

!</constants>

  !************************************************************************
  
  ! global performance configuration
  type(t_perfconfig), target, save :: ppgrd_perfconfig

  !************************************************************************

  public :: ppgrd_initPerfConfig
  public :: ppgrd_calcGradient
  public :: ppgrd_calcGradientError
  public :: ppgrd_calcGradInterpP12Q12cnf
  public :: ppgrd_calcGradSuperPatchRecov
  public :: ppgrd_calcGradLimAvgP1Q1cnf

contains

  !****************************************************************************

!<subroutine>

  subroutine ppgrd_initPerfConfig(rperfconfig)

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
      ppgrd_perfconfig = rperfconfig
    else
      ppgrd_perfconfig%NELEMSIM = PPGRD_NELEMSIM
      ppgrd_perfconfig%NPATCHSIM = PPGRD_NPATCHSIM
    end if
  
  end subroutine ppgrd_initPerfConfig

  !****************************************************************************

!<subroutine>

  subroutine ppgrd_calcGradient (rvectorScalar,rvectorGradient,cgradType,cgradSubtype)

!<description>
  ! Calculates the recovered gradient of a scalar finite element function.
  ! cgradType decides about the method to use for the calculation.
  ! This parameter is optional; if not specified, a standard method
  ! will be taken.
  !
  ! rvectorGradient receives the reconstructed gradient. For a 2D discretisation,
  ! this must be a 2D vector. For a 3D discretisation, this must be a 3D vector.
  ! The vector must provide a discretisation structure that defines the
  ! finite element space the reconstructed gradient should be calculated in.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(in) :: rvectorScalar
  
  ! OPTIONAL: Identifier for the method to use for calculating the gradient.
  ! One of the PPGRD_xxxx constants.
  ! If not specified, PPGRD_INTERPOL is taken as default.
  integer, intent(in), optional :: cgradType

  ! OPTIONAL: Identifier for the submethod to use for calculating the gradient.
  ! One of the PPGRD_xxxx constants.
  ! If not specified, PPGRD_NODEPATCH is taken as default.
  integer, intent(in), optional :: cgradSubtype
!</input>

!<inputoutput>
  ! A block vector receiving the gradient.
  ! The first subvector receives the X-gradient.
  ! In 2D/3D discretisations, the 2nd subvector recevies the Y-gradient.
  ! In 3D discretisations, the 3rd subvector receives the Z-gradient.
  ! The vector must be prepared with a discretisation structure that defines
  ! the destination finite element space for the gradient field.
  type(t_vectorBlock), intent(inout) :: rvectorGradient
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: imethod,isubmethod

    ! Some basic checks:
    
    imethod = PPGRD_INTERPOL
    if (present(cgradType)) imethod = cgradType

    isubmethod = PPGRD_NODEPATCH
    if (present(cgradSubtype)) isubmethod = cgradSubtype
    
    
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line ('No discretisation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      call sys_halt()
    end if
    
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) .and.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .ne. SPDISC_CONFORMAL)) then
      call output_line ('Only uniform and conformal discretisations supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      call sys_halt()
    end if
    
    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
      call output_line ('No triangulation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      call sys_halt()
    end if
    
    if (rvectorScalar%p_rspatialDiscr%p_rtriangulation%ndim .gt. &
        rvectorGradient%nblocks) then
      call output_line ('Dimension of destination vector not large enough!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      call sys_halt()
    end if
    
    ! Depending on the method choosen, call the appropriate gradient
    ! recovery routine.
    select case (imethod)
    case (PPGRD_INTERPOL)
      ! 1st order gradient
      call ppgrd_calcGradInterpP12Q12cnf (rvectorScalar,rvectorGradient)
      
    case (PPGRD_ZZTECHNIQUE)
      ! 2nd order gradient with ZZ.
      ! Standard method is 'nodewise'.
      call ppgrd_calcGradSuperPatchRecov (rvectorScalar,rvectorGradient,isubmethod)

    case (PPGRD_LATECHNIQUE)
      ! 1st order gradient
      call ppgrd_calcGradLimAvgP1Q1cnf(rvectorScalar, rvectorGradient)
      
    case DEFAULT
      call output_line ('Unsupported gradient recovery technique!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradient')
      call sys_halt()
    end select

  end subroutine ppgrd_calcGradient

  !****************************************************************************

!<subroutine>

  subroutine ppgrd_calcGradientError (rvectorScalar, derror, cgradType,&
                                      cgradSubtype, rerror)

!<description>
  ! Calculates the recovered gradient of a scalar finite element function
  ! and the consistent finite element gradient and compares both values
  ! to estimate the error between the approximate and the true gradient.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(in) :: rvectorScalar
  
  ! OPTIONAL: Identifier for the method to use for calculating the gradient.
  ! One of the PPGRD_xxxx constants.
  ! If not specified, PPGRD_INTERPOL is taken as default.
  integer, intent(in), optional :: cgradType

  ! OPTIONAL: Identifier for the submethod to use for calculating the gradient.
  ! One of the PPGRD_xxxx constants.
  ! If not specified, PPGRD_NODEPATCH is taken as default.
  integer, intent(in), optional :: cgradSubtype
!</input>

!<inputoutput>
  ! OPTIONAL
  type(t_vectorScalar), intent(inout), optional :: rerror
!</inputoutput>

!<output>
  ! The calculated error.
  real(DP), intent(out) :: derror
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_blockDiscretisation) :: rdiscrBlock,rdiscrBlockRef
    type(t_vectorBlock) :: rgradient,rgradientRef
    integer :: imethod,isubmethod,idim,i
    integer(I32) :: celement,ccub

    ! Some basic checks:
    
    imethod = PPGRD_INTERPOL
    if (present(cgradType)) imethod = cgradType

    isubmethod = PPGRD_NODEPATCH
    if (present(cgradSubtype)) isubmethod = cgradSubtype


    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line ('No discretisation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
      call sys_halt()
    end if
    
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) .and.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .ne. SPDISC_CONFORMAL)) then
      call output_line ('Only uniform and conformal discretisations supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
      call sys_halt()
    end if

    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
      call output_line ('No triangulation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
      call sys_halt()
    end if
    
    
    select case(imethod)
    case (PPGRD_INTERPOL, &
          PPGRD_ZZTECHNIQUE)
      
      ! Set pointer
      p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
      p_rtriangulation => p_rspatialDiscr%p_rtriangulation
      
      ! Initialise block discretisations for the consistent/reconstructed gradient vectors
      call spdiscr_initBlockDiscr(rdiscrBlock, p_rspatialDiscr%ndimension, p_rtriangulation)
      call spdiscr_initBlockDiscr(rdiscrBlockRef, p_rspatialDiscr%ndimension, p_rtriangulation)

      ! Duplicate the discretisation from the scalar vector and adjust 
      ! the FE spaces for the consistent finite element gradient
      do idim = 1, p_rspatialDiscr%ndimension
        
        call spdiscr_duplicateDiscrSc(p_rspatialDiscr, rdiscrBlock%RspatialDiscr(idim))
        call spdiscr_duplicateDiscrSc(p_rspatialDiscr, rdiscrBlockRef%RspatialDiscr(idim))
        
        ! Adjust the FE space for the consistent gradient values
        do i = 1, rdiscrBlock%RspatialDiscr(idim)%inumFESpaces
          
          select case(rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%celement)
          case (EL_P1_1D)
            celement = EL_P0_1D
          case (EL_P2_1D)
            celement = EL_P1_1D
            
          case (EL_P1_2D)
            celement = EL_P0_2D
          case (EL_P2_2D)
            celement = EL_P1_2D
          case (EL_P3_2D)
            celement = EL_P2_2D

          case (EL_Q1_2D)
            celement = EL_Q0_2D
          case (EL_Q2_2D)
            celement = EL_Q1_2D
          case (EL_Q3_2D)
            celement = EL_Q2_2D

          case (EL_P1_3D)
            celement = EL_P0_3D
          case (EL_P2_3D)
            celement = EL_P1_3D

          case (EL_Q1_3D)
            celement = EL_Q0_3D
          case (EL_Q2_3D)
            celement = EL_Q1_3D
            
          case DEFAULT
            call output_line('Unsupported element type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
            call sys_halt()
          end select
          
          ! Compute natural cubature rule
          ccub = spdiscr_getStdCubature(celement)
          
          ! Adjust element distribution
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%celement        = celement
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeBilForm = ccub
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeLinForm = ccub
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeEval    = ccub
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ctrafoType      = elem_igetTrafoType(celement)
        end do
      end do
      
      ! Create block vector for gradient values
      call lsysbl_createVecBlockByDiscr(rdiscrBlock, rgradient, .true.)
      call lsysbl_createVecBlockByDiscr(rdiscrBlockRef, rgradientRef, .true.)
      
      ! Recover consistent gradient vector
      call ppgrd_calcGradient (rvectorScalar, rgradient)

      ! Recover smoothed gradient vector
      if (imethod .eq. PPGRD_INTERPOL) then
        call ppgrd_calcGradInterpP12Q12cnf(rvectorScalar, rgradientRef)
      else
        call ppgrd_calcGradSuperPatchRecov(rvectorScalar, rgradientRef, isubmethod)
      end if

      ! Compute estimated gradient error
      call pperr_blockErrorEstimate(rgradient, rgradientRef, PPERR_L2ERROR,&
                                    derror, relementError=rerror)
          
      ! Release temporal discretizations
      call spdiscr_releaseBlockDiscr(rdiscrBlock)
      call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
      
      ! Release temporal vectors
      call lsysbl_releaseVector(rgradient)
      call lsysbl_releaseVector(rgradientRef)


    case (PPGRD_LATECHNIQUE)

      ! Set pointer
      p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
      p_rtriangulation => p_rspatialDiscr%p_rtriangulation
      
      ! Initialise block discretisations for the consistent/reconstructed gradient vectors
      call spdiscr_initBlockDiscr(rdiscrBlock, p_rspatialDiscr%ndimension, p_rtriangulation)
      call spdiscr_initBlockDiscr(rdiscrBlockRef, p_rspatialDiscr%ndimension, p_rtriangulation)

      ! Duplicate the discretisation from the scalar vector and adjust 
      ! the FE spaces for the consistent finite element gradient
      do idim = 1, p_rspatialDiscr%ndimension
        
        call spdiscr_duplicateDiscrSc(p_rspatialDiscr, rdiscrBlock%RspatialDiscr(idim))
        
        ! Adjust the FE space for the consistent gradient values
        do i = 1, rdiscrBlock%RspatialDiscr(idim)%inumFESpaces
          
          select case(rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%celement)
          case (EL_P1_1D)
            celement = EL_P0_1D
          case (EL_P2_1D)
            celement = EL_P1_1D
            
          case (EL_P1_2D)
            celement = EL_P0_2D
          case (EL_P2_2D)
            celement = EL_P1_2D
          case (EL_P3_2D)
            celement = EL_P2_2D

          case (EL_Q1_2D)
            celement = EL_Q0_2D
          case (EL_Q2_2D)
            celement = EL_Q1_2D
          case (EL_Q3_2D)
            celement = EL_Q2_2D

          case (EL_P1_3D)
            celement = EL_P0_3D
          case (EL_P2_3D)
            celement = EL_P1_3D

          case (EL_Q1_3D)
            celement = EL_Q0_3D
          case (EL_Q2_3D)
            celement = EL_Q1_3D
            
          case DEFAULT
            call output_line('Unsupported element type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
            call sys_halt()
          end select
          
          ! Compute natural cubature rule
          ccub = spdiscr_getStdCubature(celement)
          
          ! Adjust element distribution
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%celement        = celement
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeBilForm = ccub
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeLinForm = ccub
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeEval    = ccub
          rdiscrBlock%RspatialDiscr(idim)%RelementDistr(i)%ctrafoType      = elem_igetTrafoType(celement)
        end do
      end do

      ! Duplicate the discretisation from the scalar vector and adjust 
      ! the FE spaces for the reconstructed finite element gradient
      do idim = 1, p_rspatialDiscr%ndimension
        
        call spdiscr_duplicateDiscrSc(p_rspatialDiscr, rdiscrBlockRef%RspatialDiscr(idim))
        
        ! Adjust the FE space for the consistent gradient values
        do i = 1, rdiscrBlockRef%RspatialDiscr(idim)%inumFESpaces
          
          select case(rdiscrBlockRef%RspatialDiscr(idim)%RelementDistr(i)%celement)
          case (EL_P1_2D)
            celement = EL_P1T_2D
         
          case (EL_Q1_2D)
            celement = EL_Q1T_2D

          case (EL_Q1_3D)
            celement = EL_Q1T_3D
            
          case DEFAULT
            call output_line('Unsupported element type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
            call sys_halt()
          end select
          
          ! Compute natural cubature rule
          ccub = spdiscr_getStdCubature(celement)
          
          ! Adjust element distribution
          rdiscrBlockRef%RspatialDiscr(idim)%RelementDistr(i)%celement        = celement
          rdiscrBlockRef%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeBilForm = ccub
          rdiscrBlockRef%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeLinForm = ccub
          rdiscrBlockRef%RspatialDiscr(idim)%RelementDistr(i)%ccubTypeEval    = ccub
          rdiscrBlockRef%RspatialDiscr(idim)%RelementDistr(i)%ctrafoType      = elem_igetTrafoType(celement)
        end do
      end do
      
      ! Create block vector for gradient values
      call lsysbl_createVecBlockByDiscr(rdiscrBlock, rgradient, .true.)
      call lsysbl_createVecBlockByDiscr(rdiscrBlockRef, rgradientRef, .true.)
      
      ! Recover consistent gradient vector
      call ppgrd_calcGradient (rvectorScalar, rgradient)

      ! Recover smoothed gradient vector
      call ppgrd_calcGradLimAvgP1Q1cnf(rvectorScalar, rgradientRef)

      ! Compute estimated gradient error
      call pperr_blockErrorEstimate(rgradient, rgradientRef, PPERR_L2ERROR,&
                                    derror, relementError=rerror)
          
      ! Release temporal discretizations
      call spdiscr_releaseBlockDiscr(rdiscrBlock)
      call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
      
      ! Release temporal vectors
      call lsysbl_releaseVector(rgradient)
      call lsysbl_releaseVector(rgradientRef)


    case DEFAULT
      call output_line ('Unsupported gradient recovery technique!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradientError')
      call sys_halt()
    end select
      

  end subroutine ppgrd_calcGradientError

  !****************************************************************************

!<subroutine>

  subroutine ppgrd_calcGradInterpP12Q12cnf (rvectorScalar,rvectorGradient)

!<description>
  ! Calculates the recovered gradient of a scalar finite element function
  ! by standard interpolation. Supports conformal discretisations
  ! with <tex>$P_1$, $Q_1$, $P_2$ and $Q_2$</tex> mixed in the destination vector.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(in) :: rvectorScalar
!</input>

!<inputoutput>
  ! A block vector receiving the gradient.
  ! The first subvector receives the X-gradient.
  ! In 2D/3D discretisations, the 2nd subvector recevies the Y-gradient.
  ! In 3D discretisations, the 3rd subvector receives the Z-gradient.
  ! The vector must be prepared with a discretisation structure that defines
  ! the destination finite element space for the gradient field.
  type(t_vectorBlock), intent(inout) :: rvectorGradient
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,k,icurrentElementDistr,NVE,IELmax,IELset
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, 4) :: Dxi

    ! Cubature formula weights. The cotent is actually not used here.
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! Number of 'cubature points'. As we acutally do not do cubature
    ! here, this coincides with the number of DOF`s on each element
    ! in the destination space.
    integer :: nlocalDOFsDest
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofDest
    
    ! Type of transformation from the reference to the real element 
    integer(I32) :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! Current element distribution in source- and destination vector
    type(t_elementDistribution), pointer :: p_relementDistribution
    type(t_elementDistribution), pointer :: p_relementDistrDest
    
    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dderivatives
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), allocatable :: IdofsTrial
    integer, dimension(:,:), allocatable :: IdofsDest
    
    ! Pointers to the X- and Y-derivative vector
    real(DP), dimension(:), pointer :: p_DxDeriv, p_DyDeriv, p_DzDeriv
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to an array that counts the number of elements adjacent to a vertex.
    ! Ok, there is the same information in the triangulation, but that is not
    ! based on DOF`s! Actually, we will calculate how often we touch each DOF 
    ! in the destination space.
    integer :: h_IcontributionsAtDOF
    integer, dimension(:), pointer :: p_IcontributionsAtDOF
    
    ! Discretisation structures for the source- and destination vector(s)
    type(t_spatialDiscretisation), pointer :: p_rdiscrSource, p_rdiscrDest
     

    ! Loop over all blocks of the gradient and over all FE spaces
    do i = 1, min(rvectorGradient%nblocks,&
                  rvectorScalar%p_rspatialDiscr%ndimension)
      do j = 1, rvectorGradient%p_rblockDiscr%RspatialDiscr(i)%inumFESpaces
        
        select case(elem_getPrimaryElement(&
            rvectorGradient%p_rblockDiscr%RspatialDiscr(i)%RelementDistr(j)%celement))
        case (EL_Q0_2D, EL_Q0_3D,&
              EL_P0_1D, EL_P0_2D, EL_P0_3D,&
              EL_Q1_2D, EL_Q1_3D,&
              EL_P1_1D, EL_P1_2D, EL_P1_3D,&
              EL_Q2_2D, EL_Q2_3D,&
              EL_P2_1D, EL_P2_2D, EL_P2_3D)
          
        case DEFAULT
          call output_line ('Only Q0, Q1, Q2, P0, P1, and P2 supported as' // &
              ' discretisation for the destination vector!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradInterpP12Q12cnf')
          call sys_halt()
        end select
        
      end do
    end do

    ! Get the discretisation structures of the source- and destination space.
    ! Note that we assume here that all derivatives are discretised the same way!
    p_rdiscrSource => rvectorScalar%p_rspatialDiscr
    p_rdiscrDest => rvectorGradient%p_rblockDiscr%RspatialDiscr(1)
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscrSource%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPGRD_NELEMSIM,p_rtriangulation%NEL)
    
    ! Array that allows the calculation about the number of elements
    ! meeting in a vertex, based onm DOF`s.
    call storage_new ('ppgrd_calcGradInterpP1Q1cnf','DOFContrAtVertex',&
                      dof_igetNDofGlob(p_rdiscrDest),&
                      ST_INT, h_IcontributionsAtDOF, ST_NEWBLOCK_ZERO)
    call storage_getbase_int (h_IcontributionsAtDOF,p_IcontributionsAtDOF)
    
    ! Evaluate the first derivative of the FE functions.
    Bder = .false.

    ! Get pointers to the derivative destination vector.
    select case(p_rtriangulation%ndim)
    case (NDIM1D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      Bder(DER_DERIV1D_X) = .true.

    case (NDIM2D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      call lalg_clearVectorDble (p_DyDeriv)
      Bder(DER_DERIV2D_X) = .true.
      Bder(DER_DERIV2D_Y) = .true.

    case (NDIM3D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(3),p_DzDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      call lalg_clearVectorDble (p_DyDeriv)
      call lalg_clearVectorDble (p_DzDeriv)
      Bder(DER_DERIV3D_X) = .true.
      Bder(DER_DERIV3D_Y) = .true.
      Bder(DER_DERIV3D_Z) = .true.

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      call sys_halt()
    end select
    
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,p_rdiscrSource%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => p_rdiscrSource%RelementDistr(icurrentElementDistr)
      p_relementDistrDest => p_rdiscrDest%RelementDistr(icurrentElementDistr)
    
      ! If the element distribution is empty, skip it
      if (p_relementDistribution%NEL .eq. 0) cycle
    
      ! Get the number of local DOF`s for trial functions
      ! in the source and destination vector.
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      indofDest = elem_igetNDofLoc(p_relementDistrDest%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Initialise the cubature formula,
      ! That is a special trick here! The FE space of the destination vector
      ! is either P1 or Q1. We create the gradients in the corners of these points
      ! by taking the mean of the gradients of the source vector!
      !
      ! For this purpose, we initialise the 'trapezoidal rule' as cubature
      ! formula. Ok, we are not using cubature at all here (thus ignoring
      ! the cubature weights completely), but that way we get the
      ! coordinates of the corners on the reference element automatically --
      ! which coincide with the points where we want to create the gradients!
      !
      ! Note: The returned nlocalDOFsDest will coincide with the number of local DOF`s
      ! on each element indofDest!
      select case (elem_getPrimaryElement(p_relementDistrDest%celement))
      case (EL_P0_1D)
        call cub_getCubPoints(CUB_G1_1D, nlocalDOFsDest, Dxi, Domega)
      case (EL_P1_1D)
        call cub_getCubPoints(CUB_TRZ_1D, nlocalDOFsDest, Dxi, Domega)
      case (EL_P2_1D)
        ! Manually calculate the coordinates of the 
        ! corners/midpoints on the reference element.
        Dxi(1,1) = 0.0_DP
        Dxi(2,1) = 1.0_DP
        Dxi(3,1) = 0.5_DP

        nlocalDOFsDest = 3

      case (EL_P0_2D)
        call cub_getCubPoints(CUB_G1_T, nlocalDOFsDest, Dxi, Domega)
      case (EL_Q0_2D)
        call cub_getCubPoints(CUB_G1_2D, nlocalDOFsDest, Dxi, Domega)
      case (EL_P1_2D)
        call cub_getCubPoints(CUB_TRZ_T, nlocalDOFsDest, Dxi, Domega)
      case (EL_Q1_2D)
        call cub_getCubPoints(CUB_TRZ, nlocalDOFsDest, Dxi, Domega)
      case (EL_P2_2D)
        ! Manually calculate the coordinates of the 
        ! corners/midpoints on the reference element.
        Dxi(1,1)  =  1.0_DP
        Dxi(1,2)  =  0.0_DP
        Dxi(1,3)  =  0.0_DP

        Dxi(2,1)  =  0.0_DP
        Dxi(2,2)  =  1.0_DP
        Dxi(2,3)  =  0.0_DP
        
        Dxi(3,1)  =  0.0_DP
        Dxi(3,2)  =  0.0_DP
        Dxi(3,3)  =  1.0_DP
        
        Dxi(4,1)  =  0.5_DP
        Dxi(4,2)  =  0.5_DP
        Dxi(4,3)  =  0.0_DP

        Dxi(5,1)  =  0.0_DP
        Dxi(5,2)  =  0.5_DP
        Dxi(5,3)  =  0.5_DP

        Dxi(6,1)  =  0.5_DP
        Dxi(6,2)  =  0.0_DP
        Dxi(6,3)  =  0.5_DP
        
        nlocalDOFsDest = 6
        
      case (EL_Q2_2D)
        ! Manually calculate the coordinates of the 
        ! corners/midpoints on the reference element.
        Dxi(1,1)  =  -1.0_DP
        Dxi(1,2)  =  -1.0_DP

        Dxi(2,1)  =  1.0_DP
        Dxi(2,2)  =  -1.0_DP
        
        Dxi(3,1)  =  1.0_DP
        Dxi(3,2)  =  1.0_DP
        
        Dxi(4,1)  =  -1.0_DP
        Dxi(4,2)  =  1.0_DP

        Dxi(5,1)  =  0.0_DP
        Dxi(5,2)  =  -1.0_DP

        Dxi(6,1)  =  1.0_DP
        Dxi(6,2)  =  0.0_DP

        Dxi(7,1)  =  0.0_DP
        Dxi(7,2)  =  1.0_DP

        Dxi(8,1)  =  -1.0_DP
        Dxi(8,2)  =  0.0_DP
        
        Dxi(9,1)  =  0.0_DP
        Dxi(9,2)  =  0.0_DP
        
        nlocalDOFsDest = 9

      case (EL_P0_3D)
        call cub_getCubPoints(CUB_G1_3D_T, nlocalDOFsDest, Dxi, Domega)
      case (EL_Q0_3D)
        call cub_getCubPoints(CUB_G1_3D, nlocalDOFsDest, Dxi, Domega)
      case (EL_P1_3D)
        call cub_getCubPoints(CUB_TRZ_3D_T, nlocalDOFsDest, Dxi, Domega)
      case (EL_Q1_3D)
        call cub_getCubPoints(CUB_TRZ_3D, nlocalDOFsDest, Dxi, Domega)
      case (EL_P2_3D)
        ! Manually calculate the coordinates of the 
        ! corners/midpoints on the reference element.
        Dxi(1,1)  =  1.0_DP
        Dxi(1,2)  =  0.0_DP
        Dxi(1,3)  =  0.0_DP
        Dxi(1,4)  =  0.0_DP

        Dxi(2,1)  =  0.0_DP
        Dxi(2,2)  =  1.0_DP
        Dxi(2,3)  =  0.0_DP
        Dxi(2,4)  =  0.0_DP
        
        Dxi(3,1)  =  0.0_DP
        Dxi(3,2)  =  0.0_DP
        Dxi(3,3)  =  1.0_DP
        Dxi(3,4)  =  0.0_DP

        Dxi(4,1)  =  0.0_DP
        Dxi(4,2)  =  0.0_DP
        Dxi(4,3)  =  0.0_DP
        Dxi(4,4)  =  1.0_DP
        
        Dxi(5,1)  =  0.5_DP
        Dxi(5,2)  =  0.5_DP
        Dxi(5,3)  =  0.0_DP
        Dxi(5,4)  =  0.0_DP

        Dxi(6,1)  =  0.5_DP
        Dxi(6,2)  =  0.0_DP
        Dxi(6,3)  =  0.5_DP
        Dxi(6,4)  =  0.0_DP

        Dxi(7,1)  =  0.5_DP
        Dxi(7,2)  =  0.0_DP
        Dxi(7,3)  =  0.0_DP
        Dxi(7,4)  =  0.5_DP

        Dxi(8,1)  =  0.0_DP
        Dxi(8,2)  =  0.5_DP
        Dxi(8,3)  =  0.5_DP
        Dxi(8,4)  =  0.0_DP

        Dxi(9,1)  =  0.0_DP
        Dxi(9,2)  =  0.5_DP
        Dxi(9,3)  =  0.0_DP
        Dxi(9,4)  =  0.5_DP

        Dxi(10,1) =  0.0_DP
        Dxi(10,2) =  0.0_DP
        Dxi(10,3) =  0.5_DP
        Dxi(10,4) =  0.5_DP

        nlocalDOFsDest = 10

      case (EL_Q2_3D)
        ! Manually calculate the coordinates of the 
        ! corners/midpoints on the reference element.
        Dxi(1,1)  = -1.0_DP
        Dxi(1,2)  = -1.0_DP
        Dxi(1,3)  = -1.0_DP
        
        Dxi(2,1)  =  1.0_DP
        Dxi(2,2)  = -1.0_DP
        Dxi(2,3)  = -1.0_DP
        
        Dxi(3,1)  =  1.0_DP
        Dxi(3,2)  =  1.0_DP
        Dxi(3,3)  = -1.0_DP
        
        Dxi(4,1)  = -1.0_DP
        Dxi(4,2)  =  1.0_DP
        Dxi(4,3)  = -1.0_DP
        
        Dxi(5,1)  = -1.0_DP
        Dxi(5,2)  = -1.0_DP
        Dxi(5,3)  =  1.0_DP
        
        Dxi(6,1)  =  1.0_DP
        Dxi(6,2)  = -1.0_DP
        Dxi(6,3)  =  1.0_DP
        
        Dxi(7,1)  =  1.0_DP
        Dxi(7,2)  =  1.0_DP
        Dxi(7,3)  =  1.0_DP
        
        Dxi(8,1)  = -1.0_DP
        Dxi(8,2)  =  1.0_DP
        Dxi(8,3)  =  1.0_DP
                
        Dxi(9,1)  =  0.0_DP
        Dxi(9,2)  = -1.0_DP
        Dxi(9,3)  = -1.0_DP

        Dxi(10,1) =  1.0_DP
        Dxi(10,2) =  0.0_DP
        Dxi(10,3) = -1.0_DP

        Dxi(11,1)  =  0.0_DP
        Dxi(11,2)  =  1.0_DP
        Dxi(11,3)  = -1.0_DP

        Dxi(12,1)  = -1.0_DP
        Dxi(12,2)  =  0.0_DP
        Dxi(12,3)  = -1.0_DP
        
        Dxi(13,1) =  0.0_DP
        Dxi(13,2) = -1.0_DP
        Dxi(13,3) =  1.0_DP

        Dxi(14,1) =  1.0_DP
        Dxi(14,2) =  0.0_DP
        Dxi(14,3) =  1.0_DP

        Dxi(15,1)  =  0.0_DP
        Dxi(15,2)  =  1.0_DP
        Dxi(15,3)  =  1.0_DP

        Dxi(16,1)  = -1.0_DP
        Dxi(16,2)  =  0.0_DP
        Dxi(16,3)  =  1.0_DP

        Dxi(17,1)  = -1.0_DP
        Dxi(17,2)  = -1.0_DP
        Dxi(17,3)  =  0.0_DP
        
        Dxi(18,1)  =  1.0_DP
        Dxi(18,2)  = -1.0_DP
        Dxi(18,3)  =  0.0_DP

        Dxi(19,1)  =  1.0_DP
        Dxi(19,2)  =  1.0_DP
        Dxi(19,3)  =  0.0_DP

        Dxi(20,1)  = -1.0_DP
        Dxi(20,2)  =  1.0_DP
        Dxi(20,3)  =  0.0_DP

        nlocalDOFsDest = 20
              
      case DEFAULT
        call output_line('Unsupported FE space in destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradInterpP1Q1cnf')
        call sys_halt()
      end select

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,nlocalDOFsDest
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do
      
      ! Allocate memory for the DOF`s of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))
      allocate(IdofsDest(indofDest,nelementsPerBlock))

      ! Allocate memory for the values of the derivatives in the corners
      allocate(Dderivatives(nlocalDOFsDest, nelementsPerBlock, p_rtriangulation%ndim))

      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      ! Get the number of elements in the element distribution.
      NEL = p_relementDistribution%NEL
      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPGRD_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPGRD_NELEMSIM)
      
        ! Calculate the global DOF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscrSource, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Also calculate the global DOF`s in our destination vector(s)
        call dof_locGlobMapping_mult(p_rdiscrDest, p_IelementList(IELset:IELmax), &
                                     IdofsDest)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:nlocalDOFsDest))
            
        ! In the next loop, we do not have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we calculate the gradient information.
        select case(p_rtriangulation%ndim)
        case (NDIM1D)
          ! Calculate the X-derivative in the corners of the elements
          ! into Dderivatives(:,:,1).
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
              p_relementDistribution%celement, IdofsTrial, DER_DERIV1D_X,&
              Dderivatives(:,1:IELmax-IELset+1,1))

          ! Sum up the derivative values in the destination vector.
          ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
          ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
          ! local DOF`s in the destination space -- in that order!
          
          do i=1,IELmax-IELset+1
            do j=1,nlocalDOFsDest
              p_DxDeriv(IdofsDest(j,i)) = p_DxDeriv(IdofsDest(j,i)) + Dderivatives(j,i,1)
              
              ! Count how often a DOF was touched.
              p_IcontributionsAtDOF(IdofsDest(j,i)) = &
                  p_IcontributionsAtDOF(IdofsDest(j,i))+1
           end do
         end do

        case (NDIM2D)
          ! Calculate the X-derivative in the corners of the elements
          ! into Dderivatives(:,:,1) and the Y-derivative into Dderivatives(:,:,2).
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
              p_relementDistribution%celement, IdofsTrial, DER_DERIV2D_X,&
              Dderivatives(:,1:IELmax-IELset+1,1))

          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
              p_relementDistribution%celement, IdofsTrial, DER_DERIV2D_Y,&
              Dderivatives(:,1:IELmax-IELset+1,2))

          ! Sum up the derivative values in the destination vector.
          ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
          ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
          ! local DOF`s in the destination space -- in that order!
          
          do i=1,IELmax-IELset+1
            do j=1,nlocalDOFsDest
              p_DxDeriv(IdofsDest(j,i)) = p_DxDeriv(IdofsDest(j,i)) + Dderivatives(j,i,1)
              p_DyDeriv(IdofsDest(j,i)) = p_DyDeriv(IdofsDest(j,i)) + Dderivatives(j,i,2)
              
              ! Count how often a DOF was touched.
              p_IcontributionsAtDOF(IdofsDest(j,i)) = &
                  p_IcontributionsAtDOF(IdofsDest(j,i))+1
           end do
         end do

         case (NDIM3D)
           ! Calculate the X-derivative in the corners of the elements
           ! into Dderivatives(:,:,1), the Y-derivative into Dderivatives(:,:,2)
           ! and the Z-derivative into Dderivatives(:,:,3)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
              p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_X,&
              Dderivatives(:,1:IELmax-IELset+1,1))

          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
              p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_Y,&
              Dderivatives(:,1:IELmax-IELset+1,2))

          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
              p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_Z,&
              Dderivatives(:,1:IELmax-IELset+1,3))
          
          ! Sum up the derivative values in the destination vector.
          ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
          ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
          ! local DOF`s in the destination space -- in that order!
          
          do i=1,IELmax-IELset+1
            do j=1,nlocalDOFsDest
              p_DxDeriv(IdofsDest(j,i)) = p_DxDeriv(IdofsDest(j,i)) + Dderivatives(j,i,1)
              p_DyDeriv(IdofsDest(j,i)) = p_DyDeriv(IdofsDest(j,i)) + Dderivatives(j,i,2)
              p_DzDeriv(IdofsDest(j,i)) = p_DzDeriv(IdofsDest(j,i)) + Dderivatives(j,i,3)
              
              ! Count how often a DOF was touched.
              p_IcontributionsAtDOF(IdofsDest(j,i)) = &
                  p_IcontributionsAtDOF(IdofsDest(j,i))+1
           end do
         end do

       end select
       
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dderivatives)
      deallocate(IdofsDest)
      deallocate(IdofsTrial)

    end do ! icurrentElementDistr

    ! We are nearly done. The final thing: divide the calculated derivatives by 
    ! the number of elements adjacent to each vertex. That closes the calculation
    ! of the 'mean' of the derivatives.
    select case(p_rtriangulation%ndim)
    case (NDIM1D)
      do i=1,size(p_DxDeriv)
        ! Div/0 should not occur, otherwise the triangulation is 
        ! crap as there is a point not connected to any element!
        p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
      end do
      
    case (NDIM2D)
      do i=1,size(p_DxDeriv)
        ! Div/0 should not occur, otherwise the triangulation is
        ! crap as there is a point not connected to any element!
        p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
        p_DyDeriv(i) = p_DyDeriv(i) / p_IcontributionsAtDOF(i)
      end do

    case (NDIM3D)
      do i=1,size(p_DxDeriv)
        ! Div/0 should not occur, otherwise the triangulation is
        ! crap as there is a point not connected to any element!
        p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
        p_DyDeriv(i) = p_DyDeriv(i) / p_IcontributionsAtDOF(i)
        p_DzDeriv(i) = p_DzDeriv(i) / p_IcontributionsAtDOF(i)
      end do
      
    end select
    
    ! Release temp data
    call storage_free (h_IcontributionsAtDOF)

  end subroutine ppgrd_calcGradInterpP12Q12cnf

  !****************************************************************************

!<subroutine>

  subroutine ppgrd_calcGradSuperPatchRecov (rvectorScalar,rvectorGradient,cpatchType)

!<description>
    ! Calculates the recovered gradient of a scalar finite element function
    ! by means of the superconvergent patch recovery technique suggested
    ! by Zienkiewicz and Zhu. Supports conformal discretisations in arbitrary
    ! spatial dimensions with <tex>$P_1$, $Q_1$, $P_2$ and $Q_2$</tex> finite elements
    ! mixed in the source and destination vectors.
!</description>

!<input>
    ! The FE solution vector. Represents a scalar FE function.
    type(t_vectorScalar), intent(in) :: rvectorScalar
    
    ! The type of patch used to recover the gradient values
    integer, intent(in) :: cpatchType
!</input>

!<inputoutput>
    ! A block vector receiving the gradient.
    ! The first subvector receives the X-gradient.
    ! In 2D/3D discretisations, the 2nd subvector recevies the Y-gradient.
    ! In 3D discretisations, the 3rd subvector receives the Z-gradient.
    ! The vector must be prepared with a discretisation structure that defines
    ! the destination finite element space for the gradient field.
    type(t_vectorBlock), intent(inout) :: rvectorGradient
!</inputoutput>

!</subroutine>
    
    ! local variables
    real(DP), dimension(NDIM3D) :: Dval
    integer :: IVE,NVE,NVEMax,IEL,JEL,KEL,IVT,IPATCH,NPATCH,PATCHset,PATCHmax
    integer :: icurrentElementDistr,ilastElementDistr,ilocalElementDistr
    integer :: i,j,k,ipoint,idx,idx2,idxsubgroup
    integer(i32) :: icoordSystem
    logical :: bnonparTrial


    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder,BderDest

    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega

    ! Flag if discretisation is uniform
    logical :: bisuniform
    
    ! Number of cubature points on the reference element
    integer :: ncubp

    ! Maximum number of cubature points on the reference element in the block
    integer :: ncubpMax
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial

    ! Maximum numer of local degrees of freedom for test functiions in the block
    integer :: indofTrialMax

    ! Number of local degrees of freedom for destination space
    integer :: indofDest

    ! Maximum number of local degrees of freedom for destination space
    integer :: indofDestMax

    ! Number of 'cubature points'. As we acutally do not do cubature
    ! here, this coincides with the number of DOF`s on each element
    ! in the destination space.
    integer :: nlocalDOFsDest

    ! Maximum number of 'cubature points'. As we acutally do not do 
    ! cubature here, this coincides with the number of DOF`s on each 
    ! element in the destination space.
    integer :: nlocalDOFsDestMax

    ! Total number of sampling/interpolation points in set of patches
    integer :: nspoints,nipoints

    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation

    ! The spatial discretisation structure - to shorten some things...
    type(t_spatialDiscretisation), pointer :: p_rdiscrSource
    type(t_spatialDiscretisation), pointer :: p_rdiscrDest

    ! Current element distribution in use
    type(t_elementDistribution), pointer :: p_relementDistribution
    type(t_elementDistribution), pointer :: p_relementDistrDest

    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_domainIntSubset) :: rintSubsetDest

    ! An allocatable array accepting the starting positions of each patch
    integer, dimension(:), allocatable :: IelementsInPatchIdx

    ! An allocatable array accepting the element numbers of each patch
    ! Note that the first member for each patch defines the patch details, i.e.,
    ! the cubature formular that should be used for this patch, etc.
    integer, dimension(:), allocatable :: IelementsInPatch

    ! An allocatable array accepting the number of corners per element
    integer, dimension(:), allocatable :: IelementNVEInPatch

    ! An allocatable array accepting the number of cubature points per element
    integer, dimension(:), allocatable :: IelementNcubpInPatch

    ! An allocatable array accepting the number of sampling points per patch
    integer, dimension(:), allocatable :: Inpoints
    
    ! An allocatable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), allocatable :: IdofsTrial
    integer, dimension(:,:), allocatable :: IdofsDest

    ! An allocatable array accepting the coordinates of the patch bounding group
    real(DP), dimension(:,:,:), allocatable :: DpatchBound

    ! An allocatable array accepting the polynomials of the set of elements.
    real(DP), dimension(:,:,:,:), allocatable :: Dpolynomials
    real(DP), dimension(:),       allocatable :: DpolynomialsMixed

    ! An allocatable array accepting the values of the FE functions
    ! that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    real(DP), dimension(:,:),   allocatable :: DcoefficientsMixed

    ! An allocatable array accepting the averaged gradient values
    real(DP), dimension(:,:,:), allocatable :: Dderivatives

    ! Pointer to element-at-vertex index array of the triangulation
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx

    ! Pointer to element-at-vertex list of the triangulation
    integer, dimension(:), pointer :: p_IelementsAtVertex

    ! Pointer to elements-at-element list of the triangulation
    integer, dimension(:,:), pointer :: p_IneighboursAtElement

     ! Pointer to vertices-at-element list of the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList

    ! Pointer to vertex coordinates of the triangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to element distribution identifier list.
    integer, dimension(:), pointer :: p_IelementDistr

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsRef

    ! An array receiving the coordinates of cubature points on
    ! the real element for all elements in a set.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsReal

    ! Pointer to the point coordinates to pass to the element function.
    ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
    ! the trial element is parametric or not.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsTrial

    ! Array with coordinates of the corners that form the real element.
    real(DP), dimension(:,:,:), pointer :: p_Dcoords
    
    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:,:,:), pointer :: p_Djac      

    ! Number of patches in a block
    integer :: npatchesPerBlock

    ! Number of patches currently blocked
    integer :: npatchesInCurrentBlock

    ! Number of elements in a block
    integer :: nelementsPerBlock

    ! Pointer to an array that counts the number of elements adjacent to a vertex.
    ! Ok, there is the same information in the triangulation, but that is not
    ! based on DOF`s! Actually, we will calculate how often we touch each DOF 
    ! in the destination space.
    integer :: h_IcontributionsAtDOF
    integer, dimension(:), pointer :: p_IcontributionsAtDOF
    
    ! Pointers to the X-, Y- and Z-derivative vector
    real(DP), dimension(:), pointer :: p_DxDeriv, p_DyDeriv, p_DzDeriv

    ! Auxiliary integers
    integer :: icubp,idim
    

    ! Loop over all blocks of the gradient and over all FE spaces
    do i = 1, min(rvectorGradient%nblocks,&
                  rvectorScalar%p_rspatialDiscr%ndimension)
      do j = 1, rvectorGradient%p_rblockDiscr%RspatialDiscr(i)%inumFESpaces
        
        select case(elem_getPrimaryElement(&
            rvectorGradient%p_rblockDiscr%RspatialDiscr(i)%RelementDistr(j)%celement))
        case (EL_Q0_2D, EL_Q0_3D,&
              EL_P0_1D, EL_P0_2D, EL_P0_3D,&
              EL_Q1_2D, EL_Q1_3D,&
              EL_P1_1D, EL_P1_2D, EL_P1_3D,&
              EL_Q2_2D, EL_Q2_3D,&
              EL_P2_1D, EL_P2_2D, EL_P2_3D)
          
        case DEFAULT
          call output_line ('Only Q0, Q1, Q2, P0, P1, and P2 supported as' // &
              ' discretisation for the destination vector!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
          call sys_halt()
        end select
        
      end do
    end do

    !-----------------------------------------------------------------------
    ! (0)  Initialisation
    !-----------------------------------------------------------------------
       
    ! Get the discretisation structures of the source- and destination space.
    ! Note that we assume here that all derivatives are discretised the same way!
    p_rdiscrSource => rvectorScalar%p_rspatialDiscr
    p_rdiscrDest   => rvectorGradient%p_rblockDiscr%RspatialDiscr(1)
    
    ! Check if we have a non-uniform discretisation structure and if nodal 
    ! patches should be used. This does not make too much sense since it is
    ! not clear which type of element should be adopted for the patch elements.
    !
    ! Theoretically, one could allow for an optional parameter which defines
    ! the element type of the "patch" elements but this requires more complicated
    ! code. There are other possibilities for the patches that should be used !!!
    if ((p_rdiscrSource%ccomplexity .ne. SPDISC_UNIFORM .or. &
         p_rdiscrDest%ccomplexity   .ne. SPDISC_UNIFORM) .and. &
         cpatchType .eq. PPGRD_NODEPATCH) then 
      call output_line('Nodal patches are not available for non-uniform discretisations!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      call sys_halt()
    end if

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscrSource%p_rtriangulation
    
    ! Get a pointer to the vertices-at-element and vertex coordinates array
    call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_double2D (p_rtriangulation%h_DvertexCoords,  p_DvertexCoords)

    ! We only need the derivatives of trial functions
    Bder = .false.

    ! Get pointers to the derivative destination vector.
    select case(p_rtriangulation%ndim)
    case (NDIM1D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      Bder(DER_DERIV1D_X) = .true.

    case (NDIM2D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      call lalg_clearVectorDble (p_DyDeriv)
      Bder(DER_DERIV2D_X) = .true.
      Bder(DER_DERIV2D_Y) = .true.

    case (NDIM3D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(3),p_DzDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      call lalg_clearVectorDble (p_DyDeriv)
      call lalg_clearVectorDble (p_DzDeriv)
      Bder(DER_DERIV3D_X) = .true.
      Bder(DER_DERIV3D_Y) = .true.
      Bder(DER_DERIV3D_Z) = .true.

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      call sys_halt()
    end select
      
    ! Array that allows the calculation about the number of elements
    ! meeting in a vertex, based onm DOF`s.
    call storage_new ('ppgrd_calcGradSuperPatchRecov','DOFContrAtVertex',&
        dof_igetNDofGlob(p_rdiscrDest), ST_INT, h_IcontributionsAtDOF, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(h_IcontributionsAtDOF, p_IcontributionsAtDOF)
    
    ! For the recovery, we only need the function values of trial functions
    BderDest = .false.
    BderDest(DER_FUNC) = .true.
    

    ! Do we have a uniform triangulation? Would simplify a lot...
    bisUniform = (p_rdiscrSource%ccomplexity .eq. SPDISC_UNIFORM) .and. &
                 (p_rdiscrDest%ccomplexity   .eq. SPDISC_UNIFORM)

    if (.not. bisUniform) then
      ! Things are more complicated if one of the discretisations is not uniform.
      ! In this case, we always have to consider the maximum number of quadrature
      ! points, the largest number of local DOF`s, etc. 

      indofTrialMax     = 0
      indofDestMax      = 0
      NVEmax            = 0
      ncubpMax          = 0
      nlocalDOFsDestMax = 0
      
      ! Set pointer to list of element distributions which exists in this case
      call storage_getbase_int(p_rdiscrSource%h_IelementDistr, p_IelementDistr)

      ! Loop over all element distributions of the source discretisation
      do icurrentElementDistr = 1, p_rdiscrSource%inumFESpaces
        
        ! Activate the current element distribution
        p_relementDistribution => p_rdiscrSource%RelementDistr(icurrentElementDistr)
        
        ! Cancel if this element distribution is empty.
        if (p_relementDistribution%NEL .eq. 0) cycle
        
        ! Get the number of local DOF`s for trial functions
        indofTrial    = elem_igetNDofLoc(p_relementDistribution%celement)
        indofTrialMax = max(indofTrialMax,indofTrial)
        
        ! Get the number of corner vertices of the element
        NVE    = elem_igetNVE(p_relementDistribution%celement)
        NVEmax = max (NVEmax,NVE)
        
        ! Get cubature weights and point coordinates on the reference element
        call cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
        ncubpMax = max(ncubpMax,ncubp)
      end do

      ! Loop over all element distributions of the destination discretisation
      do icurrentElementDistr = 1, p_rdiscrDest%inumFESpaces
        
        ! Activate the current element distribution
        p_relementDistribution => p_rdiscrDest%RelementDistr(icurrentElementDistr)

        ! Cancel if this element distribution is empty.
        if (p_relementDistribution%NEL .eq. 0) cycle

        ! Get the number of local DOF`s for trial functions
        indofDest    = elem_igetNDofLoc(p_relementDistribution%celement)
        indofDestMax = max(indofDestMax,indofDest)

        ! Get cubature weights and point coordinates on the reference element
        call cub_getCubPoints(p_relementDistribution%ccubTypeEval, nlocalDOFsDest, Dxi, Domega)
        nlocalDOFsDestMax = max(nlocalDOFsDestMax,nlocalDOFsDest)
      end do
    end if
        

    ! Set pointers which are used to assemble patches
    select case(cpatchType)
    case (PPGRD_NODEPATCH)
      ! Get the elements-at-vertex index array.
      call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
      ! Get the elements-at-vertex array.
      call storage_getbase_int (p_rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
            
    case (PPGRD_ELEMPATCH)
      ! Get the elements-at-vertex index array.
      call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
      ! Get the elements-at-vertex array.
      call storage_getbase_int (p_rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
      ! Get the elements-at-element array.
      call storage_getbase_int2D (p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
      ! Get the vertices-at-element array.
      call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)     

    case (PPGRD_FACEPATCH)
      ! Get the elements-at-element array.
      call storage_getbase_int2D (p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)

    case DEFAULT
      call output_line('Invalid patch type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! (1)  Superconvergent patch recovery - main loop
    !---------------------------------------------------------------------------
    
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation. Note that
    ! nodal patches are only allowed for uniform discretisations, that is,
    ! we can be sure, that there is exactly one element distribution if
    ! nodal patches should be considered.

    do icurrentElementDistr = 1, p_rdiscrSource%inumFESpaces

      ! Activate the current element distribution
      p_relementDistribution => p_rdiscrSource%RelementDistr(icurrentElementDistr)
      
      ! If the element distribution is empty, skip it
      if (p_relementDistribution%NEL .eq. 0) cycle
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, p_IelementList)
      
      ! Get number of patches in current element distribution
      select case(cpatchType)
      case (PPGRD_NODEPATCH)
        ! Recall that for nodal-baes patches there MUST be exactly one element 
        ! distribution so that the number of patches equals the total number
        ! of vertices in the whole triangulation !!!
        NPATCH = p_rtriangulation%NVT
        
      case (PPGRD_ELEMPATCH,PPGRD_FACEPATCH)
        ! The number of patches equals the number of elements 
        ! in the current element distribution.
        NPATCH = p_relementDistribution%NEL
        
      case DEFAULT
        call output_line('Invalid patch type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
        call sys_halt()
      end select

      ! Determine actual number of patches in this block
      npatchesPerBlock = min(PPGRD_NPATCHSIM,NPATCH)

      ! Allocate memory for element numbers in patch index array
      allocate(IelementsInPatchIdx(npatchesPerBlock+1))
      
      ! Allocate memory for coordinates of patch bounding groups
      allocate(DpatchBound(p_rtriangulation%ndim,2,npatchesPerBlock))

      ! Allocate memory for number of sampling points
      if (.not. bisUniform) allocate(Inpoints(npatchesPerBlock))
      
      ! Loop over the patches - blockwise.
      do PATCHset = 1, NPATCH , PPGRD_NPATCHSIM
        
        !-------------------------------------------------------------------------
        ! Phase 1: Determine all elements in the set of patches
        !-------------------------------------------------------------------------
        
        ! We always handle PPGRD_NPATCHSIM patches simultaneously.
        ! How many patches have we actually here?
        ! Get the maximum patch number, such that we handle at most 
        ! PPGRD_NPATCHSIM patches simultaneously.
        PATCHmax = min(NPATCH,PATCHset-1+PPGRD_NPATCHSIM)
        
        ! Calculate the number of patches currently blocked
        npatchesInCurrentblock = PATCHmax-PATCHset+1
        
        ! Depending on the patch type, the element numbers that are contained in
        ! each patch are different and must be determined from different structures
        select case(cpatchType)
        case (PPGRD_NODEPATCH)
          
          ! We actually know how many elements are adjacent to each node, so we can
          ! allocate memory for element numbers in patch without 'dummy' loop
          nelementsPerBlock = p_IelementsAtVertexIdx(PATCHmax+1)-p_IelementsAtVertexIdx(PATCHset)
          allocate(IelementsInPatch(nelementsPerBlock))
          
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Loop over the patches in the set
          do IPATCH = 1, npatchesInCurrentBlock
            
            ! Store index of first element in this patch
            IelementsInPatchIdx(IPATCH) = nelementsPerBlock + 1
            
            ! Get global vertex number
            IVT = PATCHset-1+IPATCH
            
            ! Loop over elements adjacent to vertex
            do idx = p_IelementsAtVertexIdx(IVT),p_IelementsAtVertexIdx(IVT+1)-1
              
              ! Get element number
              IEL = p_IelementsAtVertex(idx)
              
              ! Increase element counter
              nelementsPerBlock = nelementsPerBlock + 1
              IelementsInPatch(nelementsPerBlock) = IEL
            end do
          end do
          
          ! Store index of last element in last patch increased by one
          IelementsInPatchIdx(npatchesInCurrentBlock+1) = nelementsPerBlock + 1
         
          
        case (PPGRD_ELEMPATCH)
          
          ! Unfortunately, we do not directly know how many elements are in the neighbourhood
          ! of each element. But we know, how many elements are adjacent to each corner of
          ! a particular element. If we sum up these numbers we get an upper bound for the
          ! number of elements present in each patch. Obviously, the center element is multiply
          ! counted. Moreover, all edge/face neighbours are also counted twice.
          
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Ok, let us do a dummy loop to determine the number of elements in the block
          do IPATCH = 1, npatchesInCurrentBlock
            
            ! Get the global element number from the list of elements in distribution
            IEL = p_IelementList(PATCHset+IPATCH-1)
            
            ! Loop over corner nodes
            do ive =1, NVE
              
              ! Get global vertex number
              IVT = p_IverticesAtElement(ive,IEL)
              
              ! Count number of elements surrounding corner node
              nelementsPerBlock = nelementsPerBlock + &
                  (p_IelementsAtVertexIdx(IVT+1)-p_IelementsAtVertexIdx(IVT))
            end do
            
            ! Ok, we counted element IEL NVE-times but it is only required once
            nelementsPerBlock = nelementsPerBlock - (NVE-1)
            
            ! Moreover, each adjacent element is counted twice if it is not the boundary
            nelementsPerBlock = nelementsPerBlock - count(p_IneighboursAtElement(1:NVE,IEL) > 0)
          end do
          
          ! That is it, we can allocate memory for elements numbers
          allocate(IelementsInPatch(nelementsPerBlock))
          
          ! Now, we have to fill it with the element numbers
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Loop over the patches in the set
          do IPATCH = 1, npatchesInCurrentBlock
            
            ! Store index of first element in this patch
            IelementsInPatchIdx(IPATCH) = nelementsPerBlock + 1
            
            ! Get the global element number from the list of elements in distribution
            IEL = p_IelementList(PATCHset+IPATCH-1)
            
            ! Do not forget to store the element itself
            nelementsPerBlock = nelementsPerBlock + 1
            IelementsInPatch(nelementsPerBlock) = IEL
            
            ! Loop over corner nodes
            do ive = 1, NVE
              
              ! Get the globale vertex number
              IVT = p_IverticesAtElement(ive,IEL)
              
              ! Get element number adjacent to element IEL
              JEL = p_IneighboursAtElement(ive,IEL)
              
              ! Loop over elements adjacent to vertex
              do idx = p_IelementsAtVertexIdx(IVT),p_IelementsAtVertexIdx(IVT+1)-1
                
                ! Get element number
                KEL = p_IelementsAtVertex(idx)
                
                ! Do not consider KEL = IEL and KEL = JEL to prevent the same
                ! element number to be present in a patch multiple times
                if (KEL .eq. IEL .or. KEL .eq. JEL) cycle
                
                ! Increase element counter
                nelementsPerBlock = nelementsPerBlock + 1
                IelementsInPatch(nelementsPerBlock) = KEL
              end do
            end do
          end do
          
          ! Store index of last element in last patch increased by one
          IelementsInPatchIdx(npatchesInCurrentBlock+1) = nelementsPerBlock + 1
                    
        case (PPGRD_FACEPATCH)
          
          ! We actually know how many elements are adjacent to each element, so we
          ! can allocate memory for element numbers in patch without 'dummy' loop
          nelementsPerBlock = (NVE+1)*npatchesPerBlock
          allocate(IelementsInPatch(nelementsPerBlock))
          
          ! Initialise number of elements in block
          nelementsPerBlock = 0
          
          ! Loop over the patches in the set
          do IPATCH = 1, npatchesInCurrentBlock
            
            ! Store index of first element in this patch
            IelementsInPatchIdx(IPATCH) = nelementsPerBlock + 1
            
            ! Get the global element number from the list of elements in distribution
            IEL = p_IelementList(PATCHset+IPATCH-1)
            
            ! Do not forget to store the element itself
            nelementsPerBlock = nelementsPerBlock + 1
            IelementsInPatch(nelementsPerBlock) = IEL
            
            ! Loop over adjacent elements
            do ive = 1, NVE
              
              ! Get element number adjacent to element IEL
              JEL = p_IneighboursAtElement(ive,IEL)
              
              ! Check if element neighbour is the boundary, then skip it
              if (JEL .eq. 0) cycle
              
              ! Increase element counter
              nelementsPerBlock = nelementsPerBlock + 1
              IelementsInPatch(nelementsPerBlock) = JEL
            end do
          end do
          
          ! Store index of last element in last patch increased by one
          IelementsInPatchIdx(npatchesInCurrentBlock+1) = nelementsPerBlock + 1
          
        case DEFAULT
          call output_line('Invalid patch type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
          call sys_halt()
        end select
        
        
        !-----------------------------------------------------------------------
        ! Phase 2: Perform Least-squares fitting on the set of patches
        !-----------------------------------------------------------------------
        
        ! Do we have a uniform discretisation? Would simplify a lot...
        if (bisUniform) then
          
          ! Yes, the discretisation is uniform. In this case, we have to set the 
          ! pointers, dimensions, etc. just once and can work blockwise.
          
          ! Active element distribution
          p_relementDistribution => p_rdiscrSource%RelementDistr(icurrentElementDistr)
          p_relementDistrDest    => p_rdiscrDest%RelementDistr(icurrentElementDistr)
          
          ! Get the number of local DOF`s for trial functions
          indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
          indofDest  = elem_igetNDofLoc(p_relementDistrDest%celement)
          
          ! Get the number of corner vertices of the element
          NVE = elem_igetNVE(p_relementDistribution%celement)
          

          !---------------------------------------------------------------------
          ! Step 1:  Prepare the source FE space
          !---------------------------------------------------------------------
          
          ! Allocate memory for the DOF`s of all the elements
          allocate(IdofsTrial(indofTrial,nelementsPerBlock))
          
          ! Calculate the global DOF`s into IdofsTrial.
          call dof_locGlobMapping_mult(p_rdiscrSource, &
              IelementsInPatch(1:nelementsPerBlock), IdofsTrial)
          
          ! Initialise the cubature formula for the source element distribution
          ! Get cubature weights and point coordinates on the reference element
          call cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
          
          ! Get from the trial element space the type of coordinate system
          ! that is used there:
          icoordSystem = elem_igetCoordSystem(p_relementDistribution%celement)
          
          ! Allocate memory and get local references to it. This domain integration 
          ! structure stores all information of the source FE space. 
          call domint_initIntegration (rintSubset, nelementsPerBlock, &
            ncubp, icoordSystem, p_rtriangulation%ndim, NVE)
          p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
          p_DcubPtsReal => rintSubset%p_DcubPtsReal
          p_Djac =>        rintSubset%p_Djac
          p_Ddetj =>       rintSubset%p_Ddetj
          p_Dcoords =>     rintSubset%p_DCoords
          
          ! Put the cubature point coordinates in the right format to the
          ! cubature-point array.
          ! Initialise all entries in p_DcubPtsRef with the same coordinates -
          ! as the cubature point coordinates are identical on all elements
          do j=1,size(p_DcubPtsRef,3)
            do i=1,ncubp
              do k=1,size(p_DcubPtsRef,1)
                ! Could be solved using the TRANSPOSE operator - but often it is
                ! faster this way...
                p_DcubPtsRef(k,i,j) = Dxi(i,k)
              end do
            end do
          end do
          
          ! Check if one of the trial/test elements is nonparametric
          bnonparTrial = elem_isNonparametric(p_relementDistribution%celement)
          
          ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
          ! p_DcubPtsRef - depending on whether the space is parametric or not.
          if (bnonparTrial) then
            p_DcubPtsTrial => p_DcubPtsReal
          else
            p_DcubPtsTrial => p_DcubPtsRef
          end if
          
          ! We have the coordinates of the cubature points saved in the
          ! coordinate array from above. Unfortunately for nonparametric
          ! elements, we need the real coordinate.
          ! Furthermore, we anyway need the coordinates of the element
          ! corners and the Jacobian determinants corresponding to
          ! all the points.
          
          ! At first, get the coordinates of the corners of all the
          ! elements in the current set of elements.
          call trafo_getCoords_sim (&
              elem_igetTrafoType(p_relementDistribution%celement), &
              p_rtriangulation, IelementsInPatch(1:nelementsPerBlock), p_Dcoords)
          
          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          call trafo_calctrafo_sim (p_relementDistribution%ctrafoType, &
              nelementsPerBlock, ncubp, p_Dcoords, p_DcubPtsRef, p_Djac, &
              p_Ddetj, p_DcubPtsReal)
          
          !---------------------------------------------------------------------
          ! Step 2:  Perform sampling of consistent gradient values
          !---------------------------------------------------------------------
          
          ! Allocate memory for the values of the derivaties in the corners
          allocate(Dcoefficients(ncubp,nelementsPerBlock,p_rtriangulation%ndim))
          
          ! Calculate the derivative of the FE function in the cubature
          ! points: u_h(x,y,z) and save the result to Dcoefficients(:,:,1..3)
          select case(p_rtriangulation%ndim)
          case (NDIM1D)
            call fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_relementDistribution%celement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV1D_X, Dcoefficients(:,:,1))
            
          case (NDIM2D)
            call fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_relementDistribution%celement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV2D_X, Dcoefficients(:,:,1))
            
            call fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_relementDistribution%celement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV2D_Y, Dcoefficients(:,:,2))
            
          case (NDIM3D)
            call fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_relementDistribution%celement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV3D_X, Dcoefficients(:,:,1))
            
            call fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_relementDistribution%celement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV3D_Y, Dcoefficients(:,:,2))
            
            call fevl_evaluate_sim (rvectorScalar, p_Dcoords, p_Djac, p_Ddetj, &
                p_relementDistribution%celement, IdofsTrial, ncubp, &
                nelementsPerBlock, p_DcubPtsTrial, DER_DERIV3D_Z, Dcoefficients(:,:,3))

          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
            call sys_halt()
          end select
          

          !---------------------------------------------------------------------
          ! Step 3: Prepare least-squares fitting
          !---------------------------------------------------------------------
          
          ! First, get the coordinates of the bounding group for the set of
          ! elements present in each patch. In addition, calculate the coordinates
          ! of the corners of the constant Jacobian patch "elements".
          call calc_patchBoundingGroup_sim(IelementsInPatchIdx, npatchesInCurrentBlock, &
              NVE, p_Dcoords, DpatchBound(:,:,1:npatchesInCurrentBlock))
          
          ! Next, we need to convert the physical coordinates of the curbature points
          ! to the local coordinates of the constant Jacobian "patch" elements
          call calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, &
              npatchesInCurrentBlock, NVE, p_DcubPtsReal, DpatchBound, p_DcubPtsRef)
          
          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          call trafo_calctrafo_sim (p_relementDistribution%ctrafoType, &
              nelementsPerBlock, ncubp, p_Dcoords, p_DcubPtsRef, p_Djac, p_Ddetj)
          
          ! Allocate memory for the patch interpolants matrices
          ! Note that the second dimension is DER_FUNC=1 and could be omitted. However,
          ! all routines for simultaneous element evaluation require a 4d-array so that
          ! the array Dpolynomials is artificially created as 4d-array.
          allocate(Dpolynomials(indofTrial,DER_FUNC,ncubp,nelementsPerBlock))
          
          ! Evaluate the trial functions of the constant Jacobian patch "element" for all
          ! cubature points of the elements present in the patch and store each polynomial 
          ! interpolation in the rectangular patch matrix used for least-squares fitting.
          call elem_generic_sim(p_relementDistribution%celement,&
              p_Dcoords, p_Djac, p_Ddetj, BderDest, Dpolynomials, ncubp,&
              nelementsPerBlock, p_DcubPtsTrial)
          
          
          !-----------------------------------------------------------------------
          ! Step 4: Perform least-squares fitting
          !-----------------------------------------------------------------------
          
          ! Allocate memory for the derivative values
          allocate(Dderivatives(indofTrial,npatchesInCurrentBlock,p_rtriangulation%ndim))
          
          ! Compute the patch averages by solving $(P^T * P) * x = (P^T) * b$ for x
          call calc_patchAverages_sim(IelementsInPatchIdx, npatchesInCurrentBlock, ncubp, &
              indofTrial,  Dcoefficients, Dpolynomials,Dderivatives)
          

          !-----------------------------------------------------------------------
          ! Step 5: Prepare the destination FE space
          !-----------------------------------------------------------------------
          
          ! Allocate memory for the DOF`s of all the elements
          allocate(IdofsDest(indofDest,nelementsPerBlock))
          
          ! Also calculate the global DOF`s in our destination vector(s)
          call dof_locGlobMapping_mult(p_rdiscrDest, &
              IelementsInPatch(1:nelementsPerBlock), IdofsDest)
          
          ! Initialise the cubature formula. That is a special trick here!
          ! In particular, we only need the physical coordinates of the
          ! nodal evaluation points in the source vectors.
          ! For this purpose, we initialise the 'trapezoidal rule' as cubature
          ! formula. Ok, we are not using cubature at all here (thus ignoring
          ! the cubature weights completely), but that way we get the
          ! coordinates of the corners on the reference element automatically --
          ! which coincide with the points where we want to create the gradients!
          !
          ! Note: The returned nlocalDOFsDest will coincide with the number of local DOF`s
          ! on each element indofDest!
          call calc_cubatureDest(&
              elem_getPrimaryElement(p_relementDistrDest%celement), &
              nlocalDOFsDest, Dxi, Domega)
          
          ! Get from the trial element space the type of coordinate system
          ! that is used there:
          icoordSystem = elem_igetCoordSystem(p_relementDistrDest%celement)
          
          ! Allocate memory and get local references to it. This domain integration 
          ! structure stores all information of the destination FE space.
          call domint_initIntegration (rintSubsetDest, nelementsPerBlock, &
            nlocalDOFsDest, icoordSystem, p_rtriangulation%ndim, NVE)
          p_DcubPtsRef =>  rintSubsetDest%p_DcubPtsRef
          p_DcubPtsReal => rintSubsetDest%p_DcubPtsReal
          p_Djac =>        rintSubsetDest%p_Djac
          p_Ddetj =>       rintSubsetDest%p_Ddetj
          p_Dcoords =>     rintSubsetDest%p_DCoords
          
          ! Put the cubature point coordinates in the right format to the
          ! cubature-point array.
          ! Initialise all entries in p_DcubPtsRef with the same coordinates -
          ! as the cubature point coordinates are identical on all elements
          do j=1,size(p_DcubPtsRef,3)
            do i=1,nlocalDOFsDest
              do k=1,size(p_DcubPtsRef,1)
                ! Could be solved using the TRANSPOSE operator - but often it is
                ! faster this way...
                p_DcubPtsRef(k,i,j) = Dxi(i,k)
              end do
            end do
          end do
          
          ! Check if one of the trial/test elements is nonparametric
          bnonparTrial = elem_isNonparametric(p_relementDistrDest%celement)
          
          ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
          ! p_DcubPtsRef - depending on whether the space is parametric or not.
          if (bnonparTrial) then
            p_DcubPtsTrial => p_DcubPtsReal
          else
            p_DcubPtsTrial => p_DcubPtsRef
          end if
          
          ! We have the coordinates of the cubature points saved in the
          ! coordinate array from above. Unfortunately for nonparametric
          ! elements, we need the real coordinate.
          ! Furthermore, we anyway need the coordinates of the element
          ! corners and the Jacobian determinants corresponding to
          ! all the points.
          
          ! At first, get the coordinates of the corners of all the
          ! elements in the current set of elements.
          call trafo_getCoords_sim (elem_igetTrafoType(p_relementDistrDest%celement), &
              p_rtriangulation, IelementsInPatch(1:nelementsPerBlock), p_Dcoords)
          
          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          call trafo_calctrafo_sim (p_relementDistrDest%ctrafoType, &
              nelementsPerBlock, nlocalDOFsDest, p_Dcoords, p_DcubPtsRef, p_Djac, &
              p_Ddetj, p_DcubPtsReal)
          
          ! Next, we need to convert the physical coordinates of the curvature points
          ! to the local coordinates of the constant Jacobian "patch" elements
          call calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, &
              npatchesInCurrentBlock, NVE,  p_DcubPtsReal, DpatchBound, p_DcubPtsRef)
          
          ! We do not need the corner coordinates of the elements in the destination
          ! FE space but that of the "patch" elements in the source FE space
          p_Dcoords => rintSubset%p_Dcoords

          ! Calculate the transformation from the reference elements to the real ones
          call trafo_calctrafo_sim (p_relementDistrDest%ctrafoType, &
              nelementsPerBlock, nlocalDOFsDest, p_Dcoords, p_DcubPtsRef, p_Djac, &
              p_Ddetj)
          
          ! Reallocate memory for the patch interpolants
          if (ncubp .lt. nlocalDOFsDest) then
            deallocate(Dpolynomials)
            allocate(Dpolynomials(indofTrial,DER_FUNC,nlocalDOFsDest,nelementsPerBlock))
          end if
          
          ! Evaluate the basis functions for the cubature points of the destination FE space
          call elem_generic_sim(p_relementDistribution%celement,&
              p_Dcoords, p_Djac, p_Ddetj, BderDest, Dpolynomials, &
              nlocalDOFsDest, nelementsPerBlock, p_DcubPtsTrial)
          
          !---------------------------------------------------------------------
          ! Step 6: Evaluate the averaged derivative values at the cubature
          !         points of the destination FE space and scatter them to
          !         the global degrees of freedom
          !---------------------------------------------------------------------
          
          select case (p_rtriangulation%ndim)
          case(NDIM1D) 
            ! Loop over the patches in the set
            do ipatch = 1, npatchesInCurrentBlock

              ! Loop over elements in patch
              do idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
                
                ! Loop over local degrees of freedom
                do ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  do j = 1, indofTrial
                    Dval(1) = Dval(1) + Dderivatives(j,ipatch,1) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  end do
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                end do
              end do
            end do

          case (NDIM2D)
            ! Loop over the patches in the set
            do ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              do idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
                
                ! Loop over local degrees of freedom
                do ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  do j = 1, indofTrial
                    Dval(1:2) = Dval(1:2) + Dderivatives(j,ipatch,1:2) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  end do
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                end do
              end do
            end do

          case (NDIM3D)
            ! Loop over the patches in the set
            do ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              do idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1
                
                ! Loop over local degrees of freedom
                do ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  do j = 1, indofTrial
                    Dval(1:3) = Dval(1:3) + Dderivatives(j,ipatch,1:3) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  end do
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  p_DzDeriv(IdofsDest(ipoint,idx)) = p_DzDeriv(IdofsDest(ipoint,idx)) + Dval(3)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                end do
              end do
            end do
            
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
            call sys_halt()
          end select
          
          ! Release memory
          call domint_doneIntegration(rintSubset)
          call domint_doneIntegration(rintSubsetDest)
          
          ! Deallocate temporary memory
          deallocate(IdofsDest)
          deallocate(Dderivatives)
          deallocate(Dpolynomials)
          deallocate(Dcoefficients)
          deallocate(IdofsTrial)
          
        else
          
          ! No, the discretisation is not uniform. In this case, we have to check the type
          ! of element for each individual element and adopt the FE spaces accordingly.
          !
          ! In the beginning, perform a short bubble-sort to sort the elements in each 
          ! patch for their element distribution. 
          !
          ! Loop over the patches in the set
          do ipatch = 1, npatchesInCurrentBlock

            ! Get the element distribution of the first element in the patch
            IEL = IelementsInPatch(IelementsInPatchIdx(ipatch))
            ilastElementDistr = p_IelementDistr(IEL)
            
            ! Now find all elements in the patch that are in the same element
            ! distribution. Shift them to behind IEL. 
            !
            ! Loop over elements in patch
            bigsort: do idx = IelementsInPatchIdx(ipatch)+1,IelementsInPatchIdx(ipatch+1)-1
              
              ! Get global element number
              IEL = IelementsInPatch(idx)

              ! Get number of local element distribution
              ilocalElementDistr = p_IelementDistr(IEL)
              
              ! If it is different to the current element distribution,
              ! try to find another element in the patch which has that
              ! distruibution and can be shifted to the front.
              if (ilocalElementDistr .ne. ilastElementDistr) then
                
                do idx2 = idx+1,IelementsInPatchIdx(ipatch+1)-1
                  ! Is there another element with element distribution
                  ! ilastElementDistr?

                  ilocalElementDistr = p_IelementDistr(IelementsInPatch(idx2))
                  
                  if (ilocalElementDistr .eq. ilastElementDistr) then
                    ! Yep, we found one. Shift it to the front and continue
                    ! with the next element.
                    IEL = IelementsInPatch(idx2)
                    IelementsInPatch(idx2) = IelementsInPatch(idx)
                    IelementsInPatch(idx) = IEL
                    cycle bigsort
                  end if
                  
                end do
                
                ! No other element was found that belongs to the same element
                ! group than IEL. In that case, we begin a new element distribution.
                ilastElementDistr = ilocalElementDistr
                
              end if

            end do bigsort
            
          end do

          ! Reset number of sampling points
          Inpoints = 0

          !-----------------------------------------------------------------------
          ! Step 0:  Allocate temporal memory; note that the dimensions of all
          !          arrays are set to the maximum values required in the process
          !-----------------------------------------------------------------------

          ! Allocate memory for number of corners per element
          allocate(IelementNVEInPatch(nelementsPerBlock))

          ! Allocate memory for number of cubature points per element
          allocate(IelementNcubpInPatch(nelementsPerBlock))
          
          ! Allocate memory for the DOF`s of all the elements
          allocate(IdofsTrial(indofTrialMax,nelementsPerBlock))

          ! Allocate memory for the DOF`s of all the elements
          allocate(IdofsDest(indofDestMax,nelementsPerBlock))

          ! Allocate memory for coefficient values; for mixed triangulations
          ! $ncubpMax$ is an upper bound for the number of cubature points.
          allocate(Dcoefficients(ncubpMax,nelementsPerBlock,p_rtriangulation%ndim))

          ! Calculate the global DOF`s into IdofsTrial.
          call dof_locGlobMapping_mult(p_rdiscrSource, &
              IelementsInPatch(1:nelementsPerBlock), IdofsTrial)
          
          ! Also calculate the global DOF`s in our destination vector(s)
          call dof_locGlobMapping_mult(p_rdiscrDest, &
              IelementsInPatch(1:nelementsPerBlock), IdofsDest)

          ! Allocate memory and get local references to it. This domain integration 
          ! structure stores all information of the source FE space. 
          call domint_initIntegration (rintSubset, nelementsPerBlock, &
              ncubpMax, TRAFO_CS_BARY2DTRI, p_rtriangulation%ndim, NVEMax)
          p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
          p_DcubPtsReal => rintSubset%p_DcubPtsReal
          p_Djac =>        rintSubset%p_Djac
          p_Ddetj =>       rintSubset%p_Ddetj
          p_Dcoords =>     rintSubset%p_DCoords
          
          ! Initialise arrays with zeros
          p_DcubPtsRef  = 0
          p_DcubPtsReal = 0
          p_Djac        = 0
          p_Ddetj       = 0
          p_Dcoords     = 0

          ! Allocate memory. This domain integration structure stores 
          ! all information of the destination FE space. 
          call domint_initIntegration (rintSubsetDest, nelementsPerBlock, &
              nlocalDOFsDestMax, TRAFO_CS_BARY2DTRI, p_rtriangulation%ndim, NVEMax)
          
          ! Since the discretisations are not uniform, we have to treat each
          ! element individually since it may differ from its predecessor.
          ! However, we may skip the re-initialisation of cubature points, etc.
          ! if the current elements belongs to the same element distribution as 
          ! the last one.
          ilastElementDistr = 0

          ! Loop over the patches in the set
          do ipatch = 1, npatchesInCurrentBlock
            
            ! Loop over all elements in the patch.
            idxsubgroup = IelementsInPatchIdx(ipatch)
            
            do while (idxsubgroup .lt. IelementsInPatchIdx(ipatch+1))
            
              ! Get the element distribution of the first element in the patch.
              IEL = IelementsInPatch(idxsubgroup)
              ilocalElementDistr = p_IelementDistr(IEL)
              
              ! Now, find the last element in the current patch which belongs to the
              ! same group of elements; note that we sorted the elements in the patch
              ! before, so all elements of the same element distribution are behind
              ! each other. 
              ! All these elements have the same basis functions, the same transformation,
              ! coordinate system,...
              !
              ! Start at the end of the patch and find the last element within the
              ! same element distribution as the first one. The loop ends
              ! - if such an element is found or
              ! - if no element is found; then (by Fortrag standard) idx2 points
              !   also to the first element of the patch -- our 'reference' element
              !   of the subgroup.
              ! It is better to do the loop 'from end to the beginning' instead of
              ! 'from start to the end' as there may be many patches containing only
              ! elements of the same kind; in that case, the loop is immediately left
              ! as the last element has the same kind as the first one.
              do idx2 = IelementsInPatchIdx(ipatch+1)-1, idxsubgroup+1, -1
                if (ilocalElementDistr .eq. p_IelementDistr(IelementsInPatch(idx2))) exit
              end do
              
              ! Ok, the elements of 'the same kind' can now be found between
              ! positions idxsubgroup..idx2.
              !
              ! In case the 'local' element distribution changed, activate the new one.
              ! This small check saves some time if there are some patches with completely
              ! the same element type.
              
              if (ilocalElementDistr .ne. ilastElementDistr) then

                ! Active the local element distribution of that group.
                p_relementDistribution => p_rdiscrSource%RelementDistr(ilocalElementDistr)

                ! Get the number of local DOF`s for trial functions
                indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

                ! Get the number of corner vertices of the element
                NVE = elem_igetNVE(p_relementDistribution%celement)

                !---------------------------------------------------------------------
                ! Step 1:  Prepare the source FE space
                !---------------------------------------------------------------------
                
                ! Initialise the cubature formula for the local element distribution
                ! Get cubature weights and point coordinates on the reference element
                call cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)

                ! Get from the trial element space the type of coordinate system
                ! that is used in the local element distribution
                icoordSystem = elem_igetCoordSystem(p_relementDistribution%celement)
                                
                ! Check if one of the trial/test elements is nonparametric
                bnonparTrial = elem_isNonparametric(p_relementDistribution%celement)
                
                ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                ! p_DcubPtsRef - depending on whether the space is parametric or not.
                if (bnonparTrial) then
                  p_DcubPtsTrial => p_DcubPtsReal
                else
                  p_DcubPtsTrial => p_DcubPtsRef
                end if
                
              end if
              
              ! Update number of sampling points for current patch
              Inpoints(ipatch) = Inpoints(ipatch) + ncubp*(idx2-idxsubgroup+1)

              ! Loop over elements in subgroup
              do idx = idxsubgroup,idx2

                ! Put the cubature point coordinates in the right format to the 
                ! cubature-point array.
                ! Initialise all entries in p_DcubPtsRef with the same coordinates -
                ! as the cubature point coordinates are identical on all elements.
                ! Importantly, only those entries are copied which correspond to
                ! the local element type. The remaining entries are left equal to zero.
                do i=1, ncubp
                  do k=1,size(p_DcubPtsRef,1)
                    ! Could be solved using the TRANSPOSE operator - but often it is
                    ! faster this way...
                    p_DcubPtsRef(k,i,idx) = Dxi(i,k)
                  end do
                end do

              end do
              
              ! Store number of corners per element
              IelementNVEInPatch(idxsubgroup:idx2) = NVE

              ! Store number of cubature points per element
              IelementNcubpInPatch(idxsubgroup:idx2) = ncubp

              ! We have the coordinates of the cubature points saved in the
              ! coordinate array from above. Unfortunately for nonparametric
              ! elements, we need the real coordinate.
              ! Furthermore, we anyway need the coordinates of the element
              ! corners and the Jacobian determinants corresponding to
              ! all the points.
              
              ! At first, get the coordinates of the corners of the element.
              call trafo_getCoords_sim(&
                  elem_igetTrafoType(p_relementDistribution%celement), &
                  p_rtriangulation, IelementsInPatch(idxsubgroup:idx2), &
                  p_Dcoords(:,:,idxsubgroup:idx2))

              ! Depending on the type of transformation, we must now choose
              ! the mapping between the reference and the real element.
              ! In case we use a nonparametric element as test function, we need the 
              ! coordinates of the points on the real element, too.
              ! Unfortunately, we need the real coordinates of the cubature points
              ! anyway for the function - so calculate them all.
              call trafo_calctrafo_sim (p_relementDistribution%ctrafoType, &
                  idx2-idxsubgroup+1,ncubp, &
                  p_Dcoords(:,:,idxsubgroup:idx2), p_DcubPtsRef(:,:,idxsubgroup:idx2), &
                  p_Djac(:,:,idxsubgroup:idx2), &
                  p_Ddetj(:,idxsubgroup:idx2), p_DcubPtsReal(:,:,idxsubgroup:idx2))

              !---------------------------------------------------------------------
              ! Step 2:  Perform sampling of consistent gradient values
              !---------------------------------------------------------------------
              
              ! Calculate the derivative of the FE function in the cubature
              ! points: u_h(x,y,z) and save the result to Dcoefficients(:,:,1..3)
              select case(p_rtriangulation%ndim)
              case (NDIM1D)
                call fevl_evaluate_sim (rvectorScalar, p_Dcoords(:,:,idxsubgroup:idx2), &
                    p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2), &
                    p_relementDistribution%celement, IdofsTrial(:,idxsubgroup:idx2), &
                    ncubp, idx2-idxsubgroup+1,&
                    p_DcubPtsTrial(:,:,idxsubgroup:idx2), DER_DERIV1D_X, &
                    Dcoefficients(:,idxsubgroup:idx2,1))
                
              case (NDIM2D)
                call fevl_evaluate_sim (rvectorScalar, p_Dcoords(:,:,idxsubgroup:idx2), &
                    p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2), &
                    p_relementDistribution%celement, IdofsTrial(:,idxsubgroup:idx2), &
                    ncubp, idx2-idxsubgroup+1,&
                    p_DcubPtsTrial(:,:,idxsubgroup:idx2), DER_DERIV2D_X, &
                    Dcoefficients(:,idxsubgroup:idx2,1))

                call fevl_evaluate_sim (rvectorScalar, p_Dcoords(:,:,idxsubgroup:idx2), &
                    p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2), &
                    p_relementDistribution%celement, IdofsTrial(:,idxsubgroup:idx2), &
                    ncubp, idx2-idxsubgroup+1,&
                    p_DcubPtsTrial(:,:,idxsubgroup:idx2), DER_DERIV2D_Y, &
                    Dcoefficients(:,idxsubgroup:idx2,2))
                
              case (NDIM3D)
                call fevl_evaluate_sim (rvectorScalar, p_Dcoords(:,:,idxsubgroup:idx2), &
                    p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2), &
                    p_relementDistribution%celement, IdofsTrial(:,idxsubgroup:idx2), &
                    ncubp, idx2-idxsubgroup+1,&
                    p_DcubPtsTrial(:,:,idxsubgroup:idx2), DER_DERIV3D_X, &
                    Dcoefficients(:,idxsubgroup:idx2,1))
                
                call fevl_evaluate_sim (rvectorScalar, p_Dcoords(:,:,idxsubgroup:idx2), &
                    p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2), &
                    p_relementDistribution%celement, IdofsTrial(:,idxsubgroup:idx2), &
                    ncubp, idx2-idxsubgroup+1,&
                    p_DcubPtsTrial(:,:,idxsubgroup:idx2), DER_DERIV3D_X, &
                    Dcoefficients(:,idxsubgroup:idx2,2))
                
                call fevl_evaluate_sim (rvectorScalar, p_Dcoords(:,:,idxsubgroup:idx2), &
                    p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2), &
                    p_relementDistribution%celement, IdofsTrial(:,idxsubgroup:idx2), &
                    ncubp, idx2-idxsubgroup+1,&
                    p_DcubPtsTrial(:,:,idxsubgroup:idx2), DER_DERIV3D_Z, &
                    Dcoefficients(:,idxsubgroup:idx2,3))
                
              case DEFAULT
                call output_line('Invalid spatial dimension!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
                call sys_halt()
              end select             

              ! Subgroup finished, continue with the next one.
              idxsubgroup = idx2 + 1
              
            end do
              
          end do   ! End of IPATCH loop
          
          !---------------------------------------------------------------------
          ! Step 3: Prepare least-squares fitting
          !---------------------------------------------------------------------

          ! First, get the coordinates of the bounding group for the set of
          ! elements present in each patch. In addition, calculate the coordinates
          ! of the corners of the constant Jacobian patch "elements".
          call calc_patchBoundingGroup_mult(IelementsInPatchIdx, npatchesInCurrentBlock, &
              IelementNVEInPatch, p_Dcoords, DpatchBound(:,:,1:npatchesInCurrentBlock))

          ! Reactive current element distribution of the source FE space. This is the 
          ! element distribution of those elements which make up the center of each 
          ! patch. Hence, this element distribution determine the type of element 
          ! used to construct the constant Jacobian "patch" elements.
          p_relementDistribution => p_rdiscrSource%RelementDistr(icurrentElementDistr)

          ! Initialise the cubature formula for the local element distribution
          ! Get cubature weights and point coordinates on the reference element
          call cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
          
          ! Get the number of local DOF`s for trial functions
          indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

          ! Get the number of corner vertices of the element
          NVE = elem_igetNVE(p_relementDistribution%celement)
          
          ! Get from the trial element space the type of coordinate system
          ! that is used there:
          icoordSystem = elem_igetCoordSystem(p_relementDistribution%celement)

          ! Check if one of the trial/test elements is nonparametric
          bnonparTrial = elem_isNonparametric(p_relementDistribution%celement)
          
          ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
          ! p_DcubPtsRef - depending on whether the space is parametric or not.
          if (bnonparTrial) then
            p_DcubPtsTrial => p_DcubPtsReal
          else
            p_DcubPtsTrial => p_DcubPtsRef
          end if
          
          ! Next, we need to convert the physical coordinates of the curvature points
          ! to the local coordinates of the constant Jacobian "patch" elements
          call calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem,&
              npatchesInCurrentBlock, NVE, p_DcubPtsReal, DpatchBound, p_DcubPtsRef)

          ! Depending on the type of transformation, we must now choose
          ! the mapping between the reference and the real element.
          ! In case we use a nonparametric element as test function, we need the 
          ! coordinates of the points on the real element, too.
          ! Unfortunately, we need the real coordinates of the cubature points
          ! anyway for the function - so calculate them all.
          call trafo_calctrafo_sim (p_relementDistribution%ctrafoType, &
              nelementsPerBlock, ncubpMax, p_Dcoords, p_DcubPtsRef, p_Djac, p_Ddetj)

          ! Allocate memory for the patch interpolants matrices
          ! Note that the second dimension is DER_FUNC=1 and could be omitted. However,
          ! all routines for simultaneous element evaluation require a 4d-array so that
          ! the array Dpolynomials is artificially created as 4d-array.
          allocate(Dpolynomials(indofTrial,DER_FUNC,ncubpMax,nelementsPerBlock))

          ! Evaluate the trial functions of the constant Jacobian patch "element" for all
          ! cubature points of the elements present in the patch and store each polynomial 
          ! interpolation in the rectangular patch matrix used for least-squares fitting.
          ! Note that we evaluate over the maximum number of cubature points present in
          ! the patch. Although some meaningless values may be generated, it is faster to
          ! evaluate all values simultaneously and filter the required data afterwards.
          call elem_generic_sim(p_relementDistribution%celement,&
              p_Dcoords, p_Djac, p_Ddetj, BderDest, Dpolynomials, ncubpMax,&
              nelementsPerBlock, p_DcubPtsTrial)

          ! Compute total number of sampling points
          nspoints = sum(Inpoints(1:npatchesInCurrentBlock))

          ! Allocate memory for the aligned patch data
          allocate(DcoefficientsMixed(nspoints,p_rtriangulation%ndim))
          allocate(DpolynomialsMixed(indofTrial*nspoints))
          
          ! Reset total number of sampling/interpolation points
          nspoints = 0; nipoints = 0

          ! Loop over all element indices in the set of patches
          ! and align the polynomials interpolants
          do idx = 1, IelementsInPatchIdx(npatchesInCurrentBlock+1)-1
            
            ! Get number of cubature points for current element
            ncubp = IelementNcubpInPatch(idx)

            ! Update total number of sampling points
            nspoints = nspoints+ncubp

            ! Apply coefficients for each dimension
            do idim = 1, p_rtriangulation%ndim
              DcoefficientsMixed(nspoints-ncubp+1:nspoints,idim) =&
                  Dcoefficients(1:ncubp,idx,idim)
            end do

            ! Loop over all cubature points
            do icubp = 1, ncubp
              
              ! Update total number of interpolation points
              nipoints = nipoints+indofTrial
              
              ! Align polynomial interpolats
              DpolynomialsMixed(nipoints-indofTrial+1:nipoints) =&
                   Dpolynomials(1:indofTrial,DER_FUNC,icubp,idx)
            end do
          end do

          ! Deallocate memory for the unaligned data
          deallocate(Dpolynomials, Dcoefficients)


          !-----------------------------------------------------------------------
          ! Step 4: Perform least-squares fitting
          !-----------------------------------------------------------------------

          ! Allocate memory for the derivative values
          allocate(Dderivatives(indofTrial,npatchesInCurrentBlock,p_rtriangulation%ndim))
          
          ! Compute the patch averages by solving $(P^T * P) * x = (P^T) * b$ for x
          call calc_patchAverages_mult(IelementsInPatchIdx, npatchesInCurrentBlock, &
              Inpoints, indofTrial, DcoefficientsMixed, DpolynomialsMixed, Dderivatives)

          deallocate(DpolynomialsMixed, DcoefficientsMixed)


          !-----------------------------------------------------------------------
          ! Step 5: Prepare the destination FE space
          !-----------------------------------------------------------------------

          ! Get local references to domain integration structure
          p_DcubPtsRef =>  rintSubsetDest%p_DcubPtsRef
          p_DcubPtsReal => rintSubsetDest%p_DcubPtsReal
          p_Djac =>        rintSubsetDest%p_Djac
          p_Ddetj =>       rintSubsetDest%p_Ddetj
          p_Dcoords =>     rintSubsetDest%p_DCoords
          
          ! Initialise arrays with zeros
          p_DcubPtsRef  = 0
          p_DcubPtsReal = 0
          p_Djac        = 0
          p_Ddetj       = 0
          p_Dcoords     = 0

          ! Loop over the patches in the set

          ilastElementDistr = 0

          do ipatch = 1, npatchesInCurrentBlock
            
            ! Loop over all elements in the patch.
            idxsubgroup = IelementsInPatchIdx(ipatch)
            
            do while (idxsubgroup .lt. IelementsInPatchIdx(ipatch+1))
            
              ! Get the element distribution of the first element in the patch.
              IEL = IelementsInPatch(idxsubgroup)
              ilocalElementDistr = p_IelementDistr(IEL)
              
              ! Now, find the last element in the current patch which belongs to the
              ! same group of elements; note that we sorted the elements in the patch
              ! before, so all elements of the same element distribution are behind
              ! each other. 
              do idx2 = IelementsInPatchIdx(ipatch+1)-1, idxsubgroup+1, -1
                if (ilocalElementDistr .eq. p_IelementDistr(IelementsInPatch(idx2))) exit
              end do
              
              ! Ok, the elements of 'the same kind' can now be found between
              ! positions idxsubgroup..idx2.
              !
              ! In case the 'local' element distribution changed, activate the new one.
              ! This small check saves some time if there are some patches with completely
              ! the same element type.
              if (ilocalElementDistr .ne. ilastElementDistr) then

                ! Active local element distribution
                p_relementDistrDest => p_rdiscrDest%RelementDistr(ilocalElementDistr)

                ! Get the number of local DOF`s for trial functions
                indofDest = elem_igetNDofLoc(p_relementDistrDest%celement)

                ! Initialise arrays with zeros
                Dxi    = 0
                Domega = 0

                ! Initialise the cubature formula. That is a special trick here!
                ! In particular, we only need the physical coordinates of the
                ! nodal evaluation points in the source vectors.
                ! For this purpose, we initialise the 'trapezoidal rule' as cubature
                ! formula. Ok, we are not using cubature at all here (thus ignoring
                ! the cubature weights completely), but that way we get the
                ! coordinates of the corners on the reference element automatically --
                ! which coincide with the points where we want to create the gradients!
                !
                ! Note: The returned nlocalDOFsDest will coincide with the number of local DOF`s
                ! on each element indofDest!
                call calc_cubatureDest(&
                    elem_getPrimaryElement(p_relementDistrDest%celement), &
                    nlocalDOFsDest, Dxi, Domega)
                                
                ! Check if one of the trial/test elements is nonparametric
                bnonparTrial = elem_isNonparametric(p_relementDistrDest%celement)
                
                ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                ! p_DcubPtsRef - depending on whether the space is parametric or not.
                if (bnonparTrial) then
                  p_DcubPtsTrial => p_DcubPtsReal
                else
                  p_DcubPtsTrial => p_DcubPtsRef
                end if
                
                ! Save number of last element distribution
                ilastElementDistr = ilocalElementDistr
              end if

              ! Put the cubature point coordinates in the right format to the
              ! cubature-point array.
              ! Initialise all entries in p_DcubPtsRef with the same coordinates -
              ! as the cubature point coordinates are identical on all elements
              do idx = idxsubgroup,idx2
                do i=1,nlocalDOFsDest
                  do k=1,size(p_DcubPtsRef,1)
                    ! Could be solved using the TRANSPOSE operator - but often it is 
                    ! faster this way...
                    p_DcubPtsRef(k,i,idx) = Dxi(i,k)
                  end do
                end do
              end do
              
              ! We have the coordinates of the cubature points saved in the
              ! coordinate array from above. Unfortunately for nonparametric
              ! elements, we need the real coordinate.
              ! Furthermore, we anyway need the coordinates of the element
              ! corners and the Jacobian determinants corresponding to
              ! all the points.
              
              ! At first, get the coordinates of the corners of the element
              call trafo_getCoords_sim (&
                  elem_igetTrafoType(p_relementDistrDest%celement), &
                  p_rtriangulation, IelementsInPatch(idxsubgroup:idx2) , &
                  p_Dcoords(:,:,idxsubgroup:idx2))
              
              ! Depending on the type of transformation, we must now choose
              ! the mapping between the reference and the real element.
              ! In case we use a nonparametric element as test function, we need the 
              ! coordinates of the points on the real element, too.
              ! Unfortunately, we need the real coordinates of the cubature points
              ! anyway for the function - so calculate them all.
              call trafo_calctrafo_sim (p_relementDistrDest%ctrafoType, &
                  idx2-idxsubgroup+1,nlocalDOFsDest, &
                  p_Dcoords(:,:,idxsubgroup:idx2),p_DcubPtsRef(:,:,idxsubgroup:idx2),&
                  p_Djac(:,:,idxsubgroup:idx2), &
                  p_Ddetj(:,idxsubgroup:idx2), p_DcubPtsReal(:,:,idxsubgroup:idx2))
                  
              ! Subgroup finished, continue with the next one.
              idxsubgroup = idx2 + 1
                  
            end do
          end do
          
          ! Next, we need to convert the physical coordinates of the curbature points
          ! to the local coordinates of the constant Jacobian "patch" elements. Note that 
          ! NVE and icoordSystem have not been modified and can be savely used from above.
          call calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, npatchesInCurrentBlock, &
              NVE,  p_DcubPtsReal, DpatchBound, p_DcubPtsRef)

          ! We do not need the corner coordinates of the elements in the destination
          ! FE space but that of the "patch" elements in the source FE space
          p_Dcoords => rintSubset%p_Dcoords


          ! It is mandatory that the second dimensions is DER_MAXNDER even though only the
          ! function values (DER_FUNC) are actually needed. During the evaluation of the basis
          ! functions some elements do not check if first derivatives are required. In short,
          ! IF-THEN-ELSE is more costly than just computing some superficial values ;-)
          ! Therefore, we must provide sufficient memory !!!
          allocate(Dpolynomials(indofTrial,DER_MAXNDER,nlocalDOFsDestMax,nelementsPerBlock))
          Dpolynomials = 0.0_DP

          ! Loop over the patches in the set
          ilastElementDistr = 0
          
          do ipatch = 1, npatchesInCurrentBlock
            
            ! Loop over all elements in the patch.
            idxsubgroup = IelementsInPatchIdx(ipatch)
            
            do while (idxsubgroup .lt. IelementsInPatchIdx(ipatch+1))
            
              ! Get the element distribution of the first element in the patch.
              IEL = IelementsInPatch(idxsubgroup)
              ilocalElementDistr = p_IelementDistr(IEL)
              
              ! Now, find the last element in the current patch which belongs to the
              ! same group of elements; note that we sorted the elements in the patch
              ! before, so all elements of the same element distribution are behind
              ! each other. 
              do idx2 = IelementsInPatchIdx(ipatch+1)-1, idxsubgroup+1, -1
                if (ilocalElementDistr .eq. p_IelementDistr(IelementsInPatch(idx2))) exit
              end do
              
              ! Ok, the elements of 'the same kind' can now be found between
              ! positions idxsubgroup..idx2.
              !
              ! In case the 'local' element distribution changed, activate the new one.
              ! This small check saves some time if there are some patches with completely
              ! the same element type.
              if (ilocalElementDistr .ne. ilastElementDistr) then

                ! Active local element distribution
                p_relementDistrDest => p_rdiscrDest%RelementDistr(ilocalElementDistr)
                
                ! Get the number of local DOF`s for trial functions which coincides with
                nlocalDOFsDest = elem_igetNDofLoc(p_relementDistrDest%celement)

                ! Check if one of the trial/test elements is nonparametric
                bnonparTrial = elem_isNonparametric(p_relementDistrDest%celement)
                
                ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                ! p_DcubPtsRef - depending on whether the space is parametric or not.
                if (bnonparTrial) then
                  p_DcubPtsTrial => p_DcubPtsReal
                else
                  p_DcubPtsTrial => p_DcubPtsRef
                end if

                ! Save number of last element distribution
                ilastElementDistr = ilocalElementDistr
              end if
                
              ! Calculate the transformation from the reference element to the real one
              call trafo_calctrafo_sim(p_relementDistrDest%ctrafoType, &
                  idx2-idxsubgroup+1,nlocalDOFsDest, &
                  p_Dcoords(:,:,idxsubgroup:idx2), p_DcubPtsRef(:,:,idxsubgroup:idx2),&
                  p_Djac(:,:,idxsubgroup:idx2), p_Ddetj(:,idxsubgroup:idx2))

              ! Evaluate the basis functions for the cubature points of the destination FE space
              call elem_generic_sim(p_relementDistribution%celement, &
                  p_Dcoords(:,:,idxsubgroup:idx2), p_Djac(:,:,idxsubgroup:idx2), &
                  p_Ddetj(:,idxsubgroup:idx2), &
                  BderDest, Dpolynomials(:,:,:,idxsubgroup:idx2), &
                  nlocalDOFsDest,idx2-idxsubgroup+1, &
                  p_DcubPtsTrial(:,:,idxsubgroup:idx2))
                  
              ! Subgroup finished, continue with the next one.
              idxsubgroup = idx2 + 1
                  
            end do
          end do
          ilastElementDistr = 0

          !---------------------------------------------------------------------
          ! Step 6: Evaluate the averaged derivative values at the cubature
          !         points of the destination FE space and scatter them to
          !         the global degrees of freedom
          !---------------------------------------------------------------------

          select case (p_rtriangulation%ndim)
          case(NDIM1D)
            
            ! Loop over the patches in the set
            do ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              do idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1

                ! Get global element number
                IEL = IelementsInPatch(idx)
                
                ! Get number of local element distribution
                ilocalElementDistr = p_IelementDistr(IEL)
                
                ! Check if local element distribution corresponds to the last element 
                ! distribution. Then we do not have to initialise everything again.
                if (ilocalElementDistr .ne. ilastElementDistr) then
                  
                  ! Active local element distribution
                  p_relementDistrDest => p_rdiscrDest%RelementDistr(ilocalElementDistr)
                  
                  ! Get the number of local DOF`s for trial functions which coincides with
                  nlocalDOFsDest = elem_igetNDofLoc(p_relementDistrDest%celement)

                  ! Check if one of the trial/test elements is nonparametric
                  bnonparTrial = elem_isNonparametric(p_relementDistrDest%celement)
                  
                  ! Let p_DcubPtsTrial point either to p_DcubPtsReal or
                  ! p_DcubPtsRef - depending on whether the space is parametric or not.
                  if (bnonparTrial) then
                    p_DcubPtsTrial => p_DcubPtsReal
                  else
                    p_DcubPtsTrial => p_DcubPtsRef
                  end if

                  ! Save number of last element distribution
                  ilastElementDistr = ilocalElementDistr
                end if

                ! Loop over local degrees of freedom
                do ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  do j = 1, indofTrial
                    Dval(1) = Dval(1) + Dderivatives(j,ipatch,1) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  end do
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                end do
              end do
            end do
            ilastElementDistr = 0

          case(NDIM2D)
            
            ! Loop over the patches in the set
            do ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              do idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1

                ! Get global element number
                IEL = IelementsInPatch(idx)
                
                ! Get number of local element distribution
                ilocalElementDistr = p_IelementDistr(IEL)
                
                ! Check if local element distribution corresponds to the last element 
                ! distribution. Then we do not have to initialise everything again.
                if (ilocalElementDistr .ne. ilastElementDistr) then
                  
                  ! Active local element distribution
                  p_relementDistrDest => p_rdiscrDest%RelementDistr(ilocalElementDistr)
                  
                  ! Get the number of local DOF`s for trial functions which coincides with
                  nlocalDOFsDest = elem_igetNDofLoc(p_relementDistrDest%celement)
                  
                  ! Save number of last element distribution
                  ilastElementDistr = ilocalElementDistr
                end if

                ! Loop over local degrees of freedom
                do ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  do j = 1, indofTrial
                    Dval(1:2) = Dval(1:2) + Dderivatives(j,ipatch,1:2) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  end do
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                end do
              end do
            end do
            ilastElementDistr = 0

          case(NDIM3D)
            
            ! Loop over the patches in the set
            do ipatch = 1, npatchesInCurrentBlock
              
              ! Loop over elements in patch
              do idx = IelementsInPatchIdx(ipatch),IelementsInPatchIdx(ipatch+1)-1

                ! Get global element number
                IEL = IelementsInPatch(idx)
                
                ! Get number of local element distribution
                ilocalElementDistr = p_IelementDistr(IEL)
                
                ! Check if local element distribution corresponds to the last element 
                ! distribution. Then we do not have to initialise everything again.
                if (ilocalElementDistr .ne. ilastElementDistr) then
                  
                  ! Active local element distribution
                  p_relementDistrDest => p_rdiscrDest%RelementDistr(ilocalElementDistr)
                  
                  ! Get the number of local DOF`s for trial functions which coincides with
                  nlocalDOFsDest = elem_igetNDofLoc(p_relementDistrDest%celement)
                  
                  ! Save number of last element distribution
                  ilastElementDistr = ilocalElementDistr
                end if

                ! Loop over local degrees of freedom
                do ipoint= 1, nlocalDOFsDest
                  
                  Dval = 0.0_DP
                  do j = 1, indofTrial
                    Dval(1:3) = Dval(1:3) + Dderivatives(j,ipatch,1:3) * Dpolynomials(j,DER_FUNC,ipoint,idx)
                  end do
                  
                  ! Scatter to global degrees of freedom
                  p_DxDeriv(IdofsDest(ipoint,idx)) = p_DxDeriv(IdofsDest(ipoint,idx)) + Dval(1)
                  p_DyDeriv(IdofsDest(ipoint,idx)) = p_DyDeriv(IdofsDest(ipoint,idx)) + Dval(2)
                  p_DzDeriv(IdofsDest(ipoint,idx)) = p_DzDeriv(IdofsDest(ipoint,idx)) + Dval(3)
                  
                  ! Count how often a DOF was touched.
                  p_IcontributionsAtDOF(IdofsDest(ipoint,idx)) = &
                      p_IcontributionsAtDOF(IdofsDest(ipoint,idx))+1
                end do
              end do
            end do
            ilastElementDistr = 0

          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradSuperPatchRecov')
            call sys_halt()
          end select

          ! Release memory
          call domint_doneIntegration(rintSubset)
          call domint_doneIntegration(rintSubsetDest)

          ! Deallocate temporary memory
          deallocate(Dderivatives)
          deallocate(Dpolynomials)
          deallocate(IdofsDest)
          deallocate(IdofsTrial)
          deallocate(IelementNVEInPatch)
          deallocate(IelementNcubpInPatch)
        end if
        
        ! Deallocate temporary memory
        deallocate(IelementsInPatch)

      end do   ! End of PATCHset loop
      
      ! Deallocate temporary memory
      if (.not.bisUniform) deallocate(Inpoints)
      deallocate(DpatchBound)
      deallocate(IelementsInPatchIdx)

    end do   ! End of icurrentElementDistr loop

    ! We are nearly done. The final thing: divide the calculated derivatives by the
    ! number of elements adjacent to each vertex. That closes the calculation
    ! of the 'mean' of the derivatives.
    select case(p_rtriangulation%ndim)
    case (NDIM1D)
      do i=1,size(p_DxDeriv)
        ! Div/0 should not occur, otherwise the triangulation is crap
        ! as there is a point not connected to any element!
        p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
      end do

    case (NDIM2D)
      do i=1,size(p_DxDeriv)
        ! Div/0 should not occur, otherwise the triangulation is crap
        ! as there is a point not connected to any element!
        p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
        p_DyDeriv(i) = p_DyDeriv(i) / p_IcontributionsAtDOF(i)
      end do

    case (NDIM3D)
      do i=1,size(p_DxDeriv)
        ! Div/0 should not occur, otherwise the triangulation is crap
        ! as there is a point not connected to any element!
        p_DxDeriv(i) = p_DxDeriv(i) / p_IcontributionsAtDOF(i)
        p_DyDeriv(i) = p_DyDeriv(i) / p_IcontributionsAtDOF(i)
        p_DzDeriv(i) = p_DzDeriv(i) / p_IcontributionsAtDOF(i)
      end do
    end select

    ! Free temporary memory
    call storage_free(h_IcontributionsAtDOF)
    
  contains
    
    ! Here, some auxiliary working routines follow

    !**************************************************************
    ! Calculate the "bounding group" for a set of patches
    !
    ! Each patch consists of multiple elements which are adjacent to 
    ! each other so that they cover some connected part of the domain.
    ! This routine determines the physical coordinates of the constant
    ! Jacobian patch "element" that has its local axes parallel to the
    ! global axes and completely surrounds the elements of the patch.
    ! Moreover, this routine converts the physical coordinates of the
    ! evaluation points in the destination space to the local coordinates
    ! in the constant Jacobian patch "element".

    subroutine calc_patchBoundingGroup_mult(IelementsInPatchIdx, &
        npatches, IelemNVE, DpointsReal, DpointsBound)
      
      ! Index vector and list of elements in patch
      integer, dimension(:), intent(in) :: IelementsInPatchIdx

      ! Number of elements in patch
      integer, intent(in) :: npatches

      ! Array with numbers of corners per element
      integer, dimension(:), intent(in) :: IelemNVE

      ! Physical coordinates of the corner nodes of all elements in the patch
      real(DP), dimension(:,:,:), intent(inout) :: DpointsReal

      ! Physical coordinates of the bounds of all "patch" elements
      real(DP), dimension(:,:,:), intent(out) :: DpointsBound

      ! local variables
      integer :: ipatch,idx,idxFirst,idxLast
      real(DP) :: xmin,ymin,xmax,ymax,zmin,zmax

      select case (size(DpointsReal,1))

      case (NDIM1D)
        
        ! Loop over all patches
        do ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1
          
          ! Determine minimum/maximum coordinates of the patch
          xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
          xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
          
          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(1,2,ipatch) = xmax
          
          ! Store physical coordinates of patch element corners
          DpointsReal(1,1,idxFirst:idxLast) = xmin
          DpointsReal(1,2,idxFirst:idxLast) = xmax
        end do
        
      case (NDIM2D)

        ! Loop over all patches
        do ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1
          
          ! Determine minimum/maximum coordinates of the patch
          xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
          xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
          ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
          ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
          
          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(2,1,ipatch) = ymin
          DpointsBound(1,2,ipatch) = xmax
          DpointsBound(2,2,ipatch) = ymax
          
          ! Loop over all elements
          do idx = idxFirst, idxLast
            select case (IelemNVE(idx))
            case (TRIA_NVETRI2D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(1,2,idx) = xmax+ymax-ymin
              DpointsReal(2,2,idx) = ymin
              DpointsReal(1,3,idx) = xmin
              DpointsReal(2,3,idx) = xmax+ymax-xmin

            case (TRIA_NVEQUAD2D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(1,2,idx) = xmax
              DpointsReal(2,2,idx) = ymin
              DpointsReal(1,3,idx) = xmax
              DpointsReal(2,3,idx) = ymax
              DpointsReal(1,4,idx) = xmin
              DpointsReal(2,4,idx) = ymax

            case DEFAULT
              call output_line ('Invalid number of vertices per elements!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_mult')
              call sys_halt()
            end select
          end do
        end do
        
      case (NDIM3D)

        ! Loop over all patches
        do ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Determine minimum/maximum coordinates of the patch
          xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
          xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
          ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
          ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
          zmin = minval(DpointsReal(3,:,idxFirst:idxLast))
          zmax = maxval(DpointsReal(3,:,idxFirst:idxLast))

          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(2,1,ipatch) = ymin
          DpointsBound(3,1,ipatch) = zmin
          DpointsBound(1,2,ipatch) = xmax
          DpointsBound(2,2,ipatch) = ymax
          DpointsBound(3,2,ipatch) = zmax

          ! Loop over all elements
          do idx = idxFirst, idxLast
            select case (IelemNVE(idx))
            case (TRIA_NVETET3D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(3,1,idx) = zmin
              DpointsReal(1,2,idx) = xmax+ymax-zmax-ymin+zmin
              DpointsReal(2,2,idx) = ymin
              DpointsReal(3,2,idx) = zmin
              DpointsReal(1,3,idx) = xmin
              DpointsReal(2,3,idx) = xmax+ymax-zmax-xmin+zmin
              DpointsReal(3,3,idx) = zmin
              DpointsReal(1,4,idx) = xmin
              DpointsReal(2,4,idx) = ymin
              DpointsReal(3,4,idx) = -xmax-ymax+zmax+xmin+ymin

            case (TRIA_NVEPYR3D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(3,1,idx) = zmin
              DpointsReal(1,2,idx) = xmax+zmax-zmin
              DpointsReal(2,2,idx) = ymin
              DpointsReal(3,2,idx) = zmin
              DpointsReal(1,3,idx) = xmax+zmax-zmin
              DpointsReal(2,3,idx) = ymax+zmax-zmin
              DpointsReal(3,3,idx) = zmin
              DpointsReal(1,4,idx) = xmin
              DpointsReal(2,4,idx) = ymax+zmax-zmin
              DpointsReal(3,4,idx) = zmin
              DpointsReal(1,5,idx) = xmin
              DpointsReal(2,5,idx) = ymin
              DpointsReal(3,5,idx) = xmax+zmax-xmin

            case (TRIA_NVEPRIS3D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(3,1,idx) = zmin
              DpointsReal(1,2,idx) = xmax+ymax-ymin
              DpointsReal(2,2,idx) = ymin
              DpointsReal(3,2,idx) = zmin
              DpointsReal(1,3,idx) = xmin
              DpointsReal(2,3,idx) = xmax+ymax-xmin
              DpointsReal(3,3,idx) = zmin
              DpointsReal(1,4,idx) = xmin
              DpointsReal(2,4,idx) = ymin
              DpointsReal(3,4,idx) = zmax
              DpointsReal(1,5,idx) = xmax+ymax-ymin
              DpointsReal(2,5,idx) = ymin
              DpointsReal(3,5,idx) = zmax
              DpointsReal(1,6,idx) = xmin
              DpointsReal(2,6,idx) = xmax+ymax-xmin
              DpointsReal(3,6,idx) = zmax
              
            case (TRIA_NVEHEXA3D)
              ! Store physical coordinates of patch element corners
              DpointsReal(1,1,idx) = xmin
              DpointsReal(2,1,idx) = ymin
              DpointsReal(3,1,idx) = zmin
              DpointsReal(1,2,idx) = xmax
              DpointsReal(2,2,idx) = ymin
              DpointsReal(3,2,idx) = zmin
              DpointsReal(1,3,idx) = xmax
              DpointsReal(2,3,idx) = ymax
              DpointsReal(3,3,idx) = zmin
              DpointsReal(1,4,idx) = xmin
              DpointsReal(2,4,idx) = ymax
              DpointsReal(3,4,idx) = zmin
              DpointsReal(1,5,idx) = xmin
              DpointsReal(2,5,idx) = ymin
              DpointsReal(3,5,idx) = zmax
              DpointsReal(1,6,idx) = xmax
              DpointsReal(2,6,idx) = ymin
              DpointsReal(3,6,idx) = zmax
              DpointsReal(1,7,idx) = xmax
              DpointsReal(2,7,idx) = ymax
              DpointsReal(3,7,idx) = zmax
              DpointsReal(1,8,idx) = xmin
              DpointsReal(2,8,idx) = ymax
              DpointsReal(3,8,idx) = zmax

            case DEFAULT
              call output_line ('Invalid number of vertices per elements!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_mult')
              call sys_halt()
            end select
          end do
        end do
                    
      case DEFAULT
        call output_line ('Invalid number of spatial dimensions!', &
            OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_mult')
        call sys_halt()
      end select
    end subroutine calc_patchBoundingGroup_mult
    
    !**************************************************************
    ! Calculate the "bounding group" for a set of patches
    !
    ! Each patch consists of multiple elements which are adjacent to 
    ! each other so that they cover some connected part of the domain.
    ! This routine determines the physical coordinates of the constant
    ! Jacobian patch "element" that has its local axes parallel to the
    ! global axes and completely surrounds the elements of the patch.
    ! Moreover, this routine converts the physical coordinates of the
    ! evaluation points in the destination space to the local coordinates
    ! in the constant Jacobian patch "element".

    subroutine calc_patchBoundingGroup_sim(IelementsInPatchIdx, &
        npatches, NVE, DpointsReal, DpointsBound)
      
      ! Index vector and list of elements in patch
      integer, dimension(:), intent(in) :: IelementsInPatchIdx

      ! Number of elements in patch
      integer, intent(in) :: npatches

      ! Number of vertices per element
      integer, intent(in) :: NVE

      ! Physical coordinates of the corner nodes of all elements in the patch
      real(DP), dimension(:,:,:), intent(inout) :: DpointsReal

      ! Physical coordinates of the bounds of all "patch" elements
      real(DP), dimension(:,:,:), intent(out) :: DpointsBound

      ! local variables
      integer :: ipatch,idxFirst,idxLast
      real(DP) :: xmin,ymin,xmax,ymax,zmin,zmax

      select case (size(DpointsReal,1))
        
      case (NDIM1D)

        ! Simple: Just find minimal/maximal value
        
        ! Loop over all patches
        do ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Determine minimum/maximum coordinates of the patch
          xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
          xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))

          ! Store minimum/maximum coordinates
          DpointsBound(1,1,ipatch) = xmin
          DpointsBound(1,2,ipatch) = xmax

          ! Store physical coordinates of patch element corners
          DpointsReal(1,1,idxFirst:idxLast) = xmin
          DpointsReal(1,2,idxFirst:idxLast) = xmax
        end do

      case (NDIM2D)

        select case (NVE)

        case(TRIA_NVETRI2D)
          
          ! Tricky: Find minimal/maximal value and compute the lower-right
          ! and upper-left corner by hand.
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Determine minimum/maximum coordinates of the patch
            xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
            xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
            ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
            ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
            
            ! Store minimum/maximum coordinates
            DpointsBound(1,1,ipatch) = xmin
            DpointsBound(2,1,ipatch) = ymin
            DpointsBound(1,2,ipatch) = xmax
            DpointsBound(2,2,ipatch) = ymax
            
            ! Store physical coordinates of patch element corners
            DpointsReal(1,1,idxFirst:idxLast) = xmin
            DpointsReal(2,1,idxFirst:idxLast) = ymin
            DpointsReal(1,2,idxFirst:idxLast) = xmax+ymax-ymin
            DpointsReal(2,2,idxFirst:idxLast) = ymin
            DpointsReal(1,3,idxFirst:idxLast) = xmin
            DpointsReal(2,3,idxFirst:idxLast) = xmax+ymax-xmin
          end do
          
        case (TRIA_NVEQUAD2D)
          
          ! Simple: Just find minimal/maximal value
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Determine minimum/maximum coordinates of the patch
            xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
            xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
            ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
            ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
            
            ! Store minimum/maximum coordinates
            DpointsBound(1,1,ipatch) = xmin
            DpointsBound(2,1,ipatch) = ymin
            DpointsBound(1,2,ipatch) = xmax
            DpointsBound(2,2,ipatch) = ymax
            
            ! Store physical coordinates of patch element corners
            DpointsReal(1,1,idxFirst:idxLast) = xmin
            DpointsReal(2,1,idxFirst:idxLast) = ymin
            DpointsReal(1,2,idxFirst:idxLast) = xmax
            DpointsReal(2,2,idxFirst:idxLast) = ymin
            DpointsReal(1,3,idxFirst:idxLast) = xmax
            DpointsReal(2,3,idxFirst:idxLast) = ymax
            DpointsReal(1,4,idxFirst:idxLast) = xmin
            DpointsReal(2,4,idxFirst:idxLast) = ymax
          end do
          
        case DEFAULT
          call output_line ('Invalid number of vertices per elements!', &
              OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_sim')
          call sys_halt()
        end select

      case (NDIM3D)
        
        select case (NVE)
          
        case (TRIA_NVETET3D)
          ! Tricky: Find minimal/maximal value and compute corners by hand.
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Determine minimum/maximum coordinates of the patch
            xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
            xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
            ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
            ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
            zmin = minval(DpointsReal(3,:,idxFirst:idxLast))
            zmax = maxval(DpointsReal(3,:,idxFirst:idxLast))
            
            ! Store minimum/maximum coordinates
            DpointsBound(1,1,ipatch) = xmin
            DpointsBound(2,1,ipatch) = ymin
            DpointsBound(3,1,ipatch) = zmin
            DpointsBound(1,2,ipatch) = xmax
            DpointsBound(2,2,ipatch) = ymax
            DpointsBound(3,2,ipatch) = zmax

            ! Store physical coordinates of patch element corners
            DpointsReal(1,1,idxFirst:idxLast) = xmin
            DpointsReal(2,1,idxFirst:idxLast) = ymin
            DpointsReal(3,1,idxFirst:idxLast) = zmin
            DpointsReal(1,2,idxFirst:idxLast) = xmax+ymax-zmax-ymin+zmin
            DpointsReal(2,2,idxFirst:idxLast) = ymin
            DpointsReal(3,2,idxFirst:idxLast) = zmin
            DpointsReal(1,3,idxFirst:idxLast) = xmin
            DpointsReal(2,3,idxFirst:idxLast) = xmax+ymax-zmax-xmin+zmin
            DpointsReal(3,3,idxFirst:idxLast) = zmin
            DpointsReal(1,4,idxFirst:idxLast) = xmin
            DpointsReal(2,4,idxFirst:idxLast) = ymin
            DpointsReal(3,4,idxFirst:idxLast) = -xmax-ymax+zmax+xmin+ymin
          end do
          
        case (TRIA_NVEPYR3D)
          ! Tricky: Find minimal/maximal value and compute the lower-right
          ! and upper-left corner by hand. For the z-coordinate it 
          ! suffices to find the minial/maximal value
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Determine minimum/maximum coordinates of the patch
            xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
            xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
            ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
            ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
            zmin = minval(DpointsReal(3,:,idxFirst:idxLast))
            zmax = maxval(DpointsReal(3,:,idxFirst:idxLast))
            
            ! Store minimum/maximum coordinates
            DpointsBound(1,1,ipatch) = xmin
            DpointsBound(2,1,ipatch) = ymin
            DpointsBound(3,1,ipatch) = zmin
            DpointsBound(1,2,ipatch) = xmax
            DpointsBound(2,2,ipatch) = ymax
          
            ! Store physical coordinates of patch element corners
            DpointsReal(1,1,idxFirst:idxLast) = xmin
            DpointsReal(2,1,idxFirst:idxLast) = ymin
            DpointsReal(3,1,idxFirst:idxLast) = zmin
            DpointsReal(1,2,idxFirst:idxLast) = xmax+zmax-zmin
            DpointsReal(2,2,idxFirst:idxLast) = ymin
            DpointsReal(3,2,idxFirst:idxLast) = zmin
            DpointsReal(1,3,idxFirst:idxLast) = xmax+zmax-zmin
            DpointsReal(2,3,idxFirst:idxLast) = ymax+zmax-zmin
            DpointsReal(3,3,idxFirst:idxLast) = zmin
            DpointsReal(1,4,idxFirst:idxLast) = xmin
            DpointsReal(2,4,idxFirst:idxLast) = ymax+zmax-zmin
            DpointsReal(3,4,idxFirst:idxLast) = zmin
            DpointsReal(1,5,idxFirst:idxLast) = xmin
            DpointsReal(2,5,idxFirst:idxLast) = ymin
            DpointsReal(3,5,idxFirst:idxLast) = xmax+zmax-xmin
          end do
          
        case (TRIA_NVEPRIS3D)        

          ! Tricky: Find minimal/maximal value and compute the lower-right
          ! and upper-left corner by hand. For the z-coordinate it 
          ! suffices to find the minial/maximal value
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Determine minimum/maximum coordinates of the patch
            xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
            xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
            ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
            ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
            zmin = minval(DpointsReal(3,:,idxFirst:idxLast))
            zmax = maxval(DpointsReal(3,:,idxFirst:idxLast))
            
            ! Store minimum/maximum coordinates
            DpointsBound(1,1,ipatch) = xmin
            DpointsBound(2,1,ipatch) = ymin
            DpointsBound(3,1,ipatch) = zmin
            DpointsBound(1,2,ipatch) = xmax
            DpointsBound(2,2,ipatch) = ymax
            DpointsBound(3,2,ipatch) = zmax
            
            ! Store physical coordinates of patch element corners
            DpointsReal(1,1,idxFirst:idxLast) = xmin
            DpointsReal(2,1,idxFirst:idxLast) = ymin
            DpointsReal(3,1,idxFirst:idxLast) = zmin
            DpointsReal(1,2,idxFirst:idxLast) = xmax+ymax-ymin
            DpointsReal(2,2,idxFirst:idxLast) = ymin
            DpointsReal(3,2,idxFirst:idxLast) = zmin
            DpointsReal(1,3,idxFirst:idxLast) = xmin
            DpointsReal(2,3,idxFirst:idxLast) = xmax+ymax-xmin
            DpointsReal(3,3,idxFirst:idxLast) = zmin
            DpointsReal(1,4,idxFirst:idxLast) = xmin
            DpointsReal(2,4,idxFirst:idxLast) = ymin
            DpointsReal(3,4,idxFirst:idxLast) = zmax
            DpointsReal(1,5,idxFirst:idxLast) = xmax+ymax-ymin
            DpointsReal(2,5,idxFirst:idxLast) = ymin
            DpointsReal(3,5,idxFirst:idxLast) = zmax
            DpointsReal(1,6,idxFirst:idxLast) = xmin
            DpointsReal(2,6,idxFirst:idxLast) = xmax+ymax-xmin
            DpointsReal(3,6,idxFirst:idxLast) = zmax
          end do
          
        case (TRIA_NVEHEXA3D)
          
          ! Simple: Just find minimal/maximal value
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Determine minimum/maximum coordinates of the patch
            xmin = minval(DpointsReal(1,:,idxFirst:idxLast))
            xmax = maxval(DpointsReal(1,:,idxFirst:idxLast))
            ymin = minval(DpointsReal(2,:,idxFirst:idxLast))
            ymax = maxval(DpointsReal(2,:,idxFirst:idxLast))
            zmin = minval(DpointsReal(3,:,idxFirst:idxLast))
            zmax = maxval(DpointsReal(3,:,idxFirst:idxLast))
            
            ! Store minimum/maximum coordinates
            DpointsBound(1,1,ipatch) = xmin
            DpointsBound(2,1,ipatch) = ymin
            DpointsBound(3,1,ipatch) = zmin
            DpointsBound(1,2,ipatch) = xmax
            DpointsBound(2,2,ipatch) = ymax
            DpointsBound(3,2,ipatch) = zmax
            
            ! Store physical coordinates of patch element corners
            DpointsReal(1,1,idxFirst:idxLast) = xmin
            DpointsReal(2,1,idxFirst:idxLast) = ymin
            DpointsReal(3,1,idxFirst:idxLast) = zmin
            DpointsReal(1,2,idxFirst:idxLast) = xmax
            DpointsReal(2,2,idxFirst:idxLast) = ymin
            DpointsReal(3,2,idxFirst:idxLast) = zmin
            DpointsReal(1,3,idxFirst:idxLast) = xmax
            DpointsReal(2,3,idxFirst:idxLast) = ymax
            DpointsReal(3,3,idxFirst:idxLast) = zmin
            DpointsReal(1,4,idxFirst:idxLast) = xmin
            DpointsReal(2,4,idxFirst:idxLast) = ymax
            DpointsReal(3,4,idxFirst:idxLast) = zmin
            DpointsReal(1,5,idxFirst:idxLast) = xmin
            DpointsReal(2,5,idxFirst:idxLast) = ymin
            DpointsReal(3,5,idxFirst:idxLast) = zmax
            DpointsReal(1,6,idxFirst:idxLast) = xmax
            DpointsReal(2,6,idxFirst:idxLast) = ymin
            DpointsReal(3,6,idxFirst:idxLast) = zmax
            DpointsReal(1,7,idxFirst:idxLast) = xmax
            DpointsReal(2,7,idxFirst:idxLast) = ymax
            DpointsReal(3,7,idxFirst:idxLast) = zmax
            DpointsReal(1,8,idxFirst:idxLast) = xmin
            DpointsReal(2,8,idxFirst:idxLast) = ymax
            DpointsReal(3,8,idxFirst:idxLast) = zmax
          end do
          
        case DEFAULT
          call output_line ('Invalid number of vertices per elements!', &
              OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_sim')
          call sys_halt()
        end select

      case DEFAULT
        call output_line ('Invalid number of spatial dimensions!', &
            OU_CLASS_ERROR,OU_MODE_STD,'calc_patchBoundingGroup_sim')
        call sys_halt()
      end select
    end subroutine calc_patchBoundingGroup_sim

    !**************************************************************
    ! Transform physical coordinates to local coordinates of the 
    ! constant Jacobian "patch" elements which is uniquely 
    ! determined by its minimum/maximum values.
    ! 
    ! Note that for quadrilateral/hexahedral elements the x-, y- and
    ! z-coordinates are stored whereas for triangular/tetrahedral
    ! elements barycentric coordinates are adopted.
    
    subroutine calc_localTrafo_sim(IelementsInPatchIdx, icoordSystem, &
        npatches, NVE, DpointsReal, DpointsBound, DpointsRef)

      ! Index vector and list of elements in patch
      integer, dimension(:), intent(in) :: IelementsInPatchIdx

      ! Coordinate system identifier. One of the TRAFO_CS_xxxx constants. Defines
      ! the type of the coordinate system that is used for specifying the coordinates
      ! on the reference element.
      integer(I32), intent(in) :: icoordSystem
      
      ! Number of elements in patch
      integer, intent(in) :: npatches

      ! Number of vertices per element
      integer, intent(in) :: NVE
      
      ! Coordinates of the points on the physical element.
      real(DP), dimension(:,:,:), intent(in) :: DpointsReal

      ! Coordinates of the bounding group of each patch
      real(DP), dimension(:,:,:), intent(in) :: DpointsBound

      ! Coordinates of the points on the reference elements.
      real(DP), dimension(:,:,:), intent(out) :: DpointsRef
      
      ! local variables
      integer :: ipatch,idxFirst,idxLast
      real(DP) :: xmin,ymin,xmax,ymax,daux

      ! Which coordinate system should be applied
      select case (icoordSystem)

      case (TRAFO_CS_BARY2DTRI)
        ! This coordinate system can only be applied if NVE=TRIA_NVETRI2D
        if (NVE .ne. TRIA_NVETRI2D) then
          call output_line('Invalid number of corner vertices per element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'calc_localTrafo_sim')
          call sys_halt()
        end if
        
        ! Loop over all patches
        do ipatch = 1, npatches
          idxFirst = IelementsInPatchIdx(ipatch)
          idxLast  = IelementsInPatchIdx(ipatch+1)-1

          ! Get minimum/maximum coordinates of the patch
          xmin = DpointsBound(1,1,ipatch)
          ymin = DpointsBound(2,1,ipatch)
          xmax = DpointsBound(1,2,ipatch)
          ymax = DpointsBound(2,2,ipatch)

          daux = xmax-xmin+ymax-ymin

          ! Transform physical coordinates to local ones by means
          ! of a unit coordinate transformation
          DpointsRef(3,:,idxFirst:idxLast) = (DpointsReal(2,:,idxFirst:idxLast)-ymin)/daux
          DpointsRef(2,:,idxFirst:idxLast) = (DpointsReal(1,:,idxFirst:idxLast)-xmin)/daux
          DpointsRef(1,:,idxFirst:idxLast) = 1._DP-DpointsRef(2,:,idxFirst:idxLast) &
                                                  -DpointsRef(3,:,idxFirst:idxLast)
        end do
        
      case (TRAFO_CS_REF2DTRI,TRAFO_CS_REF2DQUAD,&
            TRAFO_CS_REAL2DTRI,TRAFO_CS_REAL2DQUAD,TRAFO_CS_REF1D)
        
        ! How many corner vertices do we have?
        select case (NVE)
        case (TRIA_NVELINE1D)
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Get minimum/maximum coordinates of the patch
            xmin = DpointsBound(1,1,ipatch)
            xmax = DpointsBound(1,2,ipatch)
            
            ! Transform physical coordinates to local ones by means
            ! of a natural coordinate transformation
            DpointsRef(1,:,idxFirst:idxLast) = &
                (2.0_DP*DpointsReal(1,:,idxFirst:idxLast)-(xmax+xmin))/(xmax-xmin)
          end do

        case(TRIA_NVETRI2D)
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Get minimum/maximum coordinates of the patch
            xmin = DpointsBound(1,1,ipatch)
            ymin = DpointsBound(2,1,ipatch)
            xmax = DpointsBound(1,2,ipatch)
            ymax = DpointsBound(2,2,ipatch)
            
            ! Transform physical coordinates to local ones by means
            ! of a unit coordinate transformation
            DpointsRef(1,:,idxFirst:idxLast) = &
                (DpointsReal(1,:,idxFirst:idxLast)-(xmax+ymax-ymin))/(xmax-xmin+ymax-ymin)
            DpointsRef(2,:,idxFirst:idxLast) = &
                (DpointsReal(2,:,idxFirst:idxLast)-(xmax+ymax-xmin))/(xmax-xmin+ymax-ymin)
          end do
          
        case (TRIA_NVEQUAD2D)
          
          ! Loop over all patches
          do ipatch = 1, npatches
            idxFirst = IelementsInPatchIdx(ipatch)
            idxLast  = IelementsInPatchIdx(ipatch+1)-1
            
            ! Get minimum/maximum coordinates of the patch
            xmin = DpointsBound(1,1,ipatch)
            ymin = DpointsBound(2,1,ipatch)
            xmax = DpointsBound(1,2,ipatch)
            ymax = DpointsBound(2,2,ipatch)
            
            ! Transform physical coordinates to local ones by means
            ! of a natural coordinate transformation
            DpointsRef(1,:,idxFirst:idxLast) = &
                (2.0_DP*DpointsReal(1,:,idxFirst:idxLast)-(xmax+xmin))/(xmax-xmin)
            DpointsRef(2,:,idxFirst:idxLast) = &
                (2.0_DP*DpointsReal(2,:,idxFirst:idxLast)-(ymax+ymin))/(ymax-ymin)
          end do
          
        case DEFAULT
          call output_line ('Invalid number of vertices per elements!', &
              OU_CLASS_ERROR,OU_MODE_STD,'calc_localTrafo_sim')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Unknown coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'calc_localTrafo_sim')
        call sys_halt()
      end select
    end subroutine calc_localTrafo_sim
    
    !**************************************************************    
    ! Calculate the averaged gradient values
    ! In principal, it suffices to solve the linear system
    !   <tex> $$ (P^T * P) * x = (P^T) * b $$ </tex>
    ! for the unknown $x$. However, least-squares fitting is known
    ! to be ill-conditions. An alternative is suggested in
    !
    ! J.E. Akin, Finite Element Analysis with Error Estimation
    !
    ! Akin suggests to perform a singular value decomposition (SVD)
    ! of the matrix (<tex> $ P^T * P $ </tex>) and perform back substitution.
    ! That is exactly, what is done in this subroutine.
    
    subroutine calc_patchAverages_mult(IelementsInPatchIdx, npatches, &
        Inpoints, indofTrial, Dcoefficients, Dpolynomials, Dderivatives)
      
      ! Index vector and list of elements in patch
      integer, dimension(:), intent(in) :: IelementsInPatchIdx

      ! Number of elements in patch
      integer, intent(in) :: npatches

      ! Array with numbers of sampling points per patch
      integer, dimension(:), intent(in) :: Inpoints

      ! Number of local degrees of freedom for test functions
      integer, intent(in) :: indofTrial

      ! Vector with consistent gradient values
      real(DP), dimension(:,:), intent(in) :: Dcoefficients

      ! Vector with aligned polynomial interpolants
      real(DP), dimension(:), intent(inout) :: Dpolynomials

      ! Smoothe gradient values
      real(DP), dimension(:,:,:), intent(out) :: Dderivatives

      ! local variables
      integer :: i,ipatch,icoeffFirst,icoeffLast,ipolyFirst,ipolyLast
      real(DP), dimension(indofTrial,indofTrial) :: Dv
      real(DP), dimension(indofTrial) :: Dd
      

      ! Initialise position index
      icoeffFirst = 1
      
      ! Loop over all patches
      do ipatch = 1, npatches
        
        ! Compute absolute position of coefficient and polynomials
        icoeffLast = icoeffFirst+Inpoints(ipatch)-1
        ipolyFirst = indofTrial*(icoeffFirst-1)+1
        ipolyLast  = indofTrial*icoeffLast

        ! Compute factorisation for the singular value decomposition
        call mprim_SVD_factorise(Dpolynomials(ipolyFirst:ipolyLast),&
            indofTrial, Inpoints(ipatch), Dd, Dv, .true.)
        
        ! Perform back substitution for all componenets
        do i = 1, size(Dderivatives,3)
          call mprim_SVD_backsubst1(Dpolynomials(ipolyFirst:ipolyLast), &
              indofTrial, Inpoints(ipatch), Dd, Dv, Dderivatives(:,ipatch,i), &
              Dcoefficients(icoeffFirst:icoeffLast,i), .true.)
        end do

        ! Update position index
        icoeffFirst = icoeffFirst+Inpoints(ipatch)
      end do
    end subroutine calc_patchAverages_mult

    !**************************************************************    
    ! Calculate the averaged gradient values
    ! In principal, it suffices to solve the linear system
    ! <tex> $$ (P^T * P) * x = (P^T) * b $$ </tex>
    ! for the unknown $x$. However, least-squares fitting is known
    ! to be ill-conditions. An alternative is suggested in
    !
    ! J.E. Akin, Finite Element Analysis with Error Estimation
    !
    ! Akin suggests to perform a singular value decomposition (SVD)
    ! of the matrix (<tex> $ P^T * P $ </tex>) and perform back substitution.
    ! That is exactly, what is done in this subroutine.
    
    subroutine calc_patchAverages_sim(IelementsInPatchIdx, npatches, &
        ncubp, indofTrial, Dcoefficients, Dpolynomials, Dderivatives)

      ! Index vector and list of elements in patch
      integer, dimension(:), intent(in) :: IelementsInPatchIdx

      ! Number of elements in patch
      integer, intent(in) :: npatches

      ! Number of cubature points
      integer, intent(in) :: ncubp

      ! Number of local degrees of freedom for test functions
      integer, intent(in) :: indofTrial

      ! Vector with consistent gradient values
      !
      !   Dcoefficients(ncubp, nelemPerBlock, ndim)
      !
      real(DP), dimension(:,:,:), intent(in) :: Dcoefficients

      ! Rectangular matrix with polynomial interpolants
      !
      !   Dpolynomials(indofTrial, 1, ncubp, nelemPerBlock)
      !
      real(DP), dimension(:,:,:,:), intent(inout) :: Dpolynomials

      ! Smoothed gradient values
      !
      !   Derivatives(indofTrial, npatches, ndim)
      !
      real(DP), dimension(:,:,:), intent(out) :: Dderivatives

      ! local variables
      integer :: i,ipatch,idxFirst,idxLast,npoints

      real(DP), dimension(indofTrial,indofTrial) :: Dv
      real(DP), dimension(indofTrial) :: Dd
      
      ! Loop over all patches
      do ipatch = 1, npatches
        
        ! Get first and last element index of patch
        idxFirst = IelementsInPatchIdx(ipatch)
        idxLast  = IelementsInPatchIdx(ipatch+1)-1
        npoints  = ncubp*(idxLast-idxFirst+1)

        ! Compute factorisation for the singular value decomposition
        call mprim_SVD_factorise(Dpolynomials(:,:,:,idxFirst:idxLast),&
            indofTrial, npoints, Dd, Dv, .true.)
        
        ! Perform back substitution for all componenets
        do i = 1, size(Dderivatives,3)
          call mprim_SVD_backsubst2(Dpolynomials(:,:,:,idxFirst:idxLast), &
              indofTrial, npoints, Dd, Dv, Dderivatives(:,ipatch,i), &
              Dcoefficients(:,idxFirst:idxLast,i), .true.)
        end do
      end do
    end subroutine calc_patchAverages_sim

    !**************************************************************    
    ! Initialise the cubature formula for the destination FE space

    subroutine calc_cubatureDest(celement,ncubp, Dxi, Domega)
      
      ! Element type identifier.
      integer(I32), intent(in) :: celement

      ! Number of cubature points
      integer, intent(out) :: ncubp

      ! Coordinates of the cubature points
      real(DP), dimension(:,:), intent(out) :: Dxi
      
      ! Cubature weights of the cubature points
      real(DP), dimension(:), intent(out) :: Domega
      
      select case (celement)
      case (EL_P0)
        call cub_getCubPoints(CUB_G1_T, ncubp, Dxi, Domega)
      case (EL_Q0)
        call cub_getCubPoints(CUB_G1X1, ncubp, Dxi, Domega)
      case (EL_P1)
        call cub_getCubPoints(CUB_TRZ_T, ncubp, Dxi, Domega)
      case (EL_Q1)
        call cub_getCubPoints(CUB_TRZ, ncubp, Dxi, Domega)
      case (EL_P2)
        ! Manually calculate the coordinates of the corners/midpoints on
        ! the reference element.
        Dxi(1,1)  =  1.0_DP
        Dxi(1,2)  =  0.0_DP
        Dxi(1,3)  =  0.0_DP
        
        Dxi(2,1)  =  0.0_DP
        Dxi(2,2)  =  1.0_DP
        Dxi(2,3)  =  0.0_DP
        
        Dxi(3,1)  =  0.0_DP
        Dxi(3,2)  =  0.0_DP
        Dxi(3,3)  =  1.0_DP
        
        Dxi(4,1)  =  0.5_DP
        Dxi(4,2)  =  0.5_DP
        Dxi(4,3)  =  0.0_DP
        
        Dxi(5,1)  =  0.0_DP
        Dxi(5,2)  =  0.5_DP
        Dxi(5,3)  =  0.5_DP
        
        Dxi(6,1)  =  0.5_DP
        Dxi(6,2)  =  0.0_DP
        Dxi(6,3)  =  0.5_DP
        
        ncubp = 6
        
      case (EL_Q2)
        
        ! Manually calculate the coordinates of the corners/midpoints on
        ! the reference element.
        Dxi(1,1)  =  -1.0_DP
        Dxi(1,2)  =  -1.0_DP
        
        Dxi(2,1)  =  1.0_DP
        Dxi(2,2)  =  -1.0_DP
        
        Dxi(3,1)  =  1.0_DP
        Dxi(3,2)  =  1.0_DP
        
        Dxi(4,1)  =  -1.0_DP
        Dxi(4,2)  =  1.0_DP
        
        Dxi(5,1)  =  0.0_DP
        Dxi(5,2)  =  -1.0_DP
        
        Dxi(6,1)  =  1.0_DP
        Dxi(6,2)  =  0.0_DP
        
        Dxi(7,1)  =  0.0_DP
        Dxi(7,2)  =  1.0_DP
        
        Dxi(8,1)  =  -1.0_DP
        Dxi(8,2)  =  0.0_DP
        
        Dxi(9,1)  =  0.0_DP
        Dxi(9,2)  =  0.0_DP
        
        ncubp = 9
        
      case DEFAULT 
        call output_line ('Unsupported FE space in destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'calc_cubatureDest')
        call sys_halt()
      end select
    end subroutine calc_cubatureDest
  end subroutine ppgrd_calcGradSuperPatchRecov

  !****************************************************************************

!<subroutine>

  subroutine ppgrd_calcGradLimAvgP1Q1cnf (rvectorScalar,rvectorGradient)

!<description>
    ! Calculates the recovered gradient of a scalar finite element function
    ! by means of the limited gradient averaging technique suggested by 
    ! M. Möller and D. Kuzmin. Supports conformal discretisations in arbitrary
    ! spatial dimensions with <tex>$P_1$ and $Q_1$</tex> finite elements mixed in the
    ! source and destination vectors.
!</description>

!<input>
    ! The FE solution vector. Represents a scalar FE function.
    type(t_vectorScalar), intent(in) :: rvectorScalar
!</input>

!<inputoutput>
    ! A block vector receiving the gradient.
    ! The first subvector receives the X-gradient.
    ! In 2D/3D discretisations, the 2nd subvector recevies the Y-gradient.
    ! In 3D discretisations, the 3rd subvector receives the Z-gradient.
    ! The vector must be prepared with a discretisation structure that defines
    ! the destination finite element space for the gradient field.
    type(t_vectorBlock), intent(inout) :: rvectorGradient
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j,k,icurrentElementDistr,NVE,IELmax,IELset,idof
    logical :: bnonparTrial
    integer(I32) :: icoordinatesystem

    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! Cubature formula weights. The cotent is actually not used here.
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! Number of 'cubature points'. As we acutally do not do cubature
    ! here, this coincides with the number of DOF`s on each element
    ! in the destination space.
    integer :: nlocalDOFsDest
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofDest
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsRef

    ! An array receiving the coordinates of cubature points on
    ! the real element for all elements in a set.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsReal

    ! Pointer to the point coordinates to pass to the element function.
    ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
    ! the trial element is parametric or not.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsTrial
    
    ! Array with coordinates of the corners that form the real element.
    real(DP), dimension(:,:,:), pointer :: p_Dcoords
    
    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:,:,:), pointer :: p_Djac
    
    ! Pointer to KVERT of the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
    ! Pointer to DCORVG of the triangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    ! Current element distribution in source- and destination vector
    type(t_elementDistribution), pointer :: p_relementDistribution
    type(t_elementDistribution), pointer :: p_relementDistrDest
    
    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dderivatives
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    
    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), allocatable :: IdofsTrial
    integer, dimension(:,:), allocatable :: IdofsDest
    
    ! Pointers to the X-, Y- and Z-derivative vector
    real(DP), dimension(:), pointer :: p_DxDeriv, p_DyDeriv, p_DzDeriv
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to an array that counts if an edge has been visited.
    integer :: h_IcontributionsAtDOF
    integer, dimension(:), pointer :: p_IcontributionsAtDOF
    
    ! Discretisation structures for the source- and destination vector(s)
    type(t_spatialDiscretisation), pointer :: p_rdiscrSource, p_rdiscrDest

    ! Dimension of triangulation must be less or equal than number of subvectors
    ! in the block vector.
    
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line ('No discretisation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end if
    
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) .and.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .ne. SPDISC_CONFORMAL)) then
      call output_line ('Only uniform and conformal discretisations supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end if
    
    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
      call output_line ('No triangulation attached to the source vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end if
    
    if (rvectorScalar%p_rspatialDiscr%p_rtriangulation%ndim .gt. &
        rvectorGradient%nblocks) then
      call output_line ('Dimension of destination vector not large enough!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end if
    
    ! There must be given discretisation structures in the destination vector.
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line ('No discretisation attached to the destination vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end if

    ! The source vector must be either pure Q1, P1 or mixed Q1/P1
    p_rdiscrSource => rvectorScalar%p_rspatialDiscr
    do j = 1, p_rdiscrSource%inumFESpaces
      select case (elem_getPrimaryElement (p_rdiscrSource%RelementDistr(j)%celement))
      case (EL_Q1_2D, EL_Q1_3D, EL_P1_1D, EL_P1_2D, EL_P1_3D)
        
      case DEFAULT
        call output_line ('Only Q1, and P1 supported as' // &
            ' discretisation for the source vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
        call sys_halt()
      end select
    end do
    
    ! The destination vector must be either pure Q1T or pure P1T or mixed Q1T/P1T
    do i = 1, min(rvectorGradient%nblocks,&
                  rvectorScalar%p_rspatialDiscr%ndimension)
      p_rdiscrDest => rvectorGradient%p_rblockDiscr%RspatialDiscr(i)
      do j = 1, p_rdiscrDest%inumFESpaces
        select case(elem_getPrimaryElement(p_rdiscrDest%RelementDistr(j)%celement))
        case (EL_Q1T_2D, EL_Q1T_3D, EL_P1T_2D)
          
        case DEFAULT
          call output_line ('Only Q1T, and P1T supported as' // &
              ' discretisation for the destination vector!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
          call sys_halt()
        end select
      end do
    end do
    
    ! Get the discretisation structures of the source- and destination space.
    ! Note that we assume here that the X- and Y-derivative is discretised
    ! the same way!
    p_rdiscrSource => rvectorScalar%p_rspatialDiscr
    p_rdiscrDest   => rvectorGradient%p_rblockDiscr%RspatialDiscr(1)
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscrSource%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPGRD_NELEMSIM,p_rtriangulation%NEL)
    
    ! Get a pointer to the KVERT and DCORVG array
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                               p_IverticesAtElement)
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, &
                               p_DvertexCoords)
                     
    ! Get pointers to the derivative destination vector.
    select case(p_rtriangulation%ndim)
    case (NDIM1D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lalg_clearVectorDble (p_DxDeriv)

    case (NDIM2D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      call lalg_clearVectorDble (p_DyDeriv)

    case (NDIM3D)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(1),p_DxDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(2),p_DyDeriv)
      call lsyssc_getbase_double (rvectorGradient%RvectorBlock(3),p_DzDeriv)
      call lalg_clearVectorDble (p_DxDeriv)
      call lalg_clearVectorDble (p_DyDeriv)
      call lalg_clearVectorDble (p_DzDeriv)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end select

    ! Array that allows the calculation about the number of elements
    ! meeting in a vertex, based onm DOF`s.
    call storage_new ('ppgrd_calcGradLimAvgP1Q1cnf','p_IcontributionsAtDOF',&
                      dof_igetNDofGlob(p_rdiscrDest),&
                      ST_INT, h_IcontributionsAtDOF, ST_NEWBLOCK_ZERO)
    call storage_getbase_int (h_IcontributionsAtDOF,p_IcontributionsAtDOF)
    
    ! Evaluate the first derivative of the FE functions.

    Bder = .false.
    select case (p_rtriangulation%ndim)
    case (NDIM1D)
      Bder(DER_DERIV1D_X) = .true.

    case (NDIM2D)
      Bder(DER_DERIV2D_X) = .true.
      Bder(DER_DERIV2D_Y) = .true.

    case (NDIM3D)
      Bder(DER_DERIV3D_X) = .true.
      Bder(DER_DERIV3D_Y) = .true.
      Bder(DER_DERIV3D_Z) = .true.

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
      call sys_halt()
    end select

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,p_rdiscrSource%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => p_rdiscrSource%RelementDistr(icurrentElementDistr)
      p_relementDistrDest => p_rdiscrDest%RelementDistr(icurrentElementDistr)
    
      ! If the element distribution is empty, skip it
      if (p_relementDistribution%NEL .eq. 0) cycle
    
      ! Get the number of local DOF`s for trial functions
      ! in the source and destination vector.
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      indofDest = elem_igetNDofLoc(p_relementDistrDest%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)

      ! Initialise the cubature formula.
      ! That is a special trick here! The FE space of the destination vector
      ! is either P1 or Q1. We create the gradients in the midpoints of the edges
      ! by taking the limited average of the gradients of the source vector!

      select case (elem_getPrimaryElement(p_relementDistrDest%celement))
      case (EL_P1T)
        call cub_getCubPoints(CUB_G3_T, nlocalDOFsDest, Dxi, Domega)

      case (EL_Q1T)
        Dxi(1,1) =  0.0_DP
        Dxi(1,2) = -1.0_DP
        Dxi(2,1) =  1.0_DP
        Dxi(2,2) =  0.0_DP
        Dxi(3,1) =  0.0_DP
        Dxi(3,2) =  1.0_DP
        Dxi(4,1) = -1.0_DP
        Dxi(4,2) =  0.0_DP
        
        nlocalDOFsDest = 4

      case DEFAULT
        call output_line ('Unsupported FE space in destination vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
        call sys_halt()
      end select

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      icoordinatesystem = elem_igetCoordSystem(p_relementDistribution%celement)
      
      ! Allocate memory and get local references to it.
      ! We abuse the system of cubature points here for the evaluation.
      call domint_initIntegration (rintSubset,nelementsPerBlock,nlocalDOFsDest,&
        icoordinatesystem,p_rtriangulation%ndim,NVE)
      p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
      p_DcubPtsReal => rintSubset%p_DcubPtsReal
      p_Djac =>        rintSubset%p_Djac
      p_Ddetj =>       rintSubset%p_Ddetj
      p_Dcoords =>     rintSubset%p_DCoords

      ! Destination space is either P1 or Q1.
      ! Put the cubature point coordinates in the right format to the
      ! cubature-point array.
      ! Initialise all entries in p_DcubPtsRef with the same coordinates -
      ! as the cubature point coordinates are identical on all elements
      do j=1,size(p_DcubPtsRef,3)
        do i=1,nlocalDOFsDest
          do k=1,size(p_DcubPtsRef,1)
            ! Could be solved using the TRANSPOSE operator - but often it is
            ! faster this way...
            p_DcubPtsRef(k,i,j) = Dxi(i,k)
          end do
        end do
      end do

      ! Allocate memory for the DOF`s of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))
      allocate(IdofsDest(indofDest,nelementsPerBlock))

      ! Allocate memory for the values of the derivatives in the corners
      allocate(Dderivatives(nlocalDOFsDest,nelementsPerBlock,2))

      ! Check if one of the trial/test elements is nonparametric
      bnonparTrial  = elem_isNonparametric(p_relementDistribution%celement)
                      
      ! Let p_DcubPtsTest point either to p_DcubPtsReal or
      ! p_DcubPtsRef - depending on whether the space is parametric or not.
      if (bnonparTrial) then
        p_DcubPtsTrial => p_DcubPtsReal
      else
        p_DcubPtsTrial => p_DcubPtsRef
      end if

      ! Get the number of elements in the element distribution.
      NEL = p_relementDistribution%NEL
      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPGRD_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPGRD_NELEMSIM)
      
        ! Calculate the global DOF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscrSource, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Also calculate the global DOF`s in our destination vector(s)
        call dof_locGlobMapping_mult(p_rdiscrDest, p_IelementList(IELset:IELmax), &
                                     IdofsDest)

        ! We have the coordinates of the cubature points saved in the
        ! coordinate array from above. Unfortunately for nonparametric
        ! elements, we need the real coordinate.
        ! Furthermore, we anyway need the coordinates of the element
        ! corners and the Jacobian determinants corresponding to
        ! all the points.
        !
        ! At first, get the coordinates of the corners of all the
        ! elements in the current set. 
        
        call trafo_getCoords_sim (elem_igetTrafoType(&
            p_relementDistribution%celement),&
            p_rtriangulation,p_IelementList(IELset:IELmax),p_Dcoords)
        
        ! Depending on the type of transformation, we must now choose
        ! the mapping between the reference and the real element.
        ! In case we use a nonparametric element, we need the 
        ! coordinates of the points on the real element, too.
        call trafo_calctrafo_sim (&
              p_rdiscrSource%RelementDistr(icurrentElementDistr)%ctrafoType,&
              IELmax-IELset+1,nlocalDOFsDest,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1),&
              p_DcubPtsReal)
      
        ! Prepare the call to the evaluation routine of the analytic function.    
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        
        ! How many spatial dimensions are we?
        select case (p_rtriangulation%ndim)
        case (NDIM1D)
          ! At this point, we calculate the gradient information and
          ! Calculate the X-derivative in the corners of the elements
          ! into Dderivatives(:,:,1).
          call fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_relementDistribution%celement, IdofsTrial, &
              nlocalDOFsDest, int(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV1D_X,&
              Dderivatives(:,1:IELmax-IELset+1,1))

          ! Sum up the derivative values in the destination vector.
          ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
          ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
          ! local DOF`s in the destination space -- in that order!
          do i=1,IELmax-IELset+1
            do j=1,nlocalDOFsDest

              idof = IdofsDest(j,i)

              if (p_IcontributionsAtDOF(idof) .eq. 0) then 
                p_DxDeriv(idof) = Dderivatives(j,i,1)
                p_IcontributionsAtDOF(idof) = 1
              else
                p_DxDeriv(idof) = (sign(0.5_DP,p_DxDeriv(idof)) + &
                                   sign(0.5_DP,Dderivatives(j,i,1))) * &
                                   min(0.5_DP*abs(p_DxDeriv(idof)+Dderivatives(j,i,1)), &
                                       2._DP*abs(p_DxDeriv(idof)), &
                                       2._DP*abs(Dderivatives(j,i,1)))
              end if
            end do
          end do

        case (NDIM2D)
          ! At this point, we calculate the gradient information and
          ! Calculate the X-derivative in the corners of the elements
          ! into Dderivatives(:,:,1) and the Y-derivatives into
          ! Dderivatives(:,:,2)
          call fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_relementDistribution%celement, IdofsTrial, &
              nlocalDOFsDest, int(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV2D_X,&
              Dderivatives(:,1:IELmax-IELset+1,1))
          
          call fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_relementDistribution%celement, IdofsTrial, &
              nlocalDOFsDest, int(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV2D_Y,&
              Dderivatives(:,1:IELmax-IELset+1,2))
          
          ! Sum up the derivative values in the destination vector.
          ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
          ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
          ! local DOF`s in the destination space -- in that order!
          do i=1,IELmax-IELset+1
            do j=1,nlocalDOFsDest

              idof = IdofsDest(j,i)

              if (p_IcontributionsAtDOF(idof) .eq. 0) then 
                p_DxDeriv(idof) = Dderivatives(j,i,1)
                p_DyDeriv(idof) = Dderivatives(j,i,2)
                p_IcontributionsAtDOF(idof) = 1
              else
                p_DxDeriv(idof) = (sign(0.5_DP,p_DxDeriv(idof)) + &
                                   sign(0.5_DP,Dderivatives(j,i,1))) * &
                                   min(0.5_DP*abs(p_DxDeriv(idof)+Dderivatives(j,i,1)), &
                                       2._DP*abs(p_DxDeriv(idof)), &
                                       2._DP*abs(Dderivatives(j,i,1)))
                p_DyDeriv(idof) = (sign(0.5_DP,p_DyDeriv(idof)) + &
                                   sign(0.5_DP,Dderivatives(j,i,2))) * &
                                   min(0.5_DP*abs(p_DyDeriv(idof)+Dderivatives(j,i,2)), &
                                       2._DP*abs(p_DyDeriv(idof)), &
                                       2._DP*abs(Dderivatives(j,i,2)))
              end if
            end do
          end do

        case (NDIM3D)
          ! At this point, we calculate the gradient information and
          ! Calculate the X-derivative in the corners of the elements
          ! into Dderivatives(:,:,1), the Y-derivatives into Dderivatives(:,:,2)
          ! and the Z-derivatives into Dderivatives(:,:,3).
          call fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_relementDistribution%celement, IdofsTrial, &
              nlocalDOFsDest, int(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV3D_X,&
              Dderivatives(:,1:IELmax-IELset+1,1))
          
          call fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_relementDistribution%celement, IdofsTrial, &
              nlocalDOFsDest, int(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV3D_Y,&
              Dderivatives(:,1:IELmax-IELset+1,2))

          call fevl_evaluate_sim (rvectorScalar, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              p_relementDistribution%celement, IdofsTrial, &
              nlocalDOFsDest, int(IELmax-IELset+1), p_DcubPtsTrial, DER_DERIV3D_Z,&
              Dderivatives(:,1:IELmax-IELset+1,3))

          ! Sum up the derivative values in the destination vector.
          ! Note that we explicitly use the fact, that the each pair of nlocalDOFsDest 
          ! 'cubature points', or better to say 'corners'/'midpoints', coincides with the 
          ! local DOF`s in the destination space -- in that order!
          do i=1,IELmax-IELset+1
            do j=1,nlocalDOFsDest

              idof = IdofsDest(j,i)

              if (p_IcontributionsAtDOF(idof) .eq. 0) then 
                p_DxDeriv(idof) = Dderivatives(j,i,1)
                p_DyDeriv(idof) = Dderivatives(j,i,2)
                p_DzDeriv(idof) = Dderivatives(j,i,3)
                p_IcontributionsAtDOF(idof) = 1
              else
                p_DxDeriv(idof) = (sign(0.5_DP,p_DxDeriv(idof)) + &
                                   sign(0.5_DP,Dderivatives(j,i,1))) * &
                                   min(0.5_DP*abs(p_DxDeriv(idof)+Dderivatives(j,i,1)), &
                                       2._DP*abs(p_DxDeriv(idof)), &
                                       2._DP*abs(Dderivatives(j,i,1)))
                p_DyDeriv(idof) = (sign(0.5_DP,p_DyDeriv(idof)) + &
                                   sign(0.5_DP,Dderivatives(j,i,2))) * &
                                   min(0.5_DP*abs(p_DyDeriv(idof)+Dderivatives(j,i,2)), &
                                       2._DP*abs(p_DyDeriv(idof)), &
                                       2._DP*abs(Dderivatives(j,i,2)))
                p_DzDeriv(idof) = (sign(0.5_DP,p_DzDeriv(idof)) + &
                                   sign(0.5_DP,Dderivatives(j,i,3))) * &
                                   min(0.5_DP*abs(p_DzDeriv(idof)+Dderivatives(j,i,3)), &
                                       2._DP*abs(p_DzDeriv(idof)), &
                                       2._DP*abs(Dderivatives(j,i,3)))
              end if
            end do
          end do

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ppgrd_calcGradLimAvgP1Q1cnf')
          call sys_halt()
        end select
        
      end do ! IELset
      
      ! Release memory
      call domint_doneIntegration(rintSubset)
      
      deallocate(Dderivatives)
      deallocate(IdofsDest)
      deallocate(IdofsTrial)
      
    end do ! icurrentElementDistr

    ! Release temp data
    call storage_free (h_IcontributionsAtDOF)

  end subroutine ppgrd_calcGradLimAvgP1Q1cnf

end module

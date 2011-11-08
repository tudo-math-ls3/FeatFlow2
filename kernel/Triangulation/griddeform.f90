!##############################################################################
!# ****************************************************************************
!# <name> griddeform </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# <todo> Mal andere Elementtypen ausprobieren, eine Callbackfunction für
!#        die Monitorfunktion
!#
!# 1.) griddef_deformationInit
!#     -> Initialise a griddeform structure
!#
!# 2.) griddef_DeformationDone
!#     -> release a griddeform structure
!#
!# 3.) griddef_freeLevels
!#     -> release the resources we needed for
!#        the grid deformation levels
!#
!# 4.) griddef_calcMonitorFunction
!#     -> creates a monitor function from a given error distribution
!#
!# 5.) griddef_prepareDeformation
!#     -> allocate needed resources
!#
!# 6.) griddef_performDeformation
!#     -> main routine for the grid deformation process
!#
!# 6a.)griddef_computeAdapSteps
!#     -> computes the number of adaption steps
!#
!# 7.) griddef_performOneDefStep
!#     -> performs one deformation step of the (basic,enhanced) deformation method
!#
!# 8.) griddef_getArea
!#     -> auxilliary function to get the area in all the elements,
!#     -> put it in a Q0 vector and then project the Q0 Vector to a Q1 vector
!#
!# 9.) griddef_buildMonFuncTest
!#    -> build a test monitor function
!#
!# 10.)  griddef_normaliseFctsNum
!#     -> normalise our functions
!#
!# 11.) griddef_blendmonitor
!#     -> blend the monitor function
!#
!# 12.) griddef_normaliseFctsInv
!#     -> normalise the function inverse
!#
!# 13.) griddef_normaliseFctsInvAux
!#     -> auxilliary routine for the function inverse
!#
!# 14.) griddef_createMatrixDef
!#     -> create the matrix of the deformation poisson problem
!#
!# 15.) griddef_createRHS
!#     -> create the rhs of the deformation poisson problem
!#
!# 16.) griddef_moveMesh
!#     -> solve the ode and move the mesh
!#
!# 17.) griddef_performEE
!#     -> performs Explicit euler to solve the ODE for all
!#        innter vertices
!#
!# 18.) griddef_perform_boundary2
!#     -> performs the explicit euler on the boundary
!#        vertices.
!#
!# 19.) griddef_evalPhi_Known
!#    -> auxilliary routine for the ODE solver
!#
!# 20.) griddef_evalphi_ray
!#    -> auxilliary routine for the ODE solver
!#
!# 21.) griddef_evalphi_ray_bound
!#    -> auxilliary routine for the ODE solver
!#
!# 22.) griddef_getAreaDeformed
!#    -> compute the area distribution in the deformed, mainly
!#       for debugging purposes
!#
!# 23.) griddef_qMeasureM1
!#    -> compute the quality measure m1
!#
!# 24.) griddef_buildHGrid
!#    -> ?
!#
!# </purpose>
!##############################################################################


module griddeform

  use basicgeometry
  use bcassembly
  use bilinearformevaluation
  use boundary
  use cubature
  use derivatives
  use dofmapping
  use element
  use elementpreprocessing
  use feevaluation
  use filtersupport
  use fsystem
  use genoutput
  use geometryaux
  use linearformevaluation
  use linearsolver
  use linearsystemblock
  use linearsystemscalar
  use matrixfilters
  use perfconfig
  use pprocerror
  use pprocgradients
  use spatialdiscretisation
  use spdiscprojection
  use stdoperators
  use storage
  use transformation
  use triangulation
  use triasearch
  use ucd
  use vectorfilters
  
  implicit none
  
  private

!<constants>

  integer, parameter, public :: GRIDDEF_CLASSICAL  = 0

  integer, parameter, public :: GRIDDEF_MULTILEVEL = -1
  
  integer, parameter, public :: GRIDDEF_FIXED      = 2
  
  integer, parameter, public :: GRIDDEF_USER       = 3

!</constants>

!<types>

  !<typeblock>
  ! This type contains the different levels of a multigrid
  ! so that they can be efficiently used in the deformation routine
  type t_hgridLevels
  
    ! An pointer to the original triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! Triangulation used during the mesh defornation
    type(t_triangulation) :: rtriangulation
    
    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixScalar) :: rmatrix
    
    type(t_matrixBlock) :: rmatDeform
    
    ! A block vector to store the are distribution
    type(t_vectorBlock)  :: rvectorAreaBlockQ1
        ! A block vector where we store the monitor function
    type(t_vectorBlock)  :: rvectorMonFuncQ1
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.

    ! handle for solution of Poisson problem in deformation method
    type(t_vectorBlock) :: rSolBlock
    ! handle for rhs of Poisson problem in deformation method
    type(t_vectorBlock) :: rrhsBlock
    ! handle for rhs of Poisson problem in deformation method
    type(t_vectorBlock) :: rtempBlock
    ! block vector for the recovered gradient
    type(t_vectorBlock) :: rvecGradBlock
    ! tells whether the structure needs to be reinitialised
    logical             :: breinit = .false.

  end type
  
  public :: t_hgridLevels
  
  !<\typeblock>
  
  !<typeblock>
  ! This type contains everything which is necessary to define
  ! a certain deformation algorithm.
  type t_griddefInfo

    !
    type(t_hgridLevels), dimension(:), pointer :: p_rhLevels => null()

    ! this is a Pointer to the original triangulation structure,
    ! when our deformation was a success, we will overwrite
    ! this structure with the deformed vertex coordinates
    ! type(t_triangulation), pointer :: p_rtriangulation => NULL()
    
    ! Pointer to a triangulation structure, we use this one
    ! as our local copy, in case of a successful deformation
    ! it will be become the new triangulation
    ! type(t_triangulation) :: rDeftriangulation
    
    ! the boundary information of the undeformed grid
    type(t_boundary), pointer :: p_rboundary => null()

    ! description missing
    integer:: cdefStyle

    ! number of time steps in ODE solver
    integer:: ntimeSteps

    ! type of ODE solver
    integer:: codeMethod

    ! regularisation parameter for monitor function
    ! if dregpar .lt. 0: dregpar is weighted by <tex>$0.5^{maxMGLevel}$</tex>
    real(dp) :: dregpar

    ! form parameter for monitor function
    real(dp) ::dexponent
    
    ! Blending parameter for basic deformation
    real(dp) :: dblendPar

    ! number of adaptation steps
    integer :: nadaptionSteps

    ! admissible quality number (stopping criterion for correction)
    real(dp) :: dtolQ = ST_NOHANDLE

    ! description missing
    integer :: h_Dql2 = ST_NOHANDLE
    
    ! description missing
    integer :: h_Dqlinfty = ST_NOHANDLE

    ! number of adaption steps during deformation
    integer :: h_nadapSteps = ST_NOHANDLE

    ! description missing
    integer :: h_calcAdapSteps = ST_NOHANDLE

    ! number of adaption steps during deformation really performed
    ! (may differ from nadapSteps: if nadapSteps .lt.0, then compute
    ! nadapSteps according to the deformation problem)
    integer :: h_nadapStepsReally = ST_NOHANDLE

    ! maximal number of correction steps
    integer :: h_nmaxCorrSteps = ST_NOHANDLE

    ! number of correction steps really performed
    integer :: h_ncorrSteps = ST_NOHANDLE

    ! ODE level for deformation
    integer :: h_ilevelODE = ST_NOHANDLE

    ! ODE level for correction
    integer :: h_ilevelODECorr = ST_NOHANDLE

    ! number of smoothing steps for creating the monitor function
    ! if .lt. 0: values are multiplied by current level
    integer :: h_npreSmoothSteps = ST_NOHANDLE
    ! number of smoothing steps for creating the monitor function
    ! if .lt. 0: values are multiplied by current level
    integer :: h_npostSmoothSteps = ST_NOHANDLE
    ! number of smoothing steps for creating the monitor function
    ! if .lt. 0: values are multiplied by current level
    integer :: h_nintSmoothSteps = ST_NOHANDLE

    ! number of smoothing steps for monitor function
    integer :: h_nmonSmoothsteps = ST_NOHANDLE

    ! smoothing method used:
    integer :: h_cpreSmoothMthd = ST_NOHANDLE
    ! smoothing method used:
    integer :: h_cpostSmoothMthd = ST_NOHANDLE
    ! smoothing method used:
    integer :: h_cintSmoothMthd = ST_NOHANDLE

    ! handle for the vector containing the blending parameters
    integer:: h_DblendPar = ST_NOHANDLE

    ! method for gradient recovery (1: SPR, 2 : PPR, 3: interpolation)
    integer:: crecMethod

    ! if true, use the FEM-interpolant of the monitor function instead of the
    ! monitor function directly
    integer:: cMonFctEvalMethod

    ! if true, use the error distribution to obtain the monitor function
    logical ::buseErrorDistr

    ! description missing
    real(DP) :: damplification

    ! if true, perform relative deformation
    logical :: brelDef

    ! description missing
    integer:: iminDefLevel

    ! description missing
    integer:: idefLevIncr

    ! scaling factor for the monitor function itself such that
    ! \int dscalefmon f = |\Omega|
    real(DP) :: dscalefmon

    ! same for the reciprocal of the monitor function
    real(DP) :: dscalefmonInv

    ! Feat-vector with boundary description
    integer:: h_IboundCompDescr_handles = ST_NOHANDLE
    
    ! save the minimum level
    integer :: NLMIN
    
    ! save the maximum level
    integer :: NLMAX

  end type t_griddefInfo
  
  public :: t_griddefInfo
  !</typeblock>
  
!</types>

  !************************************************************************
  
  ! global performance configuration
  type(t_perfconfig), target, save :: griddef_perfconfig

  !************************************************************************

  public :: griddef_deformationInit
  public :: griddef_DeformationDone
  public :: griddef_calcMonitorFunction
  public :: griddef_prepareDeformation
  public :: griddef_performDeformation
  public :: griddef_computeAdapSteps
  public :: griddef_performOneDefStep
  public :: griddef_getArea
  public :: griddef_buildMonFuncTest
  public :: griddef_normaliseFctsNum
  public :: griddef_blendmonitor
  public :: griddef_normaliseFctsInv
  public :: griddef_normaliseFctsInvAux
  public :: griddef_createMatrixDef
  public :: griddef_createRHS
  public :: griddef_moveMesh
  public :: griddef_performEE
  public :: griddef_perform_boundary2
  public :: griddef_evalPhi_Known
  public :: griddef_evalphi_ray
  public :: griddef_evalphi_ray_bound
  public :: griddef_getAreaDeformed
  public :: griddef_qMeasureM1
  public :: griddef_buildHGrid
  public :: griddef_deformationInit3D
  public :: griddef_performDeformation3D

contains

!<subroutine>

  subroutine griddef_deformationInit(rgriddefInfo,NLMIN,NLMAX,rboundary,iStyle,iadaptSteps,iodesteps)
  
!<description>
  ! This subroutine initialises the rgriddefinfo structure which describes the
  ! grid deformation algorithm. The values for the certain parameters defining
  ! the deformation process are read in from the master.dat file. This routine
  ! has to be called before starting grid deformation.
!</description>

!<inputoutput>
  ! We will fill this structure with useful values
  type(t_griddefInfo), intent(inout) :: rgriddefInfo
  
  ! the boundary information
  type(t_boundary), intent(inout), target :: rboundary
  
  ! user defined number of adaptation steps
  integer,intent(in), optional :: iadaptSteps
  
  ! user defined number of adaptation steps
  integer,intent(in), optional :: iodesteps
    
!</inputoutput>

!<input>
  integer, intent(in) :: NLMAX
  
  integer, intent(in) :: NLMIN
  
  ! the deformation method to use multilevel or classical
  ! GRIDDEF_CLASSICAL or GRIDDEF_MULTILEVEL
  integer, intent(in) :: iStyle
!</input>
!</subroutine>
  
  ! Local variables
  integer :: iaux
  ! pointers to the arrays, we need to fill
  ! with useful default parameters
  integer, dimension(:), pointer :: p_Dql2
  integer, dimension(:), pointer :: p_Dqlinfty
  integer, dimension(:), pointer :: p_nadapSteps
  integer, dimension(:), pointer :: p_calcAdapSteps
  integer, dimension(:), pointer :: p_nadapStepsReally
  integer, dimension(:), pointer :: p_nmaxCorrSteps
  integer, dimension(:), pointer :: p_ncorrSteps
  integer, dimension(:), pointer :: p_ilevelODE
  integer, dimension(:), pointer :: p_ilevelODECorr
  integer, dimension(:), pointer :: p_npreSmoothSteps
  integer, dimension(:), pointer :: p_npostSmoothSteps
  integer, dimension(:), pointer :: p_nintSmoothSteps
  integer, dimension(:), pointer :: p_nmonSmoothsteps
  integer, dimension(:), pointer :: p_cpreSmoothMthd
  integer, dimension(:), pointer :: p_cpostSmoothMthd
  integer, dimension(:), pointer :: p_cintSmoothMthd
  
  rgriddefInfo%p_rboundary => rboundary
  
  rgriddefInfo%dblendPar = 1.0_dp
  
  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_Dql2',&
                   NLMAX, ST_INT, rgriddefInfo%h_Dql2,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_Dqlinfty',&
                   NLMAX, ST_INT, rgriddefInfo%h_Dqlinfty,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_nadapSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nadapSteps,&
                   ST_NEWBLOCK_ZERO)
                   
  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_calcAdapSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_calcAdapSteps,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_nadapStepsReally',&
                   NLMAX, ST_INT, rgriddefInfo%h_nadapStepsReally,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_nmaxCorrSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nmaxCorrSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_ncorrSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_ncorrSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_ilevelODE',&
                   NLMAX, ST_INT, rgriddefInfo%h_ilevelODE,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_ilevelODECorr',&
                   NLMAX, ST_INT, rgriddefInfo%h_ilevelODECorr,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_npreSmoothSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_npreSmoothSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_npostSmoothSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_npostSmoothSteps,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_nintSmoothSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nintSmoothSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_nmonSmoothsteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nmonSmoothsteps,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_cpreSmoothMthd',&
                   NLMAX, ST_INT, rgriddefInfo%h_cpreSmoothMthd,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_cpostSmoothMthd',&
                   NLMAX, ST_INT, rgriddefInfo%h_cpostSmoothMthd,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit', 'rgriddefInfo%h_cintSmoothMthd',&
                   NLMAX, ST_INT, rgriddefInfo%h_cintSmoothMthd,&
                   ST_NEWBLOCK_ZERO)

    
  
  ! Now get all the just allocated pointers
  call storage_getbase_int(rgriddefInfo%h_Dql2,p_Dql2)
  
  call storage_getbase_int(rgriddefInfo%h_Dqlinfty,p_Dqlinfty)

  call storage_getbase_int(rgriddefInfo%h_nadapSteps,p_nadapSteps)
  
  call storage_getbase_int(rgriddefInfo%h_calcAdapSteps,p_calcAdapSteps)
  
  call storage_getbase_int(rgriddefInfo%h_nadapStepsReally,p_nadapStepsReally)
  
  call storage_getbase_int(rgriddefInfo%h_nmaxCorrSteps,p_nmaxCorrSteps)
  
  call storage_getbase_int(rgriddefInfo%h_ncorrSteps,p_ncorrSteps)

  call storage_getbase_int(rgriddefInfo%h_ilevelODE,p_ilevelODE)

  call storage_getbase_int(rgriddefInfo%h_ilevelODECorr,p_ilevelODECorr)
  
  call storage_getbase_int(rgriddefInfo%h_npreSmoothSteps,p_npreSmoothSteps)

  call storage_getbase_int(rgriddefInfo%h_npostSmoothSteps,p_npostSmoothSteps)

  call storage_getbase_int(rgriddefInfo%h_nintSmoothSteps,p_nintSmoothSteps)
    
  call storage_getbase_int(rgriddefInfo%h_nmonSmoothsteps,p_nmonSmoothsteps)
      
  call storage_getbase_int(rgriddefInfo%h_cpreSmoothMthd,p_cpreSmoothMthd)

  call storage_getbase_int(rgriddefInfo%h_cpostSmoothMthd,p_cpostSmoothMthd)
  
  call storage_getbase_int(rgriddefInfo%h_cintSmoothMthd,p_cintSmoothMthd)


  ! Set the deformation style to classical
  rgriddefInfo%cdefStyle = iStyle

  ! temp variable
  iaux = 0
  
  ! set the multigrid level for the deformation PDE
  rgriddefInfo%iminDefLevel = NLMAX
  
  ! Here we initialize the structure with the standard values according to
  ! the desired grid deformation method.
  select case(rgriddefInfo%cdefStyle)
  
      case(GRIDDEF_CLASSICAL)

        ! initialise... here should be the absolute PDE Level
        ! where the deformation takes place.
        iaux = 1
      
        ! initialise the tolerance
        rgriddefInfo%dtolQ = 1.0E10_DP
        
        if(present(iodesteps))then
        ! number of ode steps
          rgriddefInfo%ntimeSteps = iodesteps
        else
          rgriddefInfo%ntimeSteps = 20
        end if
        
        ! set the multigrid level for the deformation PDE(really ?)
        p_ilevelODE(rgriddefInfo%iminDefLevel) = rgriddefInfo%iminDefLevel

        ! number of smoothing steps for creating the monitor function
        p_nintSmoothSteps(rgriddefInfo%iminDefLevel) = 0
        
        ! smoothing method used
        p_CintSmoothMthd(rgriddefInfo%iminDefLevel) = 0

        ! number of smoothing steps for monitor function
        p_NMonSmoothSteps(rgriddefInfo%iminDefLevel) = 0
        
        ! set standard value for number of correction steps
        p_nmaxCorrSteps(iaux) = 0
        
        ! set standard value for number of ODE correction steps
        p_ilevelODECorr(iaux) = 0
        
        ! initialize the blending parameter
        rgriddefInfo%dblendPar = 1.0_dp
        
        if(present(iadaptSteps))then
          ! initialize the number of adaptation steps
          rgriddefInfo%nadaptionSteps = iadaptSteps
          
          ! set Adaptive control
          ! the number of steps on all levels is user defined
          p_calcAdapSteps(:) = GRIDDEF_USER
        else
          ! initialize the number of adaptation steps
          rgriddefInfo%nadaptionSteps = 1
          
          ! the number of steps on all levels is
          ! fixed
          p_calcAdapSteps(:) = GRIDDEF_FIXED
        end if
        
     case(GRIDDEF_MULTILEVEL)
     
     ! do something here in the future
      call output_line ('MULTILEVEL method not yet implemented', &
          OU_CLASS_ERROR,OU_MODE_STD,'griddef_deformationInit')
      call sys_halt()
        
  end select
  
  ! get number of correction steps on level iaux
  p_nmaxCorrSteps(rgriddefInfo%iminDefLevel) = 0

  ! get the minimum deformation level
  p_ilevelODECorr(rgriddefInfo%iminDefLevel) = rgriddefInfo%iminDefLevel
  
  ! store the value of NLMIN
  rgriddefInfo%NLMIN = NLMIN

  ! store the value of NLMAX
  rgriddefInfo%NLMAX = NLMAX
  
  ! Allocate memory for all the levels.
  allocate(rgriddefInfo%p_rhLevels(NLMIN:NLMAX))

  end subroutine

!****************************************************************************************
  
!<subroutine>
  subroutine griddef_buildHGrid(rgriddefInfo,rtriangulation,iLevel)
!<description>
  ! In this routine we build the HGridstructure that
  ! represent the levels of the grid, this structure has an own working copy
  ! of all the vertices and boundary parameter values of the original grid.
  ! The other grid information is merely shared and not copied.
!</description>

!<inputoutput>
  ! We will fill this structure with useful values
  type(t_griddefInfo), intent(inout) :: rgriddefInfo
  
  ! The underlying triangulation
  type(t_triangulation), intent(inout), target :: rtriangulation
  
!</inputoutput>

!<input>
  integer, intent(in) :: iLevel
!</input>

!</subroutine>
  
  ! Local variables
  integer :: idupFlag
  
  
  ! assign the pointer on this level
  rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation => rtriangulation
  
  ! we dublicate all levels of the grid hierachy and
  ! share the vertices, only on the level NLMAX we
  ! actually copy the vertex coordinates
  ! Set it up to share all
  idupFlag = TR_SHARE_ALL
  if(iLevel .lt. rgriddefInfo%NLMIN)then
    !nothing
  else if((iLevel .ge. rgriddefInfo%NLMIN).and.(iLevel .lt. rgriddefInfo%NLMAX))then
    ! share all
    call tria_duplicate(rtriangulation,&
         rgriddefInfo%p_rhLevels(iLevel)%rtriangulation,idupFlag)
  else if(iLevel .eq. rgriddefInfo%NLMAX)then
    ! we only want to dublicate the vertexcoords
    ! and boundary parameters
    idupFlag = idupFlag - TR_SHARE_DVERTEXCOORDS
    idupFlag = idupFlag - TR_SHARE_DVERTEXPARAMETERVALUE

    ! We copy only the vertex coordinates and the boundary parameter values
    call tria_duplicate(rtriangulation,&
         rgriddefInfo%p_rhLevels(iLevel)%rtriangulation,idupFlag)
  else
    call output_line ('Level input not valid', &
        OU_CLASS_ERROR,OU_MODE_STD,'griddef_buildHGrid')
    call sys_halt()
  end if
 
  end subroutine ! griddef_buildHGrid

!****************************************************************************************

!<subroutine>
  subroutine griddef_DeformationDone(rgriddefInfo)
!<description>
  !
  ! Here we release the rgriddefInfo structure
  ! and all the elements that are contained
  ! within it.
  !
!</description>
  
!<inputoutput>
  ! We will fill this structure with useful values
  type(t_griddefInfo), intent(inout) :: rgriddefInfo
!</inputoutput>
!</subroutine>

  ! deallocate memory
  if (rgriddefInfo%h_Dql2 .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_Dql2)
  end if
  
  ! deallocate memory
  if (rgriddefInfo%h_Dqlinfty .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_Dqlinfty)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_calcAdapSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_calcAdapSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nadapSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nadapSteps)
  end if


  ! deallocate memory
  if (rgriddefInfo%h_nadapStepsReally .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nadapStepsReally)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nmaxCorrSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nmaxCorrSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_ncorrSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_ncorrSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_ilevelODE .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_ilevelODE)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_ilevelODECorr .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_ilevelODECorr)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_npreSmoothSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_npreSmoothSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_npostSmoothSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_npostSmoothSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nintSmoothSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nintSmoothSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nmonSmoothsteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nmonSmoothsteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_cpreSmoothMthd .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_cpreSmoothMthd)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_cpostSmoothMthd .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_cpostSmoothMthd)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_cintSmoothMthd .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_cintSmoothMthd)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_DblendPar .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_DblendPar)
  end if
  
  ! release the memory for the levels
  call griddef_freeLevels(rgriddefInfo%p_rhLevels,rgriddefInfo%NLMIN,rgriddefInfo%NLMAX)
  
  ! deallocate the levels structure
  deallocate(rgriddefInfo%p_rhLevels)
  
  end subroutine
  
!****************************************************************************************
!<subroutine>
  subroutine griddef_freeLevels(rgridLevels,inlmin,inlmax)
!<description>
  !
  ! Here we clean the memory used by the t_hgridLevels structure
  !
!</description>
  
!<inputoutput>
  ! the grid levels
  type(t_hgridLevels), dimension(:),intent(inout) :: rgridLevels
  
  integer, intent(in) :: inlmin
  
  integer, intent(in) :: inlmax
!</inputoutput>
!</subroutine>
  integer :: i
  integer :: nlmin, nlmax
  nlmin=1
  nlmax=nlmin + (inlmax - inlmin)
  ! release on all levels
  do i=nlmin,nlmax
  
    ! free the vectors
    if (rgridLevels(i)%rvectorAreaBlockQ1%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rvectorAreaBlockQ1)
    end if
    
    if (rgridLevels(i)%rvectorMonFuncQ1%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rvectorMonFuncQ1)
    end if
    
    if (rgridLevels(i)%rSolBlock%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rSolBlock)
    end if
    
    if (rgridLevels(i)%rrhsBlock%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rrhsBlock)
    end if
    
    if (rgridLevels(i)%rtempBlock%h_Ddata .ne. ST_NOHANDLE) then
       call lsysbl_releaseVector(rgridLevels(i)%rtempBlock)
    end if
    
    if (rgridLevels(i)%rvecGradBlock%h_Ddata .ne. ST_NOHANDLE) then
       call lsysbl_releaseVector(rgridLevels(i)%rvecGradBlock)
    end if

    if (rgridLevels(i)%rmatDeform%NEQ .ne. 0) then
       call lsysbl_releaseMatrix(rgridLevels(i)%rmatDeform)
    end if

    if (rgridLevels(i)%rmatrix%NEQ .ne. 0) then
       call lsyssc_releaseMatrix(rgridLevels(i)%rmatrix)
    end if
  
    call spdiscr_releaseBlockDiscr(rgridLevels(i)%rdiscretisation)
  
  end do
  
  ! release the local working copy of the grid
  call tria_done (rgridLevels(nlmax)%rtriangulation)

  end subroutine  ! end griddef_freeLevels

!****************************************************************************************
!<subroutine>
  subroutine griddef_cleanLevels(rgridLevels,inlmin,inlmax)
!<description>
  !
  ! Zero out the information that we need for the next step
  !
!</description>
  
!<inputoutput>
  ! the grid levels
  type(t_hgridLevels), dimension(:),intent(inout) :: rgridLevels
  
  integer, intent(in) :: inlmin
  
  integer, intent(in) :: inlmax
!</inputoutput>
!</subroutine>
  integer :: i
  integer :: nlmin, nlmax
  nlmin=1
  nlmax=nlmin + (inlmax - inlmin)
  ! release on all levels
  do i=nlmin,nlmax
  
    ! free the vectors
    if (rgridLevels(i)%rvectorAreaBlockQ1%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rvectorAreaBlockQ1)
    end if
    
    if (rgridLevels(i)%rvectorMonFuncQ1%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rvectorMonFuncQ1)
    end if
    
    if (rgridLevels(i)%rSolBlock%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rSolBlock)
    end if
    
    if (rgridLevels(i)%rrhsBlock%h_Ddata .ne. ST_NOHANDLE) then
      call lsysbl_releaseVector(rgridLevels(i)%rrhsBlock)
    end if
    
    if (rgridLevels(i)%rtempBlock%h_Ddata .ne. ST_NOHANDLE) then
       call lsysbl_releaseVector(rgridLevels(i)%rtempBlock)
    end if
    
    if (rgridLevels(i)%rvecGradBlock%h_Ddata .ne. ST_NOHANDLE) then
       call lsysbl_releaseVector(rgridLevels(i)%rvecGradBlock)
    end if

    if (rgridLevels(i)%rmatDeform%NEQ .ne. 0) then
       call lsysbl_releaseMatrix(rgridLevels(i)%rmatDeform)
    end if

    if (rgridLevels(i)%rmatrix%NEQ .ne. 0) then
       call lsyssc_releaseMatrix(rgridLevels(i)%rmatrix)
    end if
  
  end do
  
  end subroutine  ! end griddef_freeLevels

!****************************************************************************************

!<subroutine>
  subroutine griddef_reinitDefInfo(rgriddefInfo)
!<description>
  !
  ! Here we release the rgriddefInfo structure
  ! and all the elements that are contained
  ! within it.
  !
!</description>
  
!<inputoutput>
  ! We will fill this structure with useful values
  type(t_griddefInfo), intent(inout) :: rgriddefInfo
!</inputoutput>
!</subroutine>

  ! deallocate memory
  if (rgriddefInfo%h_Dql2 .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_Dql2)
  end if
  
  ! deallocate memory
  if (rgriddefInfo%h_Dqlinfty .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_Dqlinfty)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_calcAdapSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_calcAdapSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nadapSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nadapSteps)
  end if


  ! deallocate memory
  if (rgriddefInfo%h_nadapStepsReally .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nadapStepsReally)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nmaxCorrSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nmaxCorrSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_ncorrSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_ncorrSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_ilevelODE .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_ilevelODE)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_ilevelODECorr .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_ilevelODECorr)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_npreSmoothSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_npreSmoothSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_npostSmoothSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_npostSmoothSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nintSmoothSteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nintSmoothSteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_nmonSmoothsteps .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_nmonSmoothsteps)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_cpreSmoothMthd .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_cpreSmoothMthd)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_cpostSmoothMthd .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_cpostSmoothMthd)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_cintSmoothMthd .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_cintSmoothMthd)
  end if

  ! deallocate memory
  if (rgriddefInfo%h_DblendPar .ne. ST_NOHANDLE) then
    call storage_free(rgriddefInfo%h_DblendPar)
  end if
  
  
  end subroutine
  
!****************************************************************************************
  
!<subroutine>
  subroutine griddef_calcMonitorFunction(rgriddefInfo,iLevel,def_monitorfct)
  !<description>
    !
    ! In this subroutine we calculate the values of the monitor
    ! function by means of the given monitor function.
    !
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo

    integer, intent(in) :: iLevel

  !</inputoutput>

    ! A callback routine for the monitor function
    include 'intf_monitorfct.inc'
    optional :: def_monitorfct

!</subroutine>

  ! local variables
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  real(dp), dimension(:), pointer :: p_Dentries
  type(t_vectorScalar), pointer :: p_rvectorMonFuncQ1
  type(t_blockDiscretisation), pointer :: rdiscretisation

!  if(rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%h_Ddata .ne. ST_NOHANDLE)then
!    return
!  end if
  
  rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation
  
  ! get the pointer to the coords of the grid
  ! this only works for Q1
  call storage_getbase_double2D (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DvertexCoords,&
    p_DvertexCoords)
    
  ! Set up an empty block vector
  call lsysbl_createVecBlockByDiscr(rdiscretisation,&
       rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1,.true.)
  
  ! Get a pointer just not to write such a long name
  p_rvectorMonFuncQ1 => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)
  
  ! get the data
  call storage_getbase_double(p_rvectorMonFuncQ1%h_ddata,p_Dentries)
  
  ! calculate monitor function
  call def_monitorfct(p_DvertexCoords,p_Dentries)

   
  end subroutine
  
!****************************************************************************************

!<subroutine>
  subroutine griddef_prepareDeformation(rgriddefInfo)
  !<description>
    ! This subroutine performs all preparations for deformation which are necessary
    ! for the first deformation run in the program and after every change of the compute
    ! level.
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
  !</inputoutput>

!</subroutine>

  ! local variables
  
  ! maximum number of correction steps (local instance)
  integer:: NEL, NVT, NLMAX

  ! get these numbers
  NLMAX = rgriddefInfo%iminDefLevel
  NEL   = rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%NEL
  NVT   = rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%NVT
                                        
  end subroutine ! end griddef_prepareDeformation

! ****************************************************************************************

!<subroutine>
  subroutine griddef_performDeformation(rgriddefInfo, &
                                        h_Dcontrib,&
                                        bstartNew, blevelHasChanged, bterminate, &
                                        bdegenerated, imgLevelCalc, iiteradapt, ibcIdx,&
                                        def_monitorfct,rperfconfig)
  !<description>
    ! This subroutine is the main routine for the grid deformation process, as all
    ! necessary steps are included here. For performing grid deformation, it is sufficient
    ! to define a monitor function and call this routine or to have an error distribution
    ! at hand and call this subroutine.
  !</description>

  !<input>
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.

    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
    
    ! if true, start from scratch: new vectors, new boundary conditions structures
    logical, intent(in) :: bstartNew

    ! if true, adjust the vectors after level change since a previous deformation call
    logical, intent(in) :: blevelHasChanged

    logical, intent(in) :: bterminate
    
    ! number of adaptive iteration
    integer, intent(in) :: iiterAdapt

    ! multigrid level on which the simulation was computed
    integer, intent(in) :: imgLevelCalc

    ! index of boundary condition related to the deformation PDE
    integer, intent(in):: ibcIdx

    ! flag for grid checking: if true, the deformation process would lead to
    ! a grid with tangled elements
    logical , intent(in):: bdegenerated

    ! A callback routine for the monitor function
    include 'intf_monitorfct.inc'
    optional :: def_monitorfct

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
    
  !</input>

  !<inoutput>
    ! handle of vector with elementwise error contributions
    integer, intent(inout):: h_Dcontrib
  !</inoutput>

!</subroutine>

  ! local variables
  integer :: nmaxCorrSteps, NLMAX,NLMIN,idef,i
  ! A scalar matrix and vector. The vector accepts the RHS of the problem
  ! in scalar form.
  integer, dimension(:), pointer :: p_nadapStepsReally
  real(dp), dimension(:), pointer :: p_DblendPar
  integer, dimension(:), pointer :: p_ilevelODE
  
  ! We deform on level nlmax
  NLMAX = rgriddefInfo%iminDefLevel
  NLMIN = rgriddefInfo%NLMIN
  ! get some needed pointers
  call storage_getbase_int(rgriddefInfo%h_nadapStepsReally,p_nadapStepsReally)
  
  call storage_getbase_int(rgriddefInfo%h_ilevelODE,p_ilevelODE)

  ! necessary for boundary projection only, not for adding the Neumann boundary
  ! condition

  ! create and modify all vectors necessary for deformation steps
  call griddef_prepareDeformation(rgriddefInfo)

  ! Compute the number of deformation steps
  call griddef_computeAdapSteps(rgriddefInfo)
  ! Compute the blending parameter
  
  ! compute number of deformation steps on current level ilevel
  !call griddef_computeAdapSteps(rparBlock, rgriddefInfo, ilevel, pfmon)

  nmaxCorrSteps = 0

  call storage_getbase_double(rgriddefInfo%h_DblendPar,p_DblendPar)

  do i=NLMIN,NLMAX
      call spdiscr_initBlockDiscr (rgriddefInfo%p_rhLevels(i)%rdiscretisation,1,&
           rgriddefInfo%p_rhLevels(i)%p_rtriangulation, rgriddefInfo%p_rboundary)
      
      call spdiscr_initDiscr_simple (rgriddefInfo%p_rhLevels(i)%rdiscretisation%RspatialDiscr(1),&
                                     EL_E011,CUB_G2X2,&
                                     rgriddefInfo%p_rhLevels(i)%p_rtriangulation,&
                                     rgriddefInfo%p_rboundary)
                                     
  end do

  ! loop over adaption cycles
  do idef = 1,rgriddefInfo%nadaptionSteps

    call output_lbrk ()
    print *,"Adaptation Step: ",idef
    call output_line ('-------------------')
    call output_lbrk ()

    ! perform one deformation step: This routine appplies deformation to the
    call  griddef_performOneDefStep(rgriddefInfo,&
                                    p_DblendPar(idef), NLMAX, NLMAX,&
                                    def_monitorfct,rperfconfig)
     
    ! nullifiy where neccesary
    call griddef_cleanLevels(rgriddefInfo%p_rhLevels,NLMIN,NLMAX)
  end do ! end do

  end subroutine
 
 ! ****************************************************************************************
 
!<subroutine>
  subroutine griddef_computeAdapSteps(rgriddefInfo)

  !<description>
    ! This subroutine performs all preparations for deformation which are necessary
    ! for the first deformation run in the program and after every change of the compute
    ! level.
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo

    ! this function should get the current level, do dedide, if
    ! the number of steps should be adaptively on this level.
  !</inputoutput>
    
!</subroutine>

    ! local variables
    integer :: idef
    
    integer,dimension(:), pointer :: p_calcAdapSteps
    
    real(dp), dimension(:), pointer :: p_DblendPar
    
    call storage_getbase_int(rgriddefInfo%h_calcAdapSteps,p_calcAdapSteps)
    
    ! we prescribe a fixed number of adaptation steps
    if(p_calcAdapSteps(rgriddefInfo%iminDefLevel) .eq. GRIDDEF_FIXED)then
    
      rgriddefInfo%nadaptionSteps = 20
    
    end if
    if(p_calcAdapSteps(rgriddefInfo%iminDefLevel) .eq. GRIDDEF_USER)then
      ! the Steps are not really set adaptively in
      ! case we enter this branch, the p_calcAdapSteps have already been
      ! set.
    end if
    
    
    ! now allocate memory for the parameters
    call storage_new('griddef_computeAdapSteps','rgriddefInfo%h_DblendPar',&
          rgriddefInfo%nadaptionSteps,ST_doUBLE,rgriddefInfo%h_DblendPar,&
          ST_NEWBLOCK_ZERO)
          
    call storage_getbase_double(rgriddefInfo%h_DblendPar,p_DblendPar)
    
    ! compute the blending parameter
    do idef=1,rgriddefInfo%nadaptionSteps-1
      ! evaluate formula
      p_DblendPar(idef) = &
      sqrt(sqrt(real(idef)/real(rgriddefInfo%nadaptionSteps)))
    end do ! end idef
    
    ! set the parameter for the final step
    p_DblendPar(rgriddefInfo%nadaptionSteps) = 1.0_dp
                                      
  end subroutine ! griddef_computeAdapSteps
  
! ****************************************************************************************

!<subroutine>
  subroutine griddef_performOneDefStep(rgriddefInfo,&
                                       dblendpar, ilevelODE, ilevel,&
                                       def_monitorfct,rperfconfig)
  !<description>
    ! This subroutine performs one deformation step of the enhanced deformation method.
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
  !</inputoutput>
  
  !<input>
    ! absolute level on which the points are moved
    integer, intent(in) :: ilevelODE

    ! absolute level on which the deformation PDE is solved
    integer, intent(in) :: ilevel
    
    ! blending parameter for monitorfunction sf + (1-s)g
    real(DP), intent(inout) :: dblendpar

    ! A callback routine for the monitor function
    include 'intf_monitorfct.inc'
    optional :: def_monitorfct

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig

  !</input>

!</subroutine>

    ! local variables
    real(dp), dimension(:,:), pointer :: Dresults

    ! weighting parameter for Laplacian smoothing of the vector field
    real(dp) :: dscale1,dscale2
    logical :: bBlending
    
    ! An array of problem levels for the multigrid solver
!    type(t_level), dimension(:), pointer :: Rlevels
    
    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rDubDiscretisation
    
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(2), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    ! Error indicator during initialisation of the solver
    integer :: ierror,NLMAX, i,NLMIN
    
    NLMIN=rgriddefInfo%NLMIN
    NLMAX=rgriddefInfo%NLMAX
    
    ! initialise Dresults
    Dresults  => null()

    ! no blending
    bBlending = .true.
    
    !---------------------------------------------------------------------------------
    !                                 FOR EVERY LEVEL
    !---------------------------------------------------------------------------------
    
    do i=NLMIN,NLMAX
    
      ! compute the original area distribution g(x)
      call griddef_getArea(rgriddefInfo,i)

      if(present(def_monitorfct))then
        ! calculate the monitor/target grid area distribution f(x)
        call griddef_calcMonitorFunction(rgriddefInfo,i,def_monitorfct)
      else
        ! or build a test monitor function
        call griddef_buildMonFuncTest(rgriddefInfo, i)
      end if
      

      ! blend monitor with current area distribution, if necessary
      if(bBlending) then
        ! Compute the scaling factors dscale1 and dscale2 for f and g and scale them
        call griddef_normaliseFctsNum(rgriddefInfo,dscale1,dscale2,i)
        call griddef_blendmonitor(rgriddefInfo,dblendpar,i)
      end if

      ! normalise the reciprocal of the functions
      call griddef_normaliseFctsInv(rgriddefInfo,i,rperfconfig)
    
    
      
      ! create the matrix for the poisson problem
      call griddef_createMatrixDef(rgriddefInfo,i)

      ! create rhs for deformation problem
      call griddef_createRhs(rgriddefInfo,i)
    
    end do
    
    !---------------------------------------------------------------------------------
    !     end                          FOR EVERY LEVEL
    !---------------------------------------------------------------------------------
    
    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFreal
    RfilterChain(2)%ifilterType = FILTER_SMALLL1TO0
    RfilterChain(2)%ismallL1to0component = 1
    
    
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
    
    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      
      ! Create an ILU(0) smoother
      call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
    end do
    
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (rgriddefInfo%p_rhLevels(i)%rmatDeform,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))
    
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    p_rsolverNode%nmaxIterations = 20
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rgriddefInfo%p_rhLevels(NLMAX)%rSolBlock,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%rrhsBlock,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%rtempBlock)

    !-----------------------------------------------------------------------------
    !                            MESH MOVING PHASE
    !-----------------------------------------------------------------------------

    call spdiscr_initBlockDiscr (rDubDiscretisation,2,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation)

    call spdiscr_deriveSimpleDiscrSc (&
                 rgriddefInfo%p_rhLevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                 EL_Q1, CUB_G2X2, rDubDiscretisation%RspatialDiscr(1))
                 
    call spdiscr_deriveSimpleDiscrSc (&
                 rgriddefInfo%p_rhLevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                 EL_Q1, CUB_G2X2, rDubDiscretisation%RspatialDiscr(2))
                 
                 
    ! initialise the block vector that should hold the solution
    call lsysbl_createVecBlockByDiscr(rDubDiscretisation,&
                                      rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock,.true.)
    
    ! get the recovered gradient of the solution
    call ppgrd_calcGradient(rgriddefInfo%p_rhLevels(NLMAX)%rSolBlock%RvectorBlock(1),&
                            rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock,PPGRD_INTERPOL)

    ! Solve the ODE and move the mesh
    call griddef_moveMesh(rgriddefInfo)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! we have a dirty work structure
    !rgriddefWork%breinit= .true.
    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    
    call spdiscr_releaseBlockDiscr(rDubDiscretisation)
    
  
  end subroutine ! end griddef_performOneDefStep

  ! ***************************************************************************

!<subroutine>
  subroutine griddef_getArea(rgriddefInfo,iLevel)
  
  !<description>
    ! In this function we build the nodewise area distribution out
    ! of an elementwise distribution
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo

    integer, intent(in) :: iLevel
    
  !</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Darea
    integer :: iel
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: Dpoints
    integer :: ive
    type(t_vectorScalar) :: rvectorAreaQ0
    type(t_vectorBlock) :: rvectorAreaBlockQ0
    type(t_blockDiscretisation) :: rprjDiscretisation
    type(t_blockDiscretisation), pointer :: rdiscretisation
    
    ! Is everything here we need?
    if (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line ('h_DvertexCoords not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      call sys_halt()
    end if

    if (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      call sys_halt()
    end if
    
    ! Do we have (enough) memory for that array?
    if (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DelementVolume .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementVolume2D', 'DAREA', &
          int(rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%NEL+1,I32), ST_doUBLE, &
          rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    end if
    
    rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation
    
    ! Get the arrays
    call storage_getbase_double2D (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    
    ! Set up an empty block vector
    call lsysbl_createVecBlockByDiscr(rdiscretisation,&
         rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1,.true.)
        
    ! Create a discretisation structure for Q0, based on our
    ! previous discretisation structure:
    call spdiscr_duplicateBlockDiscr(rdiscretisation,rprjDiscretisation)
    call spdiscr_deriveSimpleDiscrSc (&
                 rdiscretisation%RspatialDiscr(1), &
                 EL_Q0, CUB_G2X2, rprjDiscretisation%RspatialDiscr(1))
                 
    ! Initialise a Q0 vector from the newly created discretisation
    call lsyssc_createVecByDiscr(rprjDiscretisation%RspatialDiscr(1), &
    rvectorAreaQ0,.true.)
    
    ! get the pointer to the entries of this vector
    call lsyssc_getbase_double(rvectorAreaQ0,p_Darea)
    
    ! Loop over all elements calculate the area
    ! and save it in our vector
    do iel=1,rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%NEL
      
      if (p_IverticesAtElement(4,iel) .eq. 0) then
        ! triangular element
        do ive=1,TRIA_NVETRI2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        end do
        p_Darea(iel) = gaux_getArea_tria2D(Dpoints)
      else
        ! quad element
        do ive=1,TRIA_NVEQUAD2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        end do
        p_Darea(iel) = gaux_getArea_quad2D(Dpoints)
      end if

    end do ! end iel
    
    ! now transform the q0 vector into a q1 vector
    ! Setup a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecFromScalar(rvectorAreaQ0,rvectorAreaBlockQ0,rprjDiscretisation)
 
    ! Take the original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution(rvectorAreaBlockQ0,rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1)
    
    call lsysbl_releaseVector(rvectorAreaBlockQ0)
    call lsyssc_releaseVector(rvectorAreaQ0)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine griddef_getArea3D(rgriddefInfo,iLevel)
  
  !<description>
    ! In this function we build the nodewise area distribution out
    ! of an elementwise distribution
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo

    integer, intent(in) :: iLevel
    
  !</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Darea
    real(DP), dimension(:), pointer :: p_DelementVolume
    integer :: iel
    type(t_vectorScalar) :: rvectorAreaQ0
    type(t_vectorBlock) :: rvectorAreaBlockQ0
    type(t_blockDiscretisation) :: rprjDiscretisation
    type(t_blockDiscretisation), pointer :: rdiscretisation
    
    rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation
    
    ! Set up an empty block vector
    call lsysbl_createVecBlockByDiscr(rdiscretisation,&
         rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1,.true.)
        
    ! Create a discretisation structure for Q0, based on our
    ! previous discretisation structure:
    call spdiscr_duplicateBlockDiscr(rdiscretisation,rprjDiscretisation)
    call spdiscr_deriveSimpleDiscrSc (&
                 rdiscretisation%RspatialDiscr(1), &
                 EL_Q0_3D, CUB_G3_3D, rprjDiscretisation%RspatialDiscr(1))
                 
    ! Initialise a Q0 vector from the newly created discretisation
    call lsyssc_createVecByDiscr(rprjDiscretisation%RspatialDiscr(1), &
    rvectorAreaQ0,.true.)
    
    ! get the pointer to the entries of this vector
    call lsyssc_getbase_double(rvectorAreaQ0,p_Darea)
    
    call tria_genElementVolume3D(rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation)
    
    call storage_getbase_double (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DelementVolume,&
        p_DelementVolume)

    do iel=1,rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%NEL
    p_Darea(iel)=p_DelementVolume(iel)
    end do
    
    ! now transform the q0 vector into a q1 vector
    ! Setup a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecFromScalar(rvectorAreaQ0,rvectorAreaBlockQ0,rprjDiscretisation)
 
    ! Take the original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution(rvectorAreaBlockQ0,rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1)
    
    call lsysbl_releaseVector(rvectorAreaBlockQ0)
    call lsyssc_releaseVector(rvectorAreaQ0)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  subroutine griddef_buildMonFuncTest(rgriddefInfo,iLevel)
  
  !<description>
    ! In this function we build the nodewise area distribution out
    ! of an elementwise distribution
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
    
    integer, intent(in) :: iLevel
  !</inputoutput>

!</subroutine>

  ! local variables
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(DP), dimension(:), pointer :: p_Dentries
  integer :: ive
  type(t_vectorScalar), pointer :: p_rvectorMonFuncQ1
  integer :: iMethod
  real(DP) :: Dist
  type(t_blockDiscretisation), pointer :: rdiscretisation
  
  rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation
  
  iMethod = 1
  
  ! get the pointer to the coords of the grid
  ! this only works for Q1
  call storage_getbase_double2D (rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DvertexCoords,&
    p_DvertexCoords)
    
  ! Set up an empty block vector
  call lsysbl_createVecBlockByDiscr(rdiscretisation,&
       rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1,.true.)

  ! Get a pointer just not to write such a long name
  p_rvectorMonFuncQ1 => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)
  
  ! get the data
  call storage_getbase_double(p_rvectorMonFuncQ1%h_ddata,p_Dentries)
  
  select case(imethod)
    case(0)
      ! loop over all vertices and compute the monitor function
      do ive=1,rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%NVT
        p_Dentries(ive) = 0.5_dp + p_DvertexCoords(1,ive)
        !p_Dentries(ive) = 1.0_dp
      end do
    case(1)
      ! loop over all vertices and compute the monitor function
      do ive=1,rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%NVT
        Dist = sqrt((0.5_dp - p_DvertexCoords(1,ive))**2 + (0.5_dp - p_DvertexCoords(2,ive))**2)
        ! Good now define the monitor function
        Dist = abs(Dist - 0.2_dp)/0.2_dp
        Dist=max(dist,0.1_dp)
        Dist=min(1.0_dp,dist)
        p_Dentries(ive)=Dist
      end do
    case default
  end select
  
  end subroutine ! end griddef_buildMonFuncTest
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine griddef_normaliseFctsNum(rgriddefInfo,dScale1,dScale2,iLevel)
  !<description>
    ! We normalize the functions f and g so that
    ! <tex>$$ dScale1 * \int_\Omega f = dScale2 * \int_\Omega g = |\Omega|  $$</tex>
    !
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo

    integer, intent(in) :: iLevel
    
  !</inputoutput>
  
  !<output>
  real(dp),intent(inout) :: dScale1
  real(dp),intent(inout) :: dScale2
  !</output>

!</subroutine>

  ! local variables
  ! a shorthand to the functions
  type(t_vectorScalar) , pointer :: p_Df1
  type(t_vectorScalar) , pointer :: p_Df2
  real(DP), dimension(:), pointer :: p_Data1
  real(DP), dimension(:), pointer :: p_Data2
  ! These will be the values of the integral
  real(DP) :: dIntF1, dIntF2,Domega
  ! Element area
  real(DP), dimension(:), pointer :: p_DelementVolume
      
  ! initialise integral value with zero
  dIntF1 = 0.0_dp
  dIntF2 = 0.0_dp
  ! we do not want to write this name
  p_Df1 => rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1%RvectorBlock(1)
  ! we do not want to write this name
  p_Df2 => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)
  
  
  call storage_getbase_double(rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DelementVolume,&
                              p_DelementVolume)
  
  ! Integrate
  call pperr_scalar (p_Df1,PPERR_L1ERROR,dIntF1)
  call pperr_scalar (p_Df2,PPERR_L1ERROR,dIntF2)
  
  ! The omega value is the total area of the domain
  Domega = p_DelementVolume(rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%NEL+1)
  
  ! compute the scaling factors
  dScale1 = Domega/dIntF1
  dScale2 = Domega/dIntF2
  
  ! get the function data
  call lsyssc_getbase_double(p_Df1,p_Data1)
  call lsyssc_getbase_double(p_Df2,p_Data2)
  
  ! scale the functions
  p_Data1(:) =  p_Data1(:) *  dScale1
  p_Data2(:) =  p_Data2(:) *  dScale2
                                     
  end subroutine  ! end griddef_normaliseFctsNum
  
  ! ***************************************************************************

!<subroutine>
  subroutine griddef_normaliseFctsInv(rgriddefInfo,iLevel,rperfconfig)
!<description>
    !
    !
!</description>

!<input>
    integer, intent(in) :: iLevel

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo),intent(inout) :: rgriddefInfo
!</inputoutput>

!</subroutine>

  ! local variables
  ! a shorthand to the functions
  type(t_vectorScalar) , pointer :: p_Df1
  type(t_vectorScalar) , pointer :: p_Df2
  real(dp), dimension(:), pointer :: p_Data1
  real(dp), dimension(:), pointer :: p_Data2
  ! These will be the values of the integral
  real(dp) :: dIntF1, dIntF2,domega,dScale1,dScale2
  ! Element area
  real(dp), dimension(:), pointer :: p_DelementVolume
      
  type(t_blockDiscretisation),pointer :: rdiscretisation

  rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation
      
  ! initialise integral value with zero
  dIntF1 = 0.0_dp
  dIntF2 = 0.0_dp
  ! we do not want to write this name
  p_Df1 => rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1%RvectorBlock(1)
  ! we do not want to write this name
  p_Df2 => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)
  
  ! get the function data
  call lsyssc_getbase_double(p_Df1,p_Data1)
  call lsyssc_getbase_double(p_Df2,p_Data2)
  
  ! get the element area pointer
  call storage_getbase_double(&
  rgriddefInfo%p_rhLevels(iLevel)%p_rtriangulation%h_DelementVolume, &
                              p_DelementVolume)
  
  ! Integrate the functions f and g
  call griddef_normaliseFctsInvAux(rgriddefInfo, iLevel, dIntF1, dIntF2,&
                                   domega, rperfconfig)
  
  ! compute the scaling factor
  dScale1 = dintF1/domega
  dScale2 = dintF2/domega
  
  ! scale the functions
  p_Data1(:) =  p_Data1(:) *  dScale1
  p_Data2(:) =  p_Data2(:) *  dScale2
  
  end subroutine ! end griddef_normaliseFctsInv
  
  !****************************************************************************

!<subroutine>

  subroutine griddef_normaliseFctsInvAux(rgriddefInfo,iLevel,&
                                         dValue1,dValue2,dOm,rperfconfig)

!<description>
  !
  ! Auxiliary function: Here we perform the actual
  ! evaluation of the function on the elements.
  !
!</description>

!<input>
  integer, intent(in) :: iLevel

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo),intent(inout) :: rgriddefInfo
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(dp), intent(out) :: dValue1,dValue2,dOm
!</output>

!</subroutine>

  ! Pointer to the vector entries
  real(dp), dimension(:), pointer :: p_DdataMon,p_DdataArea

  ! Allocateable arrays for the values of the basis functions -
  ! for test space.
  real(dp), dimension(:,:,:,:), allocatable, target :: DbasTest
  
  ! Number of local degees of freedom for test functions
  integer :: indofTest
  
  type(t_vectorScalar), pointer :: rvectorArea,rvectorMon
  
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), pointer :: rvectorScalar

  type(t_blockDiscretisation), pointer :: rdiscretisation

  integer :: i,k,icurrentElementDistr, ICUBP, NVE
  integer :: IEL, IELmax, IELset, IdoFE
  real(dp) :: OM
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(el_maxnder) :: Bder
  
  ! Cubature point coordinates on the reference element
  real(dp), dimension(cub_maxcubp, ndim3d) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(dp), dimension(cub_maxcubp) :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(dp), dimension(:,:), allocatable :: p_DcubPtsRef

  ! Arrays for saving Jacobian determinants and matrices
  real(dp), dimension(:,:), pointer :: p_Ddetj
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution
  
  ! Number of elements in the current element distribution
  integer :: NEL

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_evalElementSet) :: rintSubset
  
  ! An allocateable array accepting the doF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTest

  ! Type of transformation from the reference to the real element
  integer :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer :: cevaluationTag
  
  real(dp) :: daux1,daux2

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => griddef_perfconfig
    end if

    rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.
  
    ! we do not want to write this name every time
    rvectorArea => rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1%RvectorBlock(1)
    rvectorMon  => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)
    rvectorScalar => rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1%RvectorBlock(1)

    call lsyssc_getbase_double (rvectorMon,p_DdataMon)
    call lsyssc_getbase_double (rvectorArea,p_DdataArea)
  !-------------------------------------------------------------


    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)
    
    dValue1 = 0.0_DP
    dValue2 = 0.0_DP
    dOm     = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,rdiscretisation%RspatialDiscr(1)%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RspatialDiscr(1)%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local doF`s for trial functions
      indofTest = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistribution%ccubTypeEval, ncubp, Dxi, Domega)
      
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
      
      ! Allocate arrays for the values of the test functions.
      allocate(DbasTest(indofTest,elem_getMaxDerivative(p_relementDistribution%celement),&
               ncubp,nelementsPerBlock))
      
      ! Allocate memory for the doF`s of all the elements.
      allocate(IdofsTest(indofTest,nelementsPerBlock))

      ! Initialisation of the element set.
      call elprep_init(rintSubset)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
               
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)

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
      
        ! Calculate the global doF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global doF`s of our p_rperfconfig%NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation%RspatialDiscr(1),&
                                     p_IelementList(IELset:IELmax),IdofsTest)
                                     
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (rintSubset,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
        p_Ddetj => rintSubset%p_Ddetj

        ! In the next loop, we do not have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Calculate the values of the basis functions.
        call elem_generic_sim2 (p_relementDistribution%celement, &
            rintSubset, Bder, DbasTest)
        
        do IEL=1,IELmax-IELset+1
        
          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp
          
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            
            do IdoFE = 1,indofTest
              daux1 = daux1 + p_DdataArea(IdofsTest(IdoFE,IEL))* &
              DbasTest(IdoFE,DER_FUNC,ICUBP,IEL)
              daux2 = daux2 + p_DdataMon(IdofsTest(IdoFE,IEL))* &
              DbasTest(IdoFE,DER_FUNC,ICUBP,IEL)
              dOm = dOm + OM * DbasTest(IdoFE,DER_FUNC,ICUBP,IEL)
            end do
            
            dValue1 = dValue1 + OM/daux1
            dValue2 = dValue2 + OM/daux2

          end do ! ICUBP

        end do ! IEL
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(rintSubset)

      deallocate(p_DcubPtsRef)

    end do ! icurrentElementDistr

  end subroutine ! griddef_normaliseFctsInvAux
  
  !****************************************************************************

  
!<subroutine>
  subroutine griddef_blendmonitor(rgriddefInfo,dBlendPar,iLevel)
  
  !<description>
    ! This subroutine performs the blending between the monitor function and the area
    ! distribution function during the grid deformation process. The blending performed
    ! is a linear combination of these two functions:
    ! <tex>$ f_{mon}  \rightarrow t * f_{mon} + (1-t)*darea$</tex>.
    ! If the deformation is too harsh, the resulting grid can become invalid due to
    ! numerical errors. In this case, one performs several deformations with the blended
    ! monitor function, such that every single step results in relatively mild deformation.
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo),intent(inout) :: rgriddefInfo

  !</inputoutput>
    real(DP),intent(inout)  :: dBlendPar
    
    integer, intent(in) :: iLevel
    
!</subroutine>

    ! local variables
    ! a shorthand to the functions
    integer :: i
    type(t_vectorScalar) , pointer :: p_Df1
    type(t_vectorScalar) , pointer :: p_Df2
    real(DP), dimension(:), pointer :: p_Data1
    real(DP), dimension(:), pointer :: p_Data2
   ! blending parameter
        
    ! we do not want to write this name
    p_Df1 => rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1%RvectorBlock(1)
    ! we do not want to write this name
    p_Df2 => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)

    ! if the blending parameter is 1, nothing has to be done
    if (dblendPar .eq. 1.0_DP) then
      return
    endif

    ! ensure that the blending parameter is between 0 and 1
    if (dblendPar .gt. 1.0_DP) dblendPar = 1.0_DP
    if (dblendPar .lt. 0.0_DP) dblendPar = 0.0_DP

    ! get the function data
    call lsyssc_getbase_double(p_Df1,p_Data1)
    call lsyssc_getbase_double(p_Df2,p_Data2)

    ! scale the functions
    ! p_Data1(:) =  p_Data1(:) *  (dblendPar - 1.0_dp)
    do i=1,ubound(p_Data2,1)
      p_Data2(i) = dblendPar * p_Data2(i) + (1.0_dp - dblendPar) * p_Data1(i)
    end do
    
  end subroutine  ! end griddef_blendmonitor

  ! ***************************************************************************

!<subroutine>
  subroutine griddef_createMatrixDef(rgriddefInfo,iLevel)
!<description>
  !
  ! Here we create the matrix for the deformation problem
  !
!</description>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo),intent(inout) :: rgriddefInfo

  integer, intent(in) :: iLevel
  
!</inputoutput>

!</subroutine>

  type(t_blockDiscretisation), pointer :: rdiscretisation

  rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation

  ! We create a scalar matrix, based on the discretisation structure
  ! for our one and only solution component.
  call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                   LSYSSC_MATRIX9,&
                                   rgriddefInfo%p_rhLevels(iLevel)%rmatrix)

  call stdop_assembleLaplaceMatrix(rgriddefInfo%p_rhLevels(iLevel)%rmatrix,.true.,1.0_dp)

  ! The linear solver only works for block matrices/vectors - so make the
  ! the matrix for the deformation problem a block matrix
  call lsysbl_createMatFromScalar (rgriddefInfo%p_rhLevels(iLevel)%rmatrix,&
                                   rgriddefInfo%p_rhLevels(iLevel)%rmatDeform,&
                                   rdiscretisation)
     
  call lsysbl_createVecBlockIndMat(rgriddefInfo%p_rhLevels(iLevel)%rmatDeform,&
                                   rgriddefInfo%p_rhLevels(iLevel)%rrhsBlock,.true.)
  
  ! Now we have block vectors for the RHS and the matrix. What we
  ! need additionally is a block vector for the solution and
  ! temporary data. Create them using the RHS as template.
  ! Fill the solution vector with 0:
  call lsysbl_createVecBlockIndirect (rgriddefInfo%p_rhLevels(iLevel)%rrhsBlock,&
                                      rgriddefInfo%p_rhLevels(iLevel)%rSolBlock, .true.)
  call lsysbl_createVecBlockIndirect (rgriddefInfo%p_rhLevels(iLevel)%rrhsBlock,&
                                      rgriddefInfo%p_rhLevels(iLevel)%rtempBlock, .true.)

  end subroutine

  !****************************************************************************

!<subroutine>
  subroutine griddef_createRhs (rgriddefInfo,iLevel,rperfconfig)
!<description>
  ! This routine calculates the entries of a discretised finite element vector.
  ! The discretisation is assumed to be conformal, i.e. the doF`s
  ! of all finite elements must 'match'.
  ! The linear form is defined by
  !        (f,$phi_i$), i=1..*
  ! with $Phi_i$ being the test functions defined in the discretisation
  ! structure.
  ! In case the array for the vector entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries
  ! and initialises all necessary variables of the vector according to the
  ! parameters (NEQ, pointer to the discretisation,...)
  !
  ! Double-precision version.
!</description>

!<input>
  integer, intent(in) :: iLevel

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
!  ! structure containing all vector handles for the deformation algorithm
!  type(t_griddefWork), intent(inout)  :: rgriddefWork

  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo),intent(inout) :: rgriddefInfo
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,k,icurrentElementDistr, ICUBP
  integer :: IEL, IELmax, IELset, IdoFE
  real(dp) :: OM
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(el_maxnder) :: Bder
  
  ! Cubature point coordinates on the reference element
  real(dp), dimension(cub_maxcubp, ndim3d) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(dp), dimension(cub_maxcubp) :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp
  
  ! Pointer to the vector entries
  real(dp), dimension(:), pointer :: p_DdataMon,p_DdataArea

  ! An allocateable array accepting the doF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTest
  integer, dimension(:,:), allocatable, target :: IdofsFunc
  
  ! Allocateable arrays for the values of the basis functions -
  ! for test space.
  real(dp), dimension(:,:,:,:), allocatable, target :: DbasTest,DbasFunc
  
  ! Number of entries in the vector - for quicker access
  integer :: NEQ
  
  ! Type of transformation from the reference to the real element
  integer :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer :: cevaluationTag

  ! Number of local degees of freedom for test functions
  integer :: indofTest,indofFunc
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
    
  ! An array that takes coordinates of the cubature formula on the reference element
  real(dp), dimension(:,:), allocatable :: p_DcubPtsRef

  ! Pointer to the jacobian determinants
  real(dp), dimension(:,:), pointer :: p_Ddetj
  
  ! Entries of the right hand side
  real(dp), dimension(:), pointer :: p_Ddata
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_elementDistribution
  type(t_elementDistribution), pointer :: p_elementDistributionFunc
  type(t_blockDiscretisation), pointer :: rdiscretisation
  
  ! Number of elements in the current element distribution
  integer :: NEL

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock
  
  ! Pointer to the coefficients that are computed by the callback routine.
  real(dp) :: dcoeff,dmonVal,dareaVal
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_evalElementSet) :: revalSubset
  logical :: bcubPtsInitialised
  type(t_vectorScalar), pointer :: rvectorArea,rvectorMon,rvectorRhs

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig
  
  if (present(rperfconfig)) then
    p_rperfconfig => rperfconfig
  else
    p_rperfconfig => griddef_perfconfig
  end if

  ! we do not want to write this name every time
  rvectorArea => rgriddefInfo%p_rhLevels(iLevel)%rvectorAreaBlockQ1%RvectorBlock(1)
  rvectorMon  => rgriddefInfo%p_rhLevels(iLevel)%rvectorMonFuncQ1%RvectorBlock(1)
  rdiscretisation => rgriddefInfo%p_rhLevels(iLevel)%rdiscretisation
  
  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDER
  ! according to these.
  Bder(:) = .false.
  Bder(DER_FUNC) = .true.
  
  ! Get information about the vector:
  NEQ = rvectorMon%NEQ

  rvectorRhs  => rgriddefInfo%p_rhLevels(iLevel)%rrhsBlock%RvectorBlock(1)
  call lsyssc_getbase_double(rvectorRhs,p_Ddata)
  call lsyssc_getbase_double(rvectorMon,p_DdataMon)
  call lsyssc_getbase_double(rvectorArea,p_DdataArea)
  

  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it is at most
  ! NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  !call ZTIME(DT(2))

  do icurrentElementDistr = 1,rdiscretisation%RspatialDiscr(1)%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => &
        rvectorRhs%p_rspatialDiscr%RelementDistr(icurrentElementDistr)
    p_elementDistributionFunc => &
        rvectorMon%p_rspatialDiscr%RelementDistr(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    if (p_elementDistribution%NEL .eq. 0) cycle

    ! Get the number of local doF`s for trial and test functions
    indofTest = elem_igetNDofLoc(p_elementDistribution%celement)
    indofFunc = elem_igetNDofLoc(p_elementDistributionFunc%celement)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_elementDistribution%celement)

    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    call cub_getCubPoints(p_elementDistribution%ccubTypeBilForm, ncubp, Dxi, Domega)
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,ncubp
      do k=1,ubound(p_DcubPtsRef,1)
        p_DcubPtsRef(k,i) = Dxi(i,k)
      end do
    end do
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%celement),&
             ncubp,nelementsPerBlock))
    allocate(DbasFunc(indofTest,elem_getMaxDerivative(p_elementDistributionFunc%celement),&
             ncubp,nelementsPerBlock))

    ! Allocate memory for the doF`s of all the elements.
    allocate(IdofsTest(indofTest,nelementsPerBlock))
    allocate(IdofsFunc(indofFunc,nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(revalSubset)

    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.
    
    !call ZTIME(DT(3))
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
                              
    ! Get the number of elements there.
    NEL = p_elementDistribution%NEL
  
    ! Loop over the elements - blockwise.
    do IELset = 1, NEL, p_rperfconfig%NELEMSIM
    
      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(NEL,IELset-1+p_rperfconfig%NELEMSIM)
    
      ! Calculate the global doF`s into IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global doF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rdiscretisation%RspatialDiscr(1),&
                                   p_IelementList(IELset:IELmax), &
                                   IdofsTest)
      call dof_locGlobMapping_mult(rdiscretisation%RspatialDiscr(1),&
                                   p_IelementList(IELset:IELmax), &
                                   IdofsFunc)
                                   
      !call ZTIME(DT(4))
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the doF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_elementDistribution%celement)
      cevaluationTag = ior(elem_getEvaluationTag(p_elementDistributionFunc%celement),&
      cevaluationTag)
                      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! Note: Why not using
      !   if (IELset .EQ. 1) then
      ! here, but this strange concept with the boolean variable?
      ! Because the if-command does not work with OpenMP! bcubPtsInitialised
      ! is a local variable and will therefore ensure that every thread
      ! is initialising its local set of cubature points!
      if (.not. bcubPtsInitialised) then
        bcubPtsInitialised = .true.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalSubset,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => revalSubset%p_Ddetj
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_elementDistribution%celement, &
          revalSubset, Bder, DbasTest)
      call elem_generic_sim2 (p_elementDistributionFunc%celement, &
          revalSubset, Bder, DbasFunc)

      ! --------------------- doF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the doF`s and cubature points to calculate the
      ! integral:
      
      do IEL=1,IELmax-IELset+1
        
        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp
        
          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
          
          ! Calculate 1/f-1/g in our cubature point
          dcoeff = 0.0_DP
          dmonVal  = 0.0_dp
          dareaVal = 0.0_dp
          do IdoFE=1,indofFunc
            dmonVal = dmonVal+DbasFunc(IdoFE,DER_FUNC,ICUBP,IEL)*p_DdataMon(IdofsFunc(IdoFE,IEL))
            dareaVal = dareaVal+DbasFunc(IdoFE,DER_FUNC,ICUBP,IEL)*p_DdataArea(IdofsFunc(IdoFE,IEL))
          end do ! IdoFE
          
          dmonVal = 1.0_dp/dmonVal - 1.0_dp/dareaVal
          ! Now loop through all possible combinations of doF`s
          ! in the current cubature point. Incorporate the data to the FEM vector

          do IdoFE=1,indofTest
            p_Ddata(IdofsTest(IdoFE,IEL)) = p_Ddata(IdofsTest(IdoFE,IEL)) + &
              DbasTest(IdoFE,DER_FUNC,ICUBP,IEL)*OM*dmonVal
          end do ! IdoFE
            
        end do ! ICUBP

      end do ! IEL

    end do ! IELset
    
    ! Release memory
    deallocate(IdofsTest,IdofsFunc)
    deallocate(DbasTest,DbasFunc)

    call elprep_releaseElementSet(revalSubset)
    
    deallocate(p_DcubPtsRef)

  end do ! icurrentElementDistr
  
  end subroutine ! griddef_createRHS

  !****************************************************************************
 
!<subroutine>
  subroutine griddef_moveMesh(rgriddefInfo)
!<description>
    ! This subroutine performs the actual deformation of the mesh. To do this, the grid
    ! points are moved in a vector field represented by DphiX and DphiY, the monitor
    ! function Dfmon and the area distribution Darea.
    ! The following ODE is solved:
    !<tex>
    ! \begin{displaymath}
    !   \frac{\partial \varphi(x,t)}{\partial t} = \eta(\varphi(x,t),t),
    !   \quad 0\leq t \leq  1, \quad
    !   \varphi(x,0) = x
    ! \end{displaymath}
    ! \begin{displaymath}
    !   \eta(y,s) := \frac{v(y)}{s\tilde{f}(y) + (1-s)\tilde{g}(y)}, \quad y \in \Omega,
    !   s \in [0,1].
    ! \end{displaymath}
    !</tex>
    ! To solve the ODE which describes the movement of the grid points, several ODE
    ! solvers are avialable. The type of solver is chosen by codeMethod. Currently, we
    ! support:
    !
    ! 1) explicit Euler (GRIDDEF_EXPL_EULER)
    !
    ! 6) RK4 (GRIDDEF_RK4)
    !

!</description>

!<input>
!</input>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout)  :: rgriddefInfo

!</inputoutput>

!</subroutine>

  ! local variables
  real(dp), dimension(:,:), pointer :: p_DvertexCoordsReal
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  real(dp), dimension(:), pointer :: p_DvertexParametersReal
  real(dp), dimension(:), pointer :: p_DvertexParameters
  integer :: i,NLMAX
  
  NLMAX=rgriddefInfo%NLMAX
  
  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_dvertexCoords,&
  p_DvertexCoords)

  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%h_dvertexCoords,&
  p_DvertexCoordsReal)
  
  call storage_getbase_double(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DvertexParameterValue,&
  p_DvertexParameters)
  
  call storage_getbase_double(rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%h_DvertexParameterValue,&
  p_DvertexParametersReal)
  
  call output_lbrk ()
  call output_line('Starting Explicit Euler..."')
  call output_lbrk ()
  
  ! call the explicit euler to move the grid points
  call griddef_performEE(rgriddefInfo)

  call output_line('Finished Explicit Euler..."')
  call output_lbrk ()
      
  ! write back coordinates
  do i=1,rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%NVT
    p_DvertexCoordsReal(1,i) = p_DvertexCoords(1,i)
    p_DvertexCoordsReal(2,i) = p_DvertexCoords(2,i)
  end do

  ! write back coordinates
  do i=1,rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%NVBD
    p_DvertexParametersReal(i) = p_DvertexParameters(i)
  end do
  
  end subroutine ! end griddef_moveMesh

  !****************************************************************************
 
!<subroutine>
  subroutine griddef_moveMesh3D(rgriddefInfo)
!<description>
    ! This subroutine performs the actual deformation of the mesh. To do this, the grid
    ! points are moved in a vector field represented by DphiX and DphiY, the monitor
    ! function Dfmon and the area distribution Darea.
    ! The following ODE is solved:
    !<tex>
    ! \begin{displaymath}
    !   \frac{\partial \varphi(x,t)}{\partial t} = \eta(\varphi(x,t),t),
    !   \quad 0\leq t \leq  1, \quad
    !   \varphi(x,0) = x
    ! \end{displaymath}
    ! \begin{displaymath}
    !   \eta(y,s) := \frac{v(y)}{s\tilde{f}(y) + (1-s)\tilde{g}(y)}, \quad y \in \Omega,
    !   s \in [0,1].
    ! \end{displaymath}
    !</tex>
    ! To solve the ODE which describes the movement of the grid points, several ODE
    ! solvers are avialable. The type of solver is chosen by codeMethod. Currently, we
    ! support:
    !
    ! 1) explicit Euler (GRIDDEF_EXPL_EULER)
    !
    ! 6) RK4 (GRIDDEF_RK4)
    !

!</description>

!<input>
!</input>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout)  :: rgriddefInfo

!</inputoutput>

!</subroutine>

  ! local variables
  real(dp), dimension(:,:), pointer :: p_DvertexCoordsReal
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
!   real(dp), dimension(:), pointer :: p_DvertexParametersReal
!   real(dp), dimension(:), pointer :: p_DvertexParameters
  integer :: NLMAX
  
  NLMAX=rgriddefInfo%NLMAX
  
  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_dvertexCoords,&
  p_DvertexCoords)

  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%h_dvertexCoords,&
  p_DvertexCoordsReal)
  
!  call storage_getbase_double(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DvertexParameterValue,&
!  p_DvertexParameters)
!
!  call storage_getbase_double(rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%h_DvertexParameterValue,&
!  p_DvertexParametersReal)
  
  call output_lbrk ()
  call output_line('Starting Explicit Euler..."')
  call output_lbrk ()
  
  ! call the explicit euler to move the grid points
  call griddef_performEE3D(rgriddefInfo)

  call output_line('Finished Explicit Euler..."')
  call output_lbrk ()
      
  ! write back coordinates
!  do i=1,rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%NVT
!    p_DvertexCoordsReal(1,i) = p_DvertexCoords(1,i)
!    p_DvertexCoordsReal(2,i) = p_DvertexCoords(2,i)
!    p_DvertexCoordsReal(3,i) = p_DvertexCoords(3,i)
!  end do

!  ! write back coordinates
!  do i=1,rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%NVBD
!    p_DvertexParametersReal(i) = p_DvertexParameters(i)
!  end do
  
  end subroutine ! end griddef_moveMesh3D
  
  !****************************************************************************

!<subroutine>
  subroutine griddef_performEE(rgriddefInfo)
!<description>
  !
  !
  !
!</description>

!<input>
!</input>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout)  :: rgriddefInfo
!</inputoutput>

!</subroutine>

  ! local variables
  
  ! true, if the point is not inside the domain
  logical :: bsearchFailed

  ! coordinates of evaluation point
  real(dp) :: dx, dy, dx_old, dy_old

  ! time and step size for ODE solving
  real(dp) :: dtime, dstepSize,deps

  ! level diference between PDE and PDE level
  integer:: ilevDiff

  ! number of ODE time steps
  integer::  ntimeSteps, ive, ielement,ivbd,nlmax

  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  integer, dimension(:), pointer :: p_IelementsAtVertex
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx
  integer, dimension(:), pointer :: p_InodalProperty

  real(dp), dimension(4) :: Dvalues
  real(dp), dimension(2) :: Dpoint

  deps = 0.0000000001_dp

  NLMAX=rgriddefInfo%NLMAX
  
  ! get the elements at vertex index array
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertexIdx,&
  p_IelementsAtVertexIdx)

  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_InodalProperty,&
  p_InodalProperty)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)
  
  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_dvertexCoords,&
  p_DvertexCoords)

  ntimeSteps = rgriddefInfo%ntimeSteps
  ! difference between grid level for ODE and PDE solving
  ! write the coordinates of the moved points to the actual coordinate vector
  ! Here, only the case ilevelODE < ilevel is considered, the opposite case
  ! is treated by prolongating the vector field to ODE level.
  ilevDiff = rgriddefInfo%imindefLevel

  do ive=1, rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%NVT
  
      dx_old = p_DvertexCoords(1,ive)
      dy_old = p_DvertexCoords(2,ive)
      
      dx = p_DvertexCoords(1,ive)
      dy = p_DvertexCoords(2,ive)
      
    ! if we have a boundary node we treat it in a special routine
    if(p_InodalProperty(ive) .ne. 0)then
      call griddef_perform_boundary2(rgriddefInfo,ive)
    else
      ! inner node
      ! initialise time variable
      dtime = 0.0_DP

      ! initialise flag for failed search
      bsearchFailed = .false.
      ! initial coordinates of the vertex
      dx = p_DvertexCoords(1,ive)
      dy = p_DvertexCoords(2,ive)
      
      dx_old = p_DvertexCoords(1,ive)
      dy_old = p_DvertexCoords(2,ive)
      
      
      ! here we store the coordinates
      Dpoint(1) = dx
      Dpoint(2) = dy
      
      ! zero the Dvalues
      Dvalues(:) = 0.0_dp
      
      ! for the first element read the element index
      ielement = p_IelementsAtVertex(p_IelementsAtVertexIdx(ive))
      
      ! evaluate the functions on the element
      call griddef_evalPhi_Known(DER_FUNC, Dvalues, &
           rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(1), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(2), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvectorMonFuncQ1%RvectorBlock(1), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvectorAreaBlockQ1%RvectorBlock(1), &
           Dpoint,ielement)

      ! compute step size for next time step
      dstepSize = 0.05_dp

      ! so this means in Dvalues(1) there is the
      ! x coordinate of the recovered gradient
      ! so this means in Dvalues(2) there is the
      ! y coordinate of the recovered gradient
      !((1.0_DP - dtime)*Dvalues(4) + &
      !dtime*(Dvalues(3)))
      ! In Dvalues(4) we find the g function(area distribution)
      ! In Dvalues(3) the f function (monitor)
      ! perform the actual Euler step
      dx = dx + dstepSize* Dvalues(1)/((1.0_DP - dtime)*Dvalues(4) + &
           dtime*(Dvalues(3)))
      dy = dy + dstepSize* Dvalues(2)/((1.0_DP - dtime)*Dvalues(4) + &
           dtime*(Dvalues(3)))

      ! update time
      dtime = dtime + dstepSize

      ! While the current time is less than 1.0-eps
      if (dtime .le. 1.0_DP - deps) then

         ! for the other time steps, we have really to search
        calculationloopEE_inner : do
        
          ! zero the Dvalues
          Dvalues(:) = 0.0_dp
          
          ! here we store the coordinates
          Dpoint(1) = dx
          Dpoint(2) = dy
          
          call griddef_evalphi_ray(DER_FUNC, Dvalues, &
             rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(1), &
             rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(2), &
             rgriddefInfo%p_rhLevels(NLMAX)%rvectorMonFuncQ1%RvectorBlock(1), &
             rgriddefInfo%p_rhLevels(NLMAX)%rvectorAreaBlockQ1%RvectorBlock(1), &
             Dpoint,bsearchFailed,ielement)

          ! if the point is outside the domain, stop moving it
          if (bsearchFailed) then
            bsearchFailed = .false.
            exit calculationloopEE_inner
          endif

          ! compute step size for next time step
          ! griddef_computeStepSize(dtime, ntimeSteps)
          dstepSize = 0.05_dp

          ! perform the actual Euler step
          dx = dx + dstepSize* Dvalues(1)/((1.0_DP - dtime)*Dvalues(4) + &
               dtime*(Dvalues(3)))
          dy = dy + dstepSize* Dvalues(2)/((1.0_DP - dtime)*Dvalues(4) + &
               dtime*(Dvalues(3)))

          ! advance in time
          dtime = dtime + dstepSize
          ivbd = 0
          
          ! if time interval greater 1.0_-eps, we are finished
          if (dtime .ge. 1.0_DP - deps) then
          
!                if(dx .gt. 1.0_dp)then
!                  print *,ive
!                end if
                p_DvertexCoords(1,ive) = dx
                p_DvertexCoords(2,ive) = dy
              exit calculationloopEE_inner
            end if ! (dtime .ge. 1.0_DP - deps)
        
        enddo calculationloopEE_inner
        
      ! in case time interval exhausted in the first time step
      else
        ! write the coordinates
              if(dx .gt. 1.0_dp)then
                print *,ive
              end if
        p_DvertexCoords(1,ive) = dx
        p_DvertexCoords(2,ive) = dy
      endif ! dtime
      
  end if ! nodal_property .ne. 0
  
  end do ! ive
  
  end subroutine ! end griddef_performEE

  !****************************************************************************

!<subroutine>
  subroutine griddef_performEE3D(rgriddefInfo)
!<description>
  !
  !
  !
!</description>

!<input>
!</input>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout)  :: rgriddefInfo
!</inputoutput>

!</subroutine>

  ! local variables
  
  ! true, if the point is not inside the domain
  logical :: bsearchFailed

  ! coordinates of evaluation point
  real(dp) :: dx, dy, dz, dx_old, dy_old, dz_old,dx_old1,dy_old1,dz_old1

  ! time and step size for ODE solving
  real(dp) :: dtime, dstepSize,deps

  ! level diference between PDE and PDE level
  integer:: ilevDiff

  ! number of ODE time steps
  integer::  ntimeSteps, ive, ielement,ivbd,nlmax

  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  integer, dimension(:), pointer :: p_IelementsAtVertex
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx
  integer, dimension(:), pointer :: p_InodalProperty

  real(dp), dimension(5) :: Dvalues
  real(dp), dimension(3) :: Dpoint

  deps = 0.0000000001_dp

  NLMAX=rgriddefInfo%NLMAX
  
  ! get the elements at vertex index array
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertexIdx,&
  p_IelementsAtVertexIdx)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)

  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_InodalProperty,&
  p_InodalProperty)
  
  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_dvertexCoords,&
  p_DvertexCoords)

  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_InodalProperty,&
  p_InodalProperty)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)
  
  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_dvertexCoords,&
  p_DvertexCoords)

  ntimeSteps = rgriddefInfo%ntimeSteps
  
  dstepSize = 1.0_dp/real(ntimeSteps,dp)

  ! difference between grid level for ODE and PDE solving
  ! write the coordinates of the moved points to the actual coordinate vector
  ! Here, only the case ilevelODE < ilevel is considered, the opposite case
  ! is treated by prolongating the vector field to ODE level.
  ilevDiff = rgriddefInfo%imindefLevel

  do ive=1, rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%NVT
     
    if(p_InodalProperty(ive) .ne. 0)cycle
    ! inner node
    ! initialise time variable
    dtime = 0.0_DP

    ! initialise flag for failed search
    bsearchFailed = .false.
    ! initial coordinates of the vertex
    dx_old1 = p_DvertexCoords(1,ive)
    dy_old1 = p_DvertexCoords(2,ive)
    dz_old1 = p_DvertexCoords(3,ive)
    
    dx_old = p_DvertexCoords(1,ive)
    dy_old = p_DvertexCoords(2,ive)
    dz_old = p_DvertexCoords(3,ive)
    
    dx = p_DvertexCoords(1,ive)
    dy = p_DvertexCoords(2,ive)
    dz = p_DvertexCoords(3,ive)
    
    
    ! some routines need the coordinates in
    ! this format
    Dpoint(1) = dx
    Dpoint(2) = dy
    Dpoint(3) = dz
    
    ! zero the Dvalues
    Dvalues(:) = 0.0_dp
    
    ! for the first element read the element index
    ielement = p_IelementsAtVertex(p_IelementsAtVertexIdx(ive))
    
    ! evaluate the functions on the element
    call griddef_evalPhi_Known3D(DER_FUNC, Dvalues, &
         rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(1), &
         rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(2), &
         rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(3), &
         rgriddefInfo%p_rhLevels(NLMAX)%rvectorMonFuncQ1%RvectorBlock(1), &
         rgriddefInfo%p_rhLevels(NLMAX)%rvectorAreaBlockQ1%RvectorBlock(1), &
         Dpoint,ielement)

    ! compute step size for next time step
   !dstepSize = 0.05_dp

    ! so this means in Dvalues(1) there is the
    ! x coordinate of the recovered gradient
    ! so this means in Dvalues(2) there is the
    ! y coordinate of the recovered gradient
    ! so this means in Dvalues(3) there is the
    ! z coordinate of the recovered gradient
    ! In Dvalues(5) we find the g function(area distribution)
    ! In Dvalues(4) the f function (monitor)
    ! perform the actual Euler step
    dx = dx + dstepSize* Dvalues(1)/((1.0_DP - dtime)*Dvalues(5) + &
         dtime*(Dvalues(4)))
    dy = dy + dstepSize* Dvalues(2)/((1.0_DP - dtime)*Dvalues(5) + &
         dtime*(Dvalues(4)))
    dz = dz + dstepSize* Dvalues(3)/((1.0_DP - dtime)*Dvalues(5) + &
         dtime*(Dvalues(4)))
         
    ! update time
    dtime = dtime + dstepSize

    ! While the current time is less than 1.0-eps
    if (dtime .le. 1.0_DP - deps) then

       ! for the other time steps, we have really to search
      calculationloopEE_inner : do
      
        ! zero the Dvalues
        Dvalues(:) = 0.0_dp
        
        ! here we store the coordinates
        Dpoint(1) = dx
        Dpoint(2) = dy
        Dpoint(3) = dz
             
        ! we search for the point, we enter the element it
        ! was found in last time
        call griddef_evalphi_ray3D(DER_FUNC, Dvalues, &
           rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(1), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(2), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(3), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvectorMonFuncQ1%RvectorBlock(1), &
           rgriddefInfo%p_rhLevels(NLMAX)%rvectorAreaBlockQ1%RvectorBlock(1), &
           Dpoint,bsearchFailed,ielement)

        ! if the point is outside the domain, stop moving it
        if (bsearchFailed) then
          !print *,Dpoint(:)
          bsearchFailed = .false.
          exit calculationloopEE_inner
        endif

        ! compute step size for next time step
        ! griddef_computeStepSize(dtime, ntimeSteps)
        !dstepSize = 0.01_dp

        ! perform the actual Euler step
        dx = dx + dstepSize* Dvalues(1)/((1.0_DP - dtime)*Dvalues(5) + &
             dtime*(Dvalues(4)))
        dy = dy + dstepSize* Dvalues(2)/((1.0_DP - dtime)*Dvalues(5) + &
             dtime*(Dvalues(4)))
        dz = dz + dstepSize* Dvalues(3)/((1.0_DP - dtime)*Dvalues(5) + &
             dtime*(Dvalues(4)))

        ! advance in time
        dtime = dtime + dstepSize
        ivbd = 0
        
        ! if time interval greater 1.0_-eps, we are finished
        if (dtime .ge. 1.0_DP - deps) then
        
!                if(dx .gt. 1.0_dp)then
!                  print *,ive
!                end if
              !print *,"old coordinates: ",dx_old1,dy_old1,dz_old1
              !print *,"new coordinates: ",dx,dy,dz
              
              p_DvertexCoords(1,ive) = dx
              p_DvertexCoords(2,ive) = dy
              p_DvertexCoords(3,ive) = dz
            exit calculationloopEE_inner
          end if ! (dtime .ge. 1.0_DP - deps)
      
      enddo calculationloopEE_inner
      
    ! in case time interval exhausted in the first time step
    else
      ! write the coordinates
            if(dx .gt. 1.0_dp)then
              print *,ive
            end if
      !print *,"old coordinates: ",dx_old1,dy_old1,dz_old1
      !print *,"new coordinates: ",dx,dy,dz
            
      p_DvertexCoords(1,ive) = dx
      p_DvertexCoords(2,ive) = dy
      p_DvertexCoords(3,ive) = dz
    endif ! dtime
  
  end do ! ive
  
  ! now take care of the boundary
  call griddef_project2Boundary(rgriddefInfo)
  
  end subroutine ! end griddef_performEE
  
  !****************************************************************************

!<subroutine>
subroutine griddef_perform_boundary2(rgriddefInfo,ive)
!<description>
  !
  !
  !
!</description>

!<input>
  ! Global vertex index
  integer:: ive
!</input>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout)  :: rgriddefInfo

!</inputoutput>

!</subroutine>

  ! local variables
  
  ! true, if the point is not inside the domain
  logical :: bsearchFailed

  ! coordinates of evaluation point
  real(DP) :: dx, dy, dparam1,dparam2,dparam,dalpha_01,dnx,dny,dtmp

  ! time and step size for ODE solving
  real(DP) :: dtime, dstepSize,deps,dalpha,dalpha_old,dalpha_start

  ! level diference between PDE and PDE level
  integer:: ilevDiff

  ! number of ODE time steps
  integer::  ntimeSteps,ielement,icp1,icp2,ivbd,ibd,iboundary

  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  integer, dimension(:), pointer :: p_IelementsAtVertex
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx
  integer, dimension(:), pointer :: p_InodalProperty

  ! these arrays are needed when we treat boundary vertices
  integer, dimension(:), pointer :: p_IboundaryCpIdx
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  real(dp), dimension(:), pointer :: p_DvertexParameterValue
  real(dp), dimension(:), pointer :: p_DvertexParameterValueNew
  
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  integer, dimension(:,:), pointer :: p_IelementsAtEdge
  integer, dimension(:), pointer ::   p_IedgesAtBoundary
  
  integer, dimension(4) :: Ivalues
  real(dp), dimension(4) :: Dvalues
  real(dp), dimension(2) :: Dpoint
  
  ! make all the regions
  type(t_boundaryRegion), dimension(:),pointer :: p_rregion
  
  integer, dimension(:), allocatable :: rElements
  
  ! integer
  integer :: iregions,iedge,iupper,ivertex,icount,iinelement,iend,NLMAX

  deps = 0.0000000001_dp
  
  NLMAX=rgriddefInfo%NLMAX

  ! Get the boundary information we need
  call storage_getbase_double(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DvertexParameterValue,&
  p_DvertexParameterValueNew)
  
  call storage_getbase_double(rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%h_DvertexParameterValue,&
  p_DvertexParameterValue)
  

  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IboundaryCpIdx,&
  p_IboundaryCpIdx)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtBoundary,&
  p_IverticesAtBoundary)

  ! get the elements at vertex index array
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertexIdx,&
  p_IelementsAtVertexIdx)

  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_InodalProperty,&
  p_InodalProperty)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)
  
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
  p_IelementsAtVertex)
  
  call storage_getbase_double2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_dvertexCoords,&
  p_DvertexCoords)
  
  call storage_getbase_int(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IedgesAtBoundary,&
  p_IedgesAtBoundary)
  
  call storage_getbase_int2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtEdge,&
  p_IverticesAtEdge)
  
  call storage_getbase_int2d(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtEdge,&
  p_IelementsAtEdge)
  
  Ivalues = (/1,2,3,4/)

  ! allocate the regions
  iupper = ubound(p_IedgesAtBoundary,1)
  ! ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
  allocate(p_rregion(boundary_igetNsegments &
  (rgriddefInfo%p_rboundary,p_InodalProperty(ive))))

  allocate(rElements(iupper))

  ntimeSteps = rgriddefInfo%ntimeSteps
  ! difference between grid level for ODE and PDE solving
  ! write the coordinates of the moved points to the actual coordinate vector
  ! Here, only the case ilevelODE < ilevel is considered, the opposite case
  ! is treated by prolongating the vector field to ODE level.
  ilevDiff = rgriddefInfo%imindefLevel

  do ivbd=1,4
    if(ive .eq. Ivalues(ivbd))then
    return
    end if
  end do

  iboundary = 0
  
  ! get parameter value
  ! get the index of ive in the p_IverticesAtBoundary array
  icp1 = p_IboundaryCpIdx(p_InodalProperty(ive))
  icp2 = p_IboundaryCpIdx(p_InodalProperty(ive)+1)-1
  do ibd=icp1,icp2
    iboundary = iboundary + 1
    if(p_IverticesAtBoundary(ibd).eq.ive)then
    exit
    end if
  end do
  
  ! get the number of boundary segments in
  ! the boundary component
  iend = boundary_igetNsegments(rgriddefInfo%p_rboundary,p_InodalProperty(ive))
  
  ! create the boundary regions
  do iregions=1,iend
  ! Idea: create the regions, check in which region the parameter is
    call boundary_createRegion(rgriddefInfo%p_rboundary, p_InodalProperty(ive),&
                             iregions,p_rregion(iregions))
  end do

  ! with this information we can assign the
  ! parameter value
  dalpha = p_DvertexParameterValue(iboundary)
  
  ! this will be the start alpha
  ! we also have a copy of the 01 parameterisation
  dalpha_start = dalpha
  
  ! convert this parameter in to length parameterisation
  dalpha = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                 p_InodalProperty(ive), dalpha,&
                                 BDR_PAR_01, BDR_PAR_LENGTH)

  ! we want to stay in this boundary region
  do iregions=1,iend
  ! Idea: create the regions, check in which region the parameter is
    if(boundary_isInRegion (p_rregion(iregions),p_InodalProperty(ive),dalpha_start))then
      exit
    end if
  end do
      
  ! inner node
  ! initialise time variable
  dtime = 0.0_DP

  ! initialise flag for failed search
  bsearchFailed = .true.

  ! initial coordinates of the vertex
  dx = p_DvertexCoords(1,ive)
  dy = p_DvertexCoords(2,ive)
  
  ! here we store the coordinates
  Dpoint(1) = dx
  Dpoint(2) = dy
  
  ! zero the Dvalues
  Dvalues(:) = 0.0_dp
  
  ! for the first element read the element index
  ielement = p_IelementsAtVertex(p_IelementsAtVertexIdx(ive))
  
  ! so this means in Dvalues(1) there is the
  ! x coordinate of the recovered gradient
  ! so this means in Dvalues(2) there is the
  ! y coordinate of the recovered gradient
  !((1.0_DP - dtime)*Dvalues(4) + &
  !dtime*(Dvalues(3)))
  ! In Dvalues(4) we find the g function(area distribution)
  ! In Dvalues(3) the f function (monitor)
  ! evaluate the functions on the element
  call griddef_evalPhi_Known(DER_FUNC, Dvalues, &
       rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(1), &
       rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(2), &
       rgriddefInfo%p_rhLevels(NLMAX)%rvectorMonFuncQ1%RvectorBlock(1), &
       rgriddefInfo%p_rhLevels(NLMAX)%rvectorAreaBlockQ1%RvectorBlock(1), &
       Dpoint,Ielement)

  ! compute step size for next time step
  dstepSize = 0.05_dp

  
  ! perform the actual Euler step
  ! get the normal vector
  call boundary_getNormalVec2D(rgriddefInfo%p_rboundary,p_InodalProperty(ive),&
                             dalpha,dnx, dny, BDR_NORMAL_MEAN, BDR_PAR_LENGTH)
  ! get the tangential vector
  dtmp = dnx
  dnx  = -dny
  dny  = dtmp
  ! project the recovered gradient(x,y) on (dnx,dny)
  dtmp = Dvalues(1) * dnx + Dvalues(2) * dny
                               
  ! project the speed vector on the tangential vector
  dalpha = dalpha + dstepSize * dtmp/((1.0_DP - dtime)*Dvalues(4) + &
       dtime*(Dvalues(3)))

  ! update time
  dtime = dtime + dstepSize

  ! While we are still in the [0,1] interval(with a small tolerance)
  if (dtime .le. 1.0_DP - deps) then

    !--------------------IHATETHIStypeOFLOOP----------------------
     ! for the other time steps, we have really to search
    calculationloopEE_bdy : do
    
    ! the alpha at the beginning of the loop in lngthParam
    dalpha_old = dalpha
    
    ! zero the Dvalues
    Dvalues(:) = 0.0_dp
    
    bsearchFailed = .true.
    icount = 0

    ! we need a 01 parameterisation for this search
    dalpha_01 = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                     p_InodalProperty(ive), dalpha,&
                                     BDR_PAR_LENGTH,BDR_PAR_01)
    
    ! loop over the edges
    do iedge=1,iupper
      ivertex = p_IedgesAtBoundary(iedge)
      ivertex = p_IverticesAtEdge(1,ivertex)
      ivbd = 0
      ! get parameter value
      icp1 = p_IboundaryCpIdx(p_InodalProperty(ivertex))
      icp2 = p_IboundaryCpIdx(p_InodalProperty(ivertex)+1)-1
      do ibd=icp1,icp2
        ivbd = ivbd + 1
        if(p_IverticesAtBoundary(ibd).eq.ivertex)then
        exit
        end if
      end do
      dparam = p_DvertexParameterValue(ivbd)
      
      if(.not.(boundary_isInRegion(p_rregion(iregions),p_InodalProperty(ive),dparam)))then
        !add the element add the edge
        cycle
      end if
      
      ! Get the second vertex at this edge
      ivertex = p_IedgesAtBoundary(iedge)
      ivertex = p_IverticesAtEdge(2,ivertex)

      ivbd = 0
      ! get parameter value
      icp1 = p_IboundaryCpIdx(p_InodalProperty(ivertex))
      icp2 = p_IboundaryCpIdx(p_InodalProperty(ivertex)+1)-1
      do ibd=icp1,icp2
        ivbd = ivbd + 1
        if(p_IverticesAtBoundary(ibd).eq.ivertex)then
        exit
        end if
      end do
      dparam2 = p_DvertexParameterValue(ivbd)
      
      if(dparam2 .gt. dparam)then
        dparam1 = dparam
      else
        dparam1 = dparam2
        dparam2 = dparam
      end if
      
      !
      ! Find the boundary edge that corresponds to
      ! the current parameter
      if((dalpha_01 .ge. dparam1) .and. (dalpha_01 .le. dparam2))then
        iinelement = p_IelementsAtEdge(1,p_IedgesAtBoundary(iedge))
        bsearchFailed = .false.
        
        exit
      end if

    end do ! end do (for all boundary edges)

 !----------------------------------We either found the edge or not---------------------------------------

    ! time interval exhausted, calculation finished and exit
    if (bsearchFailed) then
    
      ! convert and write
      dalpha_old = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                     p_InodalProperty(ive), dalpha_old,&
                                     BDR_PAR_LENGTH,BDR_PAR_01)
      
      ! get the x,y coordinates for the current parameter value
      call boundary_getCoords(rgriddefInfo%p_rboundary, p_InodalProperty(ive),&
                              dalpha_old, dx, dy)
    
      ! write
      p_DvertexCoords(1,ive) = dx
      p_DvertexCoords(2,ive) = dy
      exit calculationloopEE_bdy
    end if ! (dtime .ge. 1.0_DP - deps)

    ! evaluate phi now in element iinelement
    ! convert the parameter value and evaluate
    dalpha_01 = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                     p_InodalProperty(ive), dalpha,&
                                     BDR_PAR_LENGTH,BDR_PAR_01)

    ! evaluate the parameter value and get the coordinates
    call boundary_getCoords(rgriddefInfo%p_rboundary, p_InodalProperty(ive),&
                            dalpha_01, dx, dy)

    
    ! get the x,y coordinates from the parameter value
    Dpoint(1) = dx
    Dpoint(2) = dy
   
    ! search in which element the point is
    ! iel ivbd param value
    call griddef_evalphi_ray_bound(DER_FUNC, Dvalues, &
       rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(1), &
       rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock%RvectorBlock(2), &
       rgriddefInfo%p_rhLevels(NLMAX)%rvectorMonFuncQ1%RvectorBlock(1), &
       rgriddefInfo%p_rhLevels(NLMAX)%rvectorAreaBlockQ1%RvectorBlock(1), &
       Dpoint,iinelement)


    ! compute step size for next time step
    dstepSize = 0.05_dp !griddef_computeStepSize(dtime, ntimeSteps)
    
    ! Get the normal vector on the boundary for the current parameter
    ! value.
    call boundary_getNormalVec2D(rgriddefInfo%p_rboundary,p_InodalProperty(ive),&
                               dalpha_01,dnx, dny, BDR_NORMAL_MEAN, BDR_PAR_01)
                               
    ! get the tangential vector
    dtmp = dnx
    dnx  = -dny
    dny  = dtmp
    ! project the recovered gradient(x,y) on (dnx,dny)
    dtmp = Dvalues(1) * dnx + Dvalues(2) * dny

    ! perform the actual Euler step
    dalpha = dalpha + dstepSize* dtmp/((1.0_DP - dtime)*Dvalues(4) + &
       dtime*(Dvalues(3)))

    ! advance in time
    dtime = dtime + dstepSize
    ivbd = 0
    
    ! update the 01 param value
    dalpha_01 = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                     p_InodalProperty(ive), dalpha,&
                                     BDR_PAR_LENGTH,BDR_PAR_01)
    
    ! if the point is outside the domain, stop moving it
    if ((dalpha_01 .lt. p_rregion(iregions)%dminParam) .or. &
        (dalpha_01 .gt. p_rregion(iregions)%dmaxParam)) then
        
        ! convert and write
        dalpha_old = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                       p_InodalProperty(ive), dalpha_old,&
                                       BDR_PAR_LENGTH,BDR_PAR_01)
         
        ! get the x,y coordinates of the current parameter value
        call boundary_getCoords(rgriddefInfo%p_rboundary, p_InodalProperty(ive),&
                                dalpha_old, dx, dy)
        
!        ! write back the parameter value of the vertex
        p_DvertexParameterValueNew(iboundary) = dalpha_old
                                
        ! write the coordindates
        p_DvertexCoords(1,ive) = dx
        p_DvertexCoords(2,ive) = dy
      exit calculationloopEE_bdy
    endif
    
    ! time interval exhausted, calculation finished and exit
    if (dtime .ge. 1.0_DP - deps) then

      ! convert from length parameterisation to 01
      dalpha = boundary_convertParameter(rgriddefInfo%p_rboundary, &
                                     p_InodalProperty(ive), dalpha,&
                                     BDR_PAR_LENGTH,BDR_PAR_01)
    
      ! get the x,y coordinates of the current parameter value
      call boundary_getCoords(rgriddefInfo%p_rboundary, p_InodalProperty(ive),&
                              dalpha, dx, dy)
      
      p_DvertexParameterValueNew(iboundary) = dalpha
      ! write the coordindates
      p_DvertexCoords(1,ive) = dx
      p_DvertexCoords(2,ive) = dy
      exit calculationloopEE_bdy
    end if ! (dtime .ge. 1.0_DP - deps)
    
    enddo calculationloopEE_bdy
  ! in case time interval exhausted in the first time step
  else
    ! write the coordinates
    p_DvertexCoords(1,ive) = dx
    p_DvertexCoords(2,ive) = dy
  endif ! dtime

  deallocate(p_rregion)

  end subroutine
 
 !****************************************************************************

!<subroutine>
  subroutine griddef_evalPhi_Known(iderType, Dvalues, rvecGradX, &
             rvecGradY,rvecMon,rvecArea,Dpoint, &
             Ielements, IelementsHint, cnonmeshPoints)
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvecMon in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! type of function value to evaluate. One of the DER_xxxx constants
  integer, intent(in)           :: idertype
  
  ! the scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in)              :: rvecMon
  
  type(t_vectorScalar), intent(in)              :: rvecGradX
  
  type(t_vectorScalar), intent(in)              :: rvecGradY
  
  type(t_vectorScalar), intent(in)              :: rvecArea
  
  ! a list of points where to evaluate
  ! dimension(1..ndim,1..npoints)
  real(dp), dimension(:), intent(in) :: Dpoint
  
  ! optional: a list of elements containing the points Dpoints.
  ! if this is not specified the elment numbers containing the points
  ! are determined automatically
  integer, intent(in), optional :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements
  ! containing the points. This is ignored if Ielements is specified!
  integer, intent(in), optional :: IelementsHint
  
  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
  ! if a point is located outside of the domain. May happen e.g. in
  ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default
  ! parameter if cnonmeshPoints is not specified.
  integer, intent(in), optional :: cnonmeshPoints
  
!</input>


!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(dp), dimension(:), intent(out) :: Dvalues
  

  
!</output>

  ! local variables
  integer :: cnonmesh
  integer :: ieltype,indof,nve,ibas
  integer :: iel
  integer, dimension(:), pointer :: p_IelementDistr
  logical, dimension(el_maxnder) :: Bder
  
  real(dp), dimension(:), pointer :: p_ddatamon
  real(dp), dimension(:), pointer :: p_ddataarea
  real(dp), dimension(:), pointer :: p_ddatagradx
  real(dp), dimension(:), pointer :: p_ddatagrady
  
  ! Transformation
  integer :: ctrafotype
  real(dp), dimension(trafo_maxdimrefcoord) :: dparpoint
  
  ! Values of basis functions and doF`s
  real(DP), dimension(el_maxnbas,el_maxnder) :: Dbas
  integer, dimension(el_maxnbas) :: Idofs
  
  ! List of element distributions in the discretisation structure
  type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

  ! Evaluation structure and tag
  type(t_evalElement) :: revalElement
  integer :: cevaluationTag
  
  p_RelementDistribution => rvecMon%p_rspatialDiscr%RelementDistr
  
  ! for uniform discretisations, we get the element type in advance...
  if(rvecMon%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNifORM) then
  
    ! Element type
    ieltype = rvecMon%p_rspatialDiscr%RelementDistr(1)%celement
    
    ! get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! number of vertices on the element
    nve = elem_igetNDofLoc(ieltype)
    
    ! type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! Get the element evaluation tag; neccessary for preparation
    cevaluationTag = elem_getEvaluationTag(ieltype)
    
    nullify(p_IelementDistr)
  else
     call storage_getbase_int (&
          rvecMon%p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
  end if
  
  ! get the data vector
  select case (rvecMon%cdataType)
  case (ST_doUBLE)
    call lsyssc_getbase_double(rvecMon,p_DdataMon)
    call lsyssc_getbase_double(rvecArea,p_DdataArea)
    call lsyssc_getbase_double(rvecGradX,p_DdataGradX)
    call lsyssc_getbase_double(rvecGradY,p_DdataGradY)
  case (ST_SINGLE)
  case default
    call output_line ('Unsupported vector precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
    call sys_halt()
  end select

  ! What to evaluate
  Bder = .false.
  Bder(iderType) = .true.
  
  cnonmesh = FEVL_NONMESHPTS_NONE
  if (present(cnonmeshPoints)) cnonmesh = cnonmeshPoints
  
  ! We loop over all points
  iel = 1
  
  ! Get the element number that contains the point
  if(present(Ielements))then
    ! we have it...
    iel = Ielements
  end if
  
  ! get the type of the element iel
  if(associated(p_ielementdistr))then
    ieltype = p_RelementDistribution(p_IelementDistr(iel))%celement
    
    ! Get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(ieltype)
    
    ! Type of transformation from /to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(ieltype)
  end if
  
  ! Calculate the global doF`s on that element into IdofsTest.
  call dof_locGlobMapping(rvecMon%p_rspatialDiscr, iel, Idofs)
  
  ! Get the element shape information
  call elprep_prepareForEvaluation(revalElement, EL_EVLTAG_COORDS, &
       rvecMon%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)
       
  ! calculate the transformation of the point to the reference element
  call trafo_calcRefCoords(ctrafoType, revalElement%Dcoords, &
       Dpoint(:), DparPoint)
       
  ! Now calculate everything else what is necessary for the element
  call elprep_prepareForEvaluation (revalElement, &
      iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
      rvecMon%p_rspatialDiscr%p_rtriangulation, iel, &
      ctrafoType, DparPoint, Dpoint(:))
  
  ! Call the element to calculate the values of the basis functions
  ! in the point
  call elem_generic2(ieltype, revalElement, Bder, Dbas)
  
  ! Combine the basis functions to get the function value
  Dvalues(:) = 0.0_DP
  
  ! Now that we have the basis functions, we want to have the function values.
  ! We get them by multiplying the FE-coefficients with the values of the
  ! basis functions and summing up.
  !
  ! Calculate the value in the point
  do ibas = 1, indof
    Dvalues(1) = Dvalues(1) + p_DdataGradX(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(2) = Dvalues(2) + p_DdataGradY(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(3) = Dvalues(3) + p_DdataMon(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(4) = Dvalues(4) + p_DdataArea(Idofs(ibas)) * Dbas(ibas,iderType)
  end do
  
  Dvalues(3) = 1.0_DP/Dvalues(3)
  Dvalues(4) = 1.0_DP/Dvalues(4)
  
  end subroutine ! end griddef_evalPhi_Known

 !****************************************************************************

!<subroutine>
  subroutine griddef_evalPhi_Known3D(iderType, Dvalues, rvecGradX, &
             rvecGradY,rvecGradZ,rvecMon,rvecArea,Dpoint, &
             Ielements, IelementsHint, cnonmeshPoints)
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvecMon in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! type of function value to evaluate. One of the DER_xxxx constants
  integer, intent(in)           :: idertype
  
  ! the scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in)              :: rvecMon
  
  type(t_vectorScalar), intent(in)              :: rvecGradX
  
  type(t_vectorScalar), intent(in)              :: rvecGradY
  
  type(t_vectorScalar), intent(in)              :: rvecGradZ
  
  type(t_vectorScalar), intent(in)              :: rvecArea
  
  ! a list of points where to evaluate
  ! dimension(1..ndim,1..npoints)
  real(dp), dimension(:), intent(in) :: Dpoint
  
  ! optional: a list of elements containing the points Dpoints.
  ! if this is not specified the elment numbers containing the points
  ! are determined automatically
  integer, intent(in), optional :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements
  ! containing the points. This is ignored if Ielements is specified!
  integer, intent(in), optional :: IelementsHint
  
  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
  ! if a point is located outside of the domain. May happen e.g. in
  ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default
  ! parameter if cnonmeshPoints is not specified.
  integer, intent(in), optional :: cnonmeshPoints
  
!</input>


!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(dp), dimension(:), intent(out) :: Dvalues
  

  
!</output>

  ! local variables
  integer :: cnonmesh
  integer :: ieltype,indof,nve,ibas
  integer :: iel
  integer, dimension(:), pointer :: p_IelementDistr
  logical, dimension(el_maxnder) :: Bder
  
  real(dp), dimension(:), pointer :: p_ddatamon
  real(dp), dimension(:), pointer :: p_ddataarea
  real(dp), dimension(:), pointer :: p_ddatagradx
  real(dp), dimension(:), pointer :: p_ddatagrady
  real(dp), dimension(:), pointer :: p_ddatagradz
  
  ! Transformation
  integer :: ctrafotype
  real(dp), dimension(trafo_maxdimrefcoord) :: dparpoint
  
  ! Values of basis functions and doF`s
  real(DP), dimension(el_maxnbas,el_maxnder) :: Dbas
  integer, dimension(el_maxnbas) :: Idofs
  
  ! List of element distributions in the discretisation structure
  type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

  ! Evaluation structure and tag
  type(t_evalElement) :: revalElement
  integer :: cevaluationTag
  
  p_RelementDistribution => rvecMon%p_rspatialDiscr%RelementDistr
  
  ! for uniform discretisations, we get the element type in advance...
  if(rvecMon%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNifORM) then
  
    ! Element type
    ieltype = rvecMon%p_rspatialDiscr%RelementDistr(1)%celement
    
    ! get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! number of vertices on the element
    nve = elem_igetNDofLoc(ieltype)
    
    ! type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! Get the element evaluation tag; neccessary for preparation
    cevaluationTag = elem_getEvaluationTag(ieltype)
    
    nullify(p_IelementDistr)
  else
     call storage_getbase_int (&
          rvecMon%p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
  end if
  
  ! get the data vector
  select case (rvecMon%cdataType)
  case (ST_doUBLE)
    call lsyssc_getbase_double(rvecMon,p_DdataMon)
    call lsyssc_getbase_double(rvecArea,p_DdataArea)
    call lsyssc_getbase_double(rvecGradX,p_DdataGradX)
    call lsyssc_getbase_double(rvecGradY,p_DdataGradY)
    call lsyssc_getbase_double(rvecGradZ,p_DdataGradZ)
  case (ST_SINGLE)
  case default
    call output_line ('Unsupported vector precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
    call sys_halt()
  end select

  ! What to evaluate
  Bder = .false.
  Bder(iderType) = .true.
  
  cnonmesh = FEVL_NONMESHPTS_NONE
  if (present(cnonmeshPoints)) cnonmesh = cnonmeshPoints
  
  ! We loop over all points
  iel = 1
  
  ! Get the element number that contains the point
  if(present(Ielements))then
    ! we have it...
    iel = Ielements
  end if
  
  ! get the type of the element iel
  if(associated(p_ielementdistr))then
    ieltype = p_RelementDistribution(p_IelementDistr(iel))%celement
    
    ! Get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(ieltype)
    
    ! Type of transformation from /to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(ieltype)
  end if
  
  ! Calculate the global doF`s on that element into IdofsTest.
  call dof_locGlobMapping(rvecMon%p_rspatialDiscr, iel, Idofs)
  
  ! Get the element shape information
  call elprep_prepareForEvaluation(revalElement, EL_EVLTAG_COORDS, &
       rvecMon%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)
       
  ! calculate the transformation of the point to the reference element
  call trafo_calcRefCoords(ctrafoType, revalElement%Dcoords, &
       Dpoint(:), DparPoint)
       
  ! Now calculate everything else what is necessary for the element
  call elprep_prepareForEvaluation (revalElement, &
      iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
      rvecMon%p_rspatialDiscr%p_rtriangulation, iel, &
      ctrafoType, DparPoint, Dpoint(:))
  
  ! Call the element to calculate the values of the basis functions
  ! in the point
  call elem_generic2(ieltype, revalElement, Bder, Dbas)
  
  ! Combine the basis functions to get the function value
  Dvalues(:) = 0.0_DP
  
  ! Now that we have the basis functions, we want to have the function values.
  ! We get them by multiplying the FE-coefficients with the values of the
  ! basis functions and summing up.
  !
  ! Calculate the value in the point
  do ibas = 1, indof
    Dvalues(1) = Dvalues(1) + p_DdataGradX(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(2) = Dvalues(2) + p_DdataGradY(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(3) = Dvalues(3) + p_DdataGradZ(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(4) = Dvalues(4) + p_DdataMon(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(5) = Dvalues(5) + p_DdataArea(Idofs(ibas)) * Dbas(ibas,iderType)
  end do
  
  Dvalues(4) = 1.0_DP/Dvalues(4)
  Dvalues(5) = 1.0_DP/Dvalues(5)
  
  end subroutine ! end griddef_evalPhi_Known3D

  
  !****************************************************************************
  
!<subroutine>
  subroutine griddef_evalphi_ray(iderType, Dvalues, rvecGradX, &
             rvecGradY,rvecMon,rvecArea,Dpoint,bsearchFailed, &
             ielement)
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvecMon in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! type of function value to evaluate. One of the DER_xxxx constants
  integer, intent(in)           :: iderType
  
  ! the scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in)              :: rvecMon
  
  type(t_vectorScalar), intent(in)              :: rvecGradX
  
  type(t_vectorScalar), intent(in)              :: rvecGradY
  
  type(t_vectorScalar), intent(in)              :: rvecArea
  
  ! a list of points where to evaluate
  ! dimension(1..ndim,1..npoints)
  real(dp), dimension(:), intent(in) :: Dpoint
  
  ! Previous element containing the point
  integer, intent(in) :: ielement
  
!</input>


!<inputoutput>
  logical, intent(inout) :: bsearchFailed
!</inputoutput>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(dp), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  integer :: ieltype,indof,nve,ibas
  integer :: iel
  integer, dimension(:), pointer :: p_IelementDistr
  logical, dimension(el_maxnder) :: Bder
  
  real(dp), dimension(:), pointer :: p_DdataMon
  real(dp), dimension(:), pointer :: p_DdataArea
  real(dp), dimension(:), pointer :: p_DdataGradX
  real(dp), dimension(:), pointer :: p_DdataGradY
  
  ! Transformation
  integer :: ctrafotype,iresult,ilastElement,ilastEdge
  real(dp), dimension(trafo_maxdimrefcoord) :: DparPoint
  
  ! Values of basis functions and doF`s
  real(dp), dimension(el_maxnbas,el_maxnder) :: Dbas
  integer, dimension(el_maxnbas) :: Idofs
  
  ! List of element distributions in the discretisation structure
  type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

  ! Evaluation structure and tag
  type(t_evalElement) :: revalElement
  integer :: cevaluationTag
  
  p_RelementDistribution => rvecMon%p_rspatialDiscr%RelementDistr
  
  ! for uniform discretisations, we get the element type in advance...
  if(rvecMon%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNifORM) then
  
    ! Element type
    ieltype = rvecMon%p_rspatialDiscr%RelementDistr(1)%celement
    
    ! get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! number of vertices on the element
    nve = elem_igetNDofLoc(ieltype)
    
    ! type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! Get the element evaluation tag; neccessary for preparation
    cevaluationTag = elem_getEvaluationTag(ieltype)
    
    nullify(p_IelementDistr)
  else
     call storage_getbase_int (&
          rvecMon%p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
  end if
  
  ! get the data vector
  select case (rvecMon%cdataType)
  case (ST_doUBLE)
    call lsyssc_getbase_double(rvecMon,p_DdataMon)
    call lsyssc_getbase_double(rvecArea,p_DdataArea)
    call lsyssc_getbase_double(rvecGradX,p_DdataGradX)
    call lsyssc_getbase_double(rvecGradY,p_DdataGradY)
  case (ST_SINGLE)
  case default
    call output_line ('Unsupported vector precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
    call sys_halt()
  end select

  ! What to evaluate
  Bder = .false.
  Bder(iderType) = .true.
  
  ! We loop over all points
  iel = ielement
  
  ! Use raytracing search to find the element
  ! containing the point.
  call tsrch_getElem_raytrace2D (&
    Dpoint(:),rvecMon%p_rspatialDiscr%p_rtriangulation,iel,&
    iresult,ilastElement,ilastEdge,200)
  ! Fehler, wenn die Iterationen ausgehen wird
  ! das letzte Element zurückgegeben... verkehrt...
  ! Ok, not found... Brute force search
  
  if((iel .eq. 0) .or. (iresult .le. 0))then
    call tsrch_getElem_BruteForce (Dpoint(:), &
      rvecMon%p_rspatialDiscr%p_rtriangulation,iel)
  end if
  
  if((iel .eq. 0))then
    bsearchFailed = .true.
    return
  end if
  
  ! get the type of the element iel
  if(associated(p_IelementDistr))then
    ieltype = p_RelementDistribution(p_IelementDistr(iel))%celement
    
    ! Get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(ieltype)
    
    ! Type of transformation from /to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(ieltype)
  end if
  
  ! Calculate the global doF`s on that element into IdofsTest.
  call dof_locGlobMapping(rvecMon%p_rspatialDiscr, iel, Idofs)
  
  ! Get the element shape information
  call elprep_prepareForEvaluation(revalElement, EL_EVLTAG_COORDS, &
       rvecMon%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)
       
  ! calculate the transformation of the point to the reference element
  call trafo_calcRefCoords(ctrafoType, revalElement%Dcoords, &
       Dpoint(:), DparPoint)
       
  ! Now calculate everything else what is necessary for the element
  call elprep_prepareForEvaluation (revalElement, &
      iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
      rvecMon%p_rspatialDiscr%p_rtriangulation, iel, &
      ctrafoType, DparPoint, Dpoint(:))
  
  ! Call the element to calculate the values of the basis functions
  ! in the point
  call elem_generic2(ieltype, revalElement, Bder, Dbas)
  
  ! Combine the basis functions to get the function value
  Dvalues(:) = 0.0_DP
  
  ! Now that we have the basis functions, we want to have the function values.
  ! We get them by multiplying the FE-coefficients with the values of the
  ! basis functions and summing up.
  !
  ! Calculate the value in the point
  do ibas = 1, indof
    Dvalues(1) = Dvalues(1) + p_DdataGradX(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(2) = Dvalues(2) + p_DdataGradY(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(3) = Dvalues(3) + p_DdataMon(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(4) = Dvalues(4) + p_DdataArea(Idofs(ibas)) * Dbas(ibas,iderType)
  end do
  
  Dvalues(3) = 1.0_DP/Dvalues(3)
  Dvalues(4) = 1.0_DP/Dvalues(4)
  
  
  end subroutine ! end griddef_evalphi_ray

  !****************************************************************************
  
!<subroutine>
  subroutine griddef_evalphi_ray3D(iderType, Dvalues, rvecGradX, &
             rvecGradY,rvecGradZ,rvecMon,rvecArea,Dpoint,bsearchFailed, &
             ielement)
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvecMon in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! type of function value to evaluate. One of the DER_xxxx constants
  integer, intent(in)           :: iderType
  
  ! the scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in)              :: rvecMon
  
  type(t_vectorScalar), intent(in)              :: rvecGradX
  
  type(t_vectorScalar), intent(in)              :: rvecGradY
  
  type(t_vectorScalar), intent(in)              :: rvecGradZ
  
  type(t_vectorScalar), intent(in)              :: rvecArea
  
  ! a list of points where to evaluate
  ! dimension(1..ndim,1..npoints)
  real(dp), dimension(:), intent(inout) :: Dpoint
  
  ! Previous element containing the point
  integer, intent(in) :: ielement
  
!</input>


!<inputoutput>
  logical, intent(inout) :: bsearchFailed
!</inputoutput>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(dp), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  integer :: ieltype,indof,nve,ibas
  integer :: iel
  integer, dimension(:), pointer :: p_IelementDistr
  logical, dimension(el_maxnder) :: Bder
  
  real(dp), dimension(:), pointer :: p_DdataMon
  real(dp), dimension(:), pointer :: p_DdataArea
  real(dp), dimension(:), pointer :: p_DdataGradX
  real(dp), dimension(:), pointer :: p_DdataGradY
  real(dp), dimension(:), pointer :: p_DdataGradZ
  ! Transformation
  integer :: ctrafotype,iresult,ilastElement,ilastEdge
  real(dp), dimension(trafo_maxdimrefcoord) :: DparPoint
  
  ! Values of basis functions and doF`s
  real(dp), dimension(el_maxnbas,el_maxnder) :: Dbas
  integer, dimension(el_maxnbas) :: Idofs
  
  ! List of element distributions in the discretisation structure
  type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

  ! Evaluation structure and tag
  type(t_evalElement) :: revalElement
  integer :: cevaluationTag
  
  p_RelementDistribution => rvecMon%p_rspatialDiscr%RelementDistr
  
  ! for uniform discretisations, we get the element type in advance...
  if(rvecMon%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNifORM) then
  
    ! Element type
    ieltype = rvecMon%p_rspatialDiscr%RelementDistr(1)%celement
    
    ! get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! number of vertices on the element
    nve = elem_igetNDofLoc(ieltype)
    
    ! type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! Get the element evaluation tag; neccessary for preparation
    cevaluationTag = elem_getEvaluationTag(ieltype)
    
    nullify(p_IelementDistr)
  else
     call storage_getbase_int (&
          rvecMon%p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
  end if
  
  ! get the data vector
  select case (rvecMon%cdataType)
  case (ST_doUBLE)
    call lsyssc_getbase_double(rvecMon,p_DdataMon)
    call lsyssc_getbase_double(rvecArea,p_DdataArea)
    call lsyssc_getbase_double(rvecGradX,p_DdataGradX)
    call lsyssc_getbase_double(rvecGradY,p_DdataGradY)
    call lsyssc_getbase_double(rvecGradZ,p_DdataGradZ)
  case (ST_SINGLE)
  case default
    call output_line ('Unsupported vector precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
    call sys_halt()
  end select

  ! What to evaluate
  Bder = .false.
  Bder(iderType) = .true.
  
  ! We loop over all points
  iel = ielement
  ! Use raytracing search to find the element
  ! containing the point.
  
  call tsrch_getElem_raytrace3D (&
    Dpoint(:),rvecMon%p_rspatialDiscr%p_rtriangulation,iel,&
    iresult,ilastElement,ilastEdge,200)
  ! Fehler, wenn die Iterationen ausgehen wird
  ! das letzte Element zurückgegeben... verkehrt...
  ! Ok, not found... Brute force search

  if((iel .eq. 0) .or. (iresult .le. 0))then
    print *,"Brute force"
    call tsrch_getElem_BruteForce (Dpoint(:), &
      rvecMon%p_rspatialDiscr%p_rtriangulation,iel)
  end if
  
  if((iel .eq. 0))then
    print *,"daneben"
    bsearchFailed = .true.
    return
  end if
  
  ! get the type of the element iel
  if(associated(p_IelementDistr))then
    ieltype = p_RelementDistribution(p_IelementDistr(iel))%celement
    
    ! Get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(ieltype)
    
    ! Type of transformation from /to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(ieltype)
  end if
  
  ! Calculate the global doF`s on that element into IdofsTest.
  call dof_locGlobMapping(rvecMon%p_rspatialDiscr, iel, Idofs)
  
  ! Get the element shape information
  call elprep_prepareForEvaluation(revalElement, EL_EVLTAG_COORDS, &
       rvecMon%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)
       
  ! calculate the transformation of the point to the reference element
  call trafo_calcRefCoords(ctrafoType, revalElement%Dcoords, &
       Dpoint(:), DparPoint)
       
  ! Now calculate everything else what is necessary for the element
  call elprep_prepareForEvaluation (revalElement, &
      iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
      rvecMon%p_rspatialDiscr%p_rtriangulation, iel, &
      ctrafoType, DparPoint, Dpoint(:))
  
  ! Call the element to calculate the values of the basis functions
  ! in the point
  call elem_generic2(ieltype, revalElement, Bder, Dbas)
  
  ! Combine the basis functions to get the function value
  Dvalues(:) = 0.0_DP
  
  ! Now that we have the basis functions, we want to have the function values.
  ! We get them by multiplying the FE-coefficients with the values of the
  ! basis functions and summing up.
  !
  ! Calculate the value in the point
  do ibas = 1, indof
    Dvalues(1) = Dvalues(1) + p_DdataGradX(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(2) = Dvalues(2) + p_DdataGradY(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(3) = Dvalues(3) + p_DdataGradZ(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(4) = Dvalues(4) + p_DdataMon(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(5) = Dvalues(5) + p_DdataArea(Idofs(ibas)) * Dbas(ibas,iderType)
  end do
  
  Dvalues(4) = 1.0_DP/Dvalues(4)
  Dvalues(5) = 1.0_DP/Dvalues(5)
  
  
  end subroutine ! end griddef_evalphi_ray

! ***************************************************************************
 
!<subroutine>
  subroutine griddef_evalphi_ray_bound(iderType, Dvalues, rvecGradX, &
             rvecGradY,rvecMon,rvecArea,Dpoint,iinelement)
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvecMon in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! type of function value to evaluate. One of the DER_xxxx constants
  integer, intent(in)           :: iderType
  integer, intent(in)           :: iinelement
  
  ! the scalar solution vector that is to be evaluated.
  type(t_vectorScalar), intent(in)              :: rvecMon
  
  type(t_vectorScalar), intent(in)              :: rvecGradX
  
  type(t_vectorScalar), intent(in)              :: rvecGradY
  
  type(t_vectorScalar), intent(in)              :: rvecArea
  
  ! a list of points where to evaluate
  ! dimension(1..ndim,1..npoints)
  real(dp), dimension(:), intent(in) :: Dpoint
  
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  real(dp), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  integer :: ieltype,indof,nve,ibas
  integer :: iel
  integer, dimension(:), pointer :: p_IelementDistr
  logical, dimension(el_maxnder) :: Bder
  
  real(dp), dimension(:), pointer :: p_DdataMon
  real(dp), dimension(:), pointer :: p_DdataArea
  real(dp), dimension(:), pointer :: p_DdataGradX
  real(dp), dimension(:), pointer :: p_DdataGradY
  
  ! Transformation
  integer :: ctrafoType
  real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DparPoint
  
  ! Values of basis functions and doF`s
  real(DP), dimension(EL_MAXNBAS,EL_MAXNDER) :: Dbas
  integer, dimension(EL_MAXNBAS) :: Idofs
  
  ! List of element distributions in the discretisation structure
  type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

  ! Evaluation structure and tag
  type(t_evalElement) :: revalElement
  integer :: cevaluationTag
  
  p_RelementDistribution => rvecMon%p_rspatialDiscr%RelementDistr
  
  ! for uniform discretisations, we get the element type in advance...
  if(rvecMon%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNifORM) then
  
    ! Element type
    ieltype = rvecMon%p_rspatialDiscr%RelementDistr(1)%celement
    
    ! get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! number of vertices on the element
    nve = elem_igetNDofLoc(ieltype)
    
    ! type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! Get the element evaluation tag; neccessary for preparation
    cevaluationTag = elem_getEvaluationTag(ieltype)
    
    nullify(p_IelementDistr)
  else
     call storage_getbase_int (&
          rvecMon%p_rspatialDiscr%h_IelementDistr,p_IelementDistr)
  end if
  
  ! get the data vector
  select case (rvecMon%cdataType)
  case (ST_doUBLE)
    call lsyssc_getbase_double(rvecMon,p_DdataMon)
    call lsyssc_getbase_double(rvecArea,p_DdataArea)
    call lsyssc_getbase_double(rvecGradX,p_DdataGradX)
    call lsyssc_getbase_double(rvecGradY,p_DdataGradY)
  case (ST_SINGLE)
  case default
    call output_line ('Unsupported vector precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
    call sys_halt()
  end select

  ! What to evaluate
  Bder = .false.
  Bder(iderType) = .true.
  
  ! We loop over all points
  iel = iinelement
  
  ! get the type of the element iel
  if(associated(p_IelementDistr))then
    ieltype = p_RelementDistribution(p_IelementDistr(iel))%celement
    
    ! Get the number of local doF`s for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(ieltype)
    
    ! Type of transformation from /to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(ieltype)
  end if
  
  ! Calculate the global doF`s on that element into IdofsTest.
  call dof_locGlobMapping(rvecMon%p_rspatialDiscr, iel, Idofs)
  
  ! Get the element shape information
  call elprep_prepareForEvaluation(revalElement, EL_EVLTAG_COORDS, &
       rvecMon%p_rspatialDiscr%p_rtriangulation, iel, ctrafoType)
       
  ! calculate the transformation of the point to the reference element
  call trafo_calcRefCoords(ctrafoType, revalElement%Dcoords, &
       Dpoint(:), DparPoint)
       
  ! Now calculate everything else what is necessary for the element
  call elprep_prepareForEvaluation (revalElement, &
      iand(cevaluationTag,not(EL_EVLTAG_COORDS)), &
      rvecMon%p_rspatialDiscr%p_rtriangulation, iel, &
      ctrafoType, DparPoint, Dpoint(:))
  
  ! Call the element to calculate the values of the basis functions
  ! in the point
  call elem_generic2(ieltype, revalElement, Bder, Dbas)
  
  ! Combine the basis functions to get the function value
  Dvalues(:) = 0.0_DP
  
  ! Now that we have the basis functions, we want to have the function values.
  ! We get them by multiplying the FE-coefficients with the values of the
  ! basis functions and summing up.
  !
  ! Calculate the value in the point
  do ibas = 1, indof
    Dvalues(1) = Dvalues(1) + p_DdataGradX(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(2) = Dvalues(2) + p_DdataGradY(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(3) = Dvalues(3) + p_DdataMon(Idofs(ibas)) * Dbas(ibas,iderType)
    Dvalues(4) = Dvalues(4) + p_DdataArea(Idofs(ibas)) * Dbas(ibas,iderType)
  end do
  
  Dvalues(3) = 1.0_DP/Dvalues(3)
  Dvalues(4) = 1.0_DP/Dvalues(4)

  end subroutine ! griddef_evalphi_ray_bound
 
 !****************************************************************************
  
!<subroutine>
  subroutine griddef_getAreaDeformed(rgriddefInfo,&
             rvectorAreaBlockQ0,rvectorAreaBlockQ1,rvectorAreaQ0)
  
  !<description>
    ! In this function we build the nodewise area distribution and
    ! the elementwise distribution, for debugging purposes
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo

    type(t_vectorBlock),intent(inout)  :: rvectorAreaBlockQ0
    type(t_vectorBlock),intent(inout)  :: rvectorAreaBlockQ1
    type(t_vectorScalar),intent(inout) :: rvectorAreaQ0
  !</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(dp), dimension(:,:), pointer :: p_DvertexCoords
    real(dp), dimension(:), pointer :: p_Darea
    integer :: iel
    real(dp), dimension(ndim2d,tria_maxnve2d) :: Dpoints
    integer :: ive,NLMAX

    type(t_blockDiscretisation) :: rprjDiscretisation
    type(t_blockDiscretisation) :: rdiscretisation
    
    
    NLMAX=rgriddefInfo%NLMAX
        
    ! Is everything here we need?
    if (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line ('h_DvertexCoords not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      call sys_halt()
    end if

    if (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      call sys_halt()
    end if
    
    ! Do we have (enough) memory for that array?
    if (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DelementVolume .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementVolume2D', 'DAREA', &
          int(rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation   %NEL+1,I32), ST_doUBLE, &
          rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    end if
    
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation   )
                                   
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_E011,CUB_G2X2,rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation)
                                   
    
    ! Get the arrays
    call storage_getbase_double2D (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    
    ! Set up an empty block vector
    call lsysbl_createVecBlockByDiscr(rdiscretisation,rvectorAreaBlockQ1,.true.)
        
    ! Create a discretisation structure for Q0, based on our
    ! previous discretisation structure:
    call spdiscr_duplicateBlockDiscr(rdiscretisation,rprjDiscretisation)
    call spdiscr_deriveSimpleDiscrSc (&
                 rdiscretisation%RspatialDiscr(1), &
                 EL_Q0, CUB_G2X2, rprjDiscretisation%RspatialDiscr(1))
                 
    ! Initialise a Q0 vector from the newly created discretisation
    call lsyssc_createVecByDiscr(rprjDiscretisation%RspatialDiscr(1), &
    rvectorAreaQ0,.true.)
    
    ! get the pointer to the entries of this vector
    call lsyssc_getbase_double(rvectorAreaQ0,p_Darea)
    
    ! Loop over all elements calculate the area
    ! and save it in our vector
    do iel=1,rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%NEL
      
      if (p_IverticesAtElement(4,iel) .eq. 0) then
        ! triangular element
        do ive=1,TRIA_NVETRI2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        end do
        p_Darea(iel) = gaux_getArea_tria2D(Dpoints)
      else
        ! quad element
        do ive=1,TRIA_NVEQUAD2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        end do
        p_Darea(iel) = gaux_getArea_quad2D(Dpoints)
      end if

    end do ! end iel
    ! now transform the q0 vector into a q1 vector
    ! Setup a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecFromScalar(rvectorAreaQ0,rvectorAreaBlockQ0,rprjDiscretisation)
 
 
    ! Take the original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution(rvectorAreaBlockQ0,rvectorAreaBlockQ1)

    !call lsyssc_releaseVector(rvectorAreaQ0)
    call spdiscr_releaseBlockDiscr(rdiscretisation)

  
  end subroutine ! griddef_getAreaDeformed
  
  ! ***************************************************************************
 
!<subroutine>
 subroutine griddef_qMeasureM1(rgriddefInfo)
!<description>
  !
  !
  !
!</description>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout)  :: rgriddefInfo

!</inputoutput>
!</subroutine>
 
 end subroutine ! griddef_qMeasureM1
 
! ***************************************************************************
  
!<subroutine>
  subroutine griddef_deformationInit3D(rgriddefInfo,NLMIN,NLMAX,iStyle,iadaptSteps,iodesteps)
  
!<description>
  ! This subroutine initialises the rgriddefinfo structure which describes the
  ! grid deformation algorithm. The values for the certain parameters defining
  ! the deformation process are read in from the master.dat file. This routine
  ! has to be called before starting grid deformation.
!</description>

!<inputoutput>
  ! We will fill this structure with useful values
  type(t_griddefInfo), intent(inout) :: rgriddefInfo
  
  ! user defined number of adaptation steps
  integer,intent(in), optional :: iadaptSteps
  
  ! user defined number of adaptation steps
  integer,intent(in), optional :: iodesteps
  
!</inputoutput>

!<input>
  integer, intent(in) :: NLMAX
  
  integer, intent(in) :: NLMIN
  
  ! the deformation method to use multilevel or classical
  ! GRIDDEF_CLASSICAL or GRIDDEF_MULTILEVEL
  integer, intent(in) :: iStyle
!</input>
!</subroutine>
  
  ! Local variables
  integer :: iaux
  ! pointers to the arrays, we need to fill
  ! with useful default parameters
  integer, dimension(:), pointer :: p_Dql2
  integer, dimension(:), pointer :: p_Dqlinfty
  integer, dimension(:), pointer :: p_nadapSteps
  integer, dimension(:), pointer :: p_calcAdapSteps
  integer, dimension(:), pointer :: p_nadapStepsReally
  integer, dimension(:), pointer :: p_nmaxCorrSteps
  integer, dimension(:), pointer :: p_ncorrSteps
  integer, dimension(:), pointer :: p_ilevelODE
  integer, dimension(:), pointer :: p_ilevelODECorr
  integer, dimension(:), pointer :: p_npreSmoothSteps
  integer, dimension(:), pointer :: p_npostSmoothSteps
  integer, dimension(:), pointer :: p_nintSmoothSteps
  integer, dimension(:), pointer :: p_nmonSmoothsteps
  integer, dimension(:), pointer :: p_cpreSmoothMthd
  integer, dimension(:), pointer :: p_cpostSmoothMthd
  integer, dimension(:), pointer :: p_cintSmoothMthd
  
  rgriddefInfo%dblendPar = 1.0_dp
  
  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_Dql2',&
                   NLMAX, ST_INT, rgriddefInfo%h_Dql2,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_Dqlinfty',&
                   NLMAX, ST_INT, rgriddefInfo%h_Dqlinfty,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_nadapSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nadapSteps,&
                   ST_NEWBLOCK_ZERO)
                   
  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_calcAdapSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_calcAdapSteps,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_nadapStepsReally',&
                   NLMAX, ST_INT, rgriddefInfo%h_nadapStepsReally,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_nmaxCorrSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nmaxCorrSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_ncorrSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_ncorrSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_ilevelODE',&
                   NLMAX, ST_INT, rgriddefInfo%h_ilevelODE,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_ilevelODECorr',&
                   NLMAX, ST_INT, rgriddefInfo%h_ilevelODECorr,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_npreSmoothSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_npreSmoothSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_npostSmoothSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_npostSmoothSteps,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_nintSmoothSteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nintSmoothSteps,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_nmonSmoothsteps',&
                   NLMAX, ST_INT, rgriddefInfo%h_nmonSmoothsteps,&
                   ST_NEWBLOCK_ZERO)


  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_cpreSmoothMthd',&
                   NLMAX, ST_INT, rgriddefInfo%h_cpreSmoothMthd,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_cpostSmoothMthd',&
                   NLMAX, ST_INT, rgriddefInfo%h_cpostSmoothMthd,&
                   ST_NEWBLOCK_ZERO)

  ! now allocate memory for the arrays
  call storage_new('griddef_deformationInit3D', 'rgriddefInfo%h_cintSmoothMthd',&
                   NLMAX, ST_INT, rgriddefInfo%h_cintSmoothMthd,&
                   ST_NEWBLOCK_ZERO)

    
  
  ! Now get all the just allocated pointers
  call storage_getbase_int(rgriddefInfo%h_Dql2,p_Dql2)
  
  call storage_getbase_int(rgriddefInfo%h_Dqlinfty,p_Dqlinfty)

  call storage_getbase_int(rgriddefInfo%h_nadapSteps,p_nadapSteps)
  
  call storage_getbase_int(rgriddefInfo%h_calcAdapSteps,p_calcAdapSteps)
  
  call storage_getbase_int(rgriddefInfo%h_nadapStepsReally,p_nadapStepsReally)
  
  call storage_getbase_int(rgriddefInfo%h_nmaxCorrSteps,p_nmaxCorrSteps)
  
  call storage_getbase_int(rgriddefInfo%h_ncorrSteps,p_ncorrSteps)

  call storage_getbase_int(rgriddefInfo%h_ilevelODE,p_ilevelODE)

  call storage_getbase_int(rgriddefInfo%h_ilevelODECorr,p_ilevelODECorr)
  
  call storage_getbase_int(rgriddefInfo%h_npreSmoothSteps,p_npreSmoothSteps)

  call storage_getbase_int(rgriddefInfo%h_npostSmoothSteps,p_npostSmoothSteps)

  call storage_getbase_int(rgriddefInfo%h_nintSmoothSteps,p_nintSmoothSteps)
    
  call storage_getbase_int(rgriddefInfo%h_nmonSmoothsteps,p_nmonSmoothsteps)
      
  call storage_getbase_int(rgriddefInfo%h_cpreSmoothMthd,p_cpreSmoothMthd)

  call storage_getbase_int(rgriddefInfo%h_cpostSmoothMthd,p_cpostSmoothMthd)
  
  call storage_getbase_int(rgriddefInfo%h_cintSmoothMthd,p_cintSmoothMthd)


  ! Set the deformation style to classical
  rgriddefInfo%cdefStyle = iStyle

  ! temp variable
  iaux = 0
  
  ! set the multigrid level for the deformation PDE
  rgriddefInfo%iminDefLevel = NLMAX
  
  ! Here we initialize the structure with the standard values according to
  ! the desired grid deformation method.
  select case(rgriddefInfo%cdefStyle)
  
      case(GRIDDEF_CLASSICAL)

        ! initialise... here should be the absolute PDE Level
        ! where the deformation takes place.
        iaux = 1
      
        ! initialise the tolerance
        rgriddefInfo%dtolQ = 1.0E10_DP
        
        if(present(iodesteps))then
        ! number of ode steps
          rgriddefInfo%ntimeSteps = iodesteps
        else
          rgriddefInfo%ntimeSteps = 20
        end if
       
        ! set the multigrid level for the deformation PDE(really ?)
        p_ilevelODE(rgriddefInfo%iminDefLevel) = rgriddefInfo%iminDefLevel

        ! number of smoothing steps for creating the monitor function
        p_nintSmoothSteps(rgriddefInfo%iminDefLevel) = 0
        
        ! smoothing method used
        p_CintSmoothMthd(rgriddefInfo%iminDefLevel) = 0

        ! number of smoothing steps for monitor function
        p_NMonSmoothSteps(rgriddefInfo%iminDefLevel) = 0
        
        ! set standard value for number of correction steps
        p_nmaxCorrSteps(iaux) = 0
        
        ! set standard value for number of ODE correction steps
        p_ilevelODECorr(iaux) = 0
        
        ! initialize the blending parameter
        rgriddefInfo%dblendPar = 1.0_dp
        
        if(present(iadaptSteps))then
          ! initialize the number of adaptation steps
          rgriddefInfo%nadaptionSteps = iadaptSteps
          
          ! set Adaptive control
          ! the number of steps on all levels is user defined
          p_calcAdapSteps(:) = GRIDDEF_USER
        else
          ! initialize the number of adaptation steps
          rgriddefInfo%nadaptionSteps = 1
          
          ! the number of steps on all levels is
          ! fixed
          p_calcAdapSteps(:) = GRIDDEF_FIXED
        end if
        
     case(GRIDDEF_MULTILEVEL)
     
     ! do something here in the future
      call output_line ('MULTILEVEL method not yet implemented', &
          OU_CLASS_ERROR,OU_MODE_STD,'griddef_deformationInit')
      call sys_halt()
        
  end select
  
  ! get number of correction steps on level iaux
  p_nmaxCorrSteps(rgriddefInfo%iminDefLevel) = 0

  ! get the minimum deformation level
  p_ilevelODECorr(rgriddefInfo%iminDefLevel) = rgriddefInfo%iminDefLevel
  
  ! store the value of NLMIN
  rgriddefInfo%NLMIN = NLMIN

  ! store the value of NLMAX
  rgriddefInfo%NLMAX = NLMAX
  
  ! Allocate memory for all the levels.
  allocate(rgriddefInfo%p_rhLevels(NLMIN:NLMAX))

  end subroutine

!****************************************************************************************

!<subroutine>
  subroutine griddef_performDeformation3D(rgriddefInfo, &
                                        h_Dcontrib,&
                                        bstartNew, blevelHasChanged, bterminate, &
                                        bdegenerated, imgLevelCalc, iiteradapt, ibcIdx,&
                                        def_monitorfct,rperfconfig)
  !<description>
    ! This subroutine is the main routine for the grid deformation process, as all
    ! necessary steps are included here. For performing grid deformation, it is sufficient
    ! to define a monitor function and call this routine or to have an error distribution
    ! at hand and call this subroutine.
  !</description>

  !<input>
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.

    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
    
    ! if true, start from scratch: new vectors, new boundary conditions structures
    logical, intent(in) :: bstartNew

    ! if true, adjust the vectors after level change since a previous deformation call
    logical, intent(in) :: blevelHasChanged

    logical, intent(in) :: bterminate
    
    ! number of adaptive iteration
    integer, intent(in) :: iiterAdapt

    ! multigrid level on which the simulation was computed
    integer, intent(in) :: imgLevelCalc

    ! index of boundary condition related to the deformation PDE
    integer, intent(in):: ibcIdx

    ! flag for grid checking: if true, the deformation process would lead to
    ! a grid with tangled elements
    logical , intent(in):: bdegenerated

    ! A callback routine for the monitor function
    include 'intf_monitorfct.inc'
    optional :: def_monitorfct

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig

  !</input>

  !<inoutput>
    ! handle of vector with elementwise error contributions
    integer, intent(inout):: h_Dcontrib
  !</inoutput>

!</subroutine>

  ! local variables
  integer :: nmaxCorrSteps, NLMAX,NLMIN,idef,i
  ! A scalar matrix and vector. The vector accepts the RHS of the problem
  ! in scalar form.
  integer, dimension(:), pointer :: p_nadapStepsReally
  real(dp), dimension(:), pointer :: p_DblendPar
  integer, dimension(:), pointer :: p_ilevelODE
  
  ! We deform on level nlmax
  NLMAX = rgriddefInfo%iminDefLevel
  NLMIN = rgriddefInfo%NLMIN
  ! get some needed pointers
  call storage_getbase_int(rgriddefInfo%h_nadapStepsReally,p_nadapStepsReally)
  
  call storage_getbase_int(rgriddefInfo%h_ilevelODE,p_ilevelODE)

  ! necessary for boundary projection only, not for adding the Neumann boundary
  ! condition

  ! create and modify all vectors necessary for deformation steps
  call griddef_prepareDeformation(rgriddefInfo)

  ! Compute the number of deformation steps
  call griddef_computeAdapSteps(rgriddefInfo)
  ! Compute the blending parameter
  
  ! compute number of deformation steps on current level ilevel
  !call griddef_computeAdapSteps(rparBlock, rgriddefInfo, ilevel, pfmon)

  nmaxCorrSteps = 0

  call storage_getbase_double(rgriddefInfo%h_DblendPar,p_DblendPar)

  do i=NLMIN,NLMAX
      call spdiscr_initBlockDiscr (rgriddefInfo%p_rhLevels(i)%rdiscretisation,1,&
           rgriddefInfo%p_rhLevels(i)%p_rtriangulation)
      
      call spdiscr_initDiscr_simple (rgriddefInfo%p_rhLevels(i)%rdiscretisation%RspatialDiscr(1),&
                                     EL_Q1_3D,CUB_G3_3D,&
                                     rgriddefInfo%p_rhLevels(i)%p_rtriangulation)
                                     
  end do

  ! loop over adaption cycles
  do idef = 1,rgriddefInfo%nadaptionSteps

    call output_lbrk ()
    print *,"Adaptation Step: ",idef
    call output_line ('-------------------')
    call output_lbrk ()

!    ! perform one deformation step: This routine appplies deformation to the
    call griddef_performOneDefStep3D(rgriddefInfo,&
                                     p_DblendPar(idef), NLMAX, NLMAX,&
                                     def_monitorfct,rperfconfig)
     
    ! nullifiy where neccesary
    if(idef .ne. rgriddefInfo%nadaptionSteps)then
      call griddef_cleanLevels(rgriddefInfo%p_rhLevels,NLMIN,NLMAX)
    end if
  end do ! end do

  end subroutine
 
! ****************************************************************************************

!<subroutine>
  subroutine griddef_performOneDefStep3D(rgriddefInfo,&
                                         dblendpar, ilevelODE, ilevel,&
                                         def_monitorfct,rperfconfig)
  !<description>
    ! This subroutine performs one deformation step of the enhanced deformation method.
  !</description>

  !<inputoutput>
    ! structure containing all parameter settings for grid deformation
    type(t_griddefInfo), intent(inout) :: rgriddefInfo
  !</inputoutput>
  
  !<input>
    ! absolute level on which the points are moved
    integer, intent(in) :: ilevelODE

    ! absolute level on which the deformation PDE is solved
    integer, intent(in) :: ilevel
    
    ! blending parameter for monitorfunction sf + (1-s)g
    real(DP), intent(inout) :: dblendpar

    ! A callback routine for the monitor function
    include 'intf_monitorfct.inc'
    optional :: def_monitorfct
    
    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
  !</input>

!</subroutine>

    ! local variables
    real(dp), dimension(:,:), pointer :: Dresults

    ! weighting parameter for Laplacian smoothing of the vector field
    real(dp) :: dscale1,dscale2
    logical :: bBlending
    
    ! An array of problem levels for the multigrid solver
!    type(t_level), dimension(:), pointer :: Rlevels
    
    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rDubDiscretisation
    
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(2), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    ! Error indicator during initialisation of the solver
    integer :: ierror,NLMAX, i,NLMIN
    
    NLMIN=rgriddefInfo%NLMIN
    NLMAX=rgriddefInfo%NLMAX
    
    ! initialise Dresults
    Dresults  => null()

    ! no blending
    bBlending = .true.
    
    !---------------------------------------------------------------------------------
    !                                 FOR EVERY LEVEL
    !---------------------------------------------------------------------------------
    
    do i=NLMIN,NLMAX
    
      ! compute the original area distribution g(x)
      call griddef_getArea3D(rgriddefInfo,i)

      if(present(def_monitorfct))then
        ! calculate the monitor/target grid area distribution f(x)
        call griddef_calcMonitorFunction(rgriddefInfo,i,def_monitorfct)
      else
        ! or build a test monitor function
        call griddef_buildMonFuncTest(rgriddefInfo, i)
      end if

      ! blend monitor with current area distribution, if necessary
      if(bBlending) then
        ! Compute the scaling factors dscale1 and dscale2 for f and g and scale them
        call griddef_normaliseFctsNum(rgriddefInfo,dscale1,dscale2,i)
        call griddef_blendmonitor(rgriddefInfo,dblendpar,i)
      end if

      ! normalise the reciprocal of the functions
      call griddef_normaliseFctsInv(rgriddefInfo,i,rperfconfig)
      
      ! create the matrix for the poisson problem
      call griddef_createMatrixDef(rgriddefInfo,i)

      ! create rhs for deformation problem
      call griddef_createRhs(rgriddefInfo,i,rperfconfig)
!
    end do
    
    !---------------------------------------------------------------------------------
    !     end                          FOR EVERY LEVEL
    !---------------------------------------------------------------------------------
    
    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFreal
    RfilterChain(2)%ifilterType = FILTER_SMALLL1TO0
    RfilterChain(2)%ismallL1to0component = 1
    
    
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)
    
    ! Set up a coarse grid solver.
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
    call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
    
    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      
      ! Create an ILU(0) smoother
      call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
    end do
    
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (rgriddefInfo%p_rhLevels(i)%rmatDeform,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))
    
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    p_rsolverNode%nmaxIterations = 20
    
    p_rsolverNode%depsRel = 1e-13
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rgriddefInfo%p_rhLevels(NLMAX)%rSolBlock,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%rrhsBlock,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%rtempBlock)

    !-----------------------------------------------------------------------------
    !                            MESH MOVING PHASE
    !-----------------------------------------------------------------------------

    call spdiscr_initBlockDiscr (rDubDiscretisation,3,&
                                 rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation)

    call spdiscr_deriveSimpleDiscrSc (&
                 rgriddefInfo%p_rhLevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                 EL_Q1_3D, CUB_G3_3D, rDubDiscretisation%RspatialDiscr(1))
                 
    call spdiscr_deriveSimpleDiscrSc (&
                 rgriddefInfo%p_rhLevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                 EL_Q1_3D, CUB_G3_3D, rDubDiscretisation%RspatialDiscr(2))
                 
    call spdiscr_deriveSimpleDiscrSc (&
                 rgriddefInfo%p_rhLevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                 EL_Q1_3D, CUB_G3_3D, rDubDiscretisation%RspatialDiscr(3))
                 
    ! initialise the block vector that should hold the solution
    call lsysbl_createVecBlockByDiscr(rDubDiscretisation,&
                                      rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock,.true.)
    
    ! get the recovered gradient of the solution
    call ppgrd_calcGradient(rgriddefInfo%p_rhLevels(NLMAX)%rSolBlock%RvectorBlock(1),&
                            rgriddefInfo%p_rhLevels(NLMAX)%rvecGradBlock,PPGRD_INTERPOL)

    ! Solve the ODE and move the mesh
    call griddef_moveMesh3D(rgriddefInfo)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    
    call spdiscr_releaseBlockDiscr(rDubDiscretisation)
    
  
  end subroutine ! end griddef_performOneDefStep
  
! ****************************************************************************************

!<subroutine>
  subroutine griddef_project2Boundary(rgriddefInfo)
!<description>
  ! This subroutine handles the boundary nodes in 3d griddeformation
  ! we do simple projection of the neighbouring inner Vertex to the
  ! the nearest boundary face
!</description>

!<inputoutput>
  ! structure containing all parameter settings for grid deformation
  type(t_griddefInfo), intent(inout) :: rgriddefInfo
!</inputoutput>

!</subroutine>
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IvertAtEdge
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx
  integer, dimension(:), pointer :: p_IelementsAtVertex
  integer, dimension(:), pointer :: p_InodalProperty
  integer, dimension(:,:), pointer :: p_IfacesAtElement
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(DP), dimension(:,:), pointer :: p_DvertexCoordsNew
  integer, dimension(:,:), pointer :: p_IverticesAtFace
  
  
!-------------------------------------------------------
  real(DP), dimension(3) :: DPoint
  real(DP), dimension(3) :: DR1
  real(DP), dimension(3) :: DR2
  real(DP), dimension(3) :: DQ
  real(DP), dimension(3) :: DQP,DP2
!-------------------------------------------------------
  
  integer :: NLMAX,i,iVindex,ied,iel,iedge,v1,v2,iVert,iface
  integer :: nmt,nvt
  
  NLMAX = rgriddefInfo%NLMAX
  nmt=rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%NMT
  nvt=rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%NVT
  
  call storage_getbase_int2D(&
       rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IedgesAtElement,&
       p_IedgesAtElement)

  call storage_getbase_int2D(&
       rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtEdge,&
       p_IvertAtEdge)
       
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtBoundary,&
  p_IverticesAtBoundary)
  
  call storage_getbase_int(&
      rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertexIdx,&
      p_IelementsAtVertexIdx)
      
  call storage_getbase_int(&
      rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IelementsAtVertex,&
      p_IelementsAtVertex)
      
  call storage_getbase_int (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_InodalProperty,&
  p_InodalProperty)

  call storage_getbase_double2D (rgriddefInfo%p_rhLevels(NLMAX)%p_rtriangulation%h_DvertexCoords,&
      p_DvertexCoords)

  call storage_getbase_double2D (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_DvertexCoords,&
      p_DvertexCoordsNew)

      
      

  ! get the pointer
  call storage_getbase_int2D(&
      rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IverticesAtFace,&
      p_IverticesAtFace)
  
  call storage_getbase_int2D(&
      rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%h_IfacesAtElement,&
      p_IfacesAtElement)
  
  ! take care of all boundary nodes
  do i=1,rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation%NVBD
    iVindex=p_IverticesAtBoundary(i)
    iel=p_IelementsAtVertex(p_IelementsAtVertexIdx(iVindex))
    iVert=-1
    ! find the neighbouring vertex
    ! that is an inner vertex(if it exists..)
    do ied=1,12
      iedge=p_IedgesAtElement(ied,iel)
      v1=p_IvertAtEdge(1,iedge)
      v2=p_IvertAtEdge(2,iedge)
      ! check whether the edge contains the vertex
      if((v1 .eq. iVindex) .or.(v2 .eq. iVindex))then
        ! which of the edge vertices is the current vertex
        if(v1 .eq. iVindex)then
          ! if it is an inner node we can use
          ! it to project
          if(p_InodalProperty(v2) .eq. 0)then
            iVert=v2
            exit
          end if
        else
          ! if it is an inner node we can use
          ! it to project
          if(p_InodalProperty(v1) .eq. 0)then
            iVert=v1
            exit
          end if
        end if
        
      end if ! end if v1 or v2 ...
      
    end do ! end for all edges
    
    ! if there is no neigbour that
    ! is an inner vertex, we cannot do
    ! the projection, skip to next vertex
    if(iVert .eq. -1) cycle
    
    ! get the plane from the hexahedron and project...
    ! look for the face that is a boundary face
    ! AND contains the vertex iVert
    call findFace(p_IverticesAtFace,p_InodalProperty,p_IfacesAtElement,nvt,nmt,iVindex,iface)
!    do iae=1,6
!      iface=p_IfacesAtElement(iae,iel)
!      ! found inner face... skip
!      if(p_InodalProperty(nvt+nmt+iface) .eq. 0)cycle
!      do ive=1,4
!        if(p_IverticesAtFace(ive,iface) .eq. iVert)exit
!      end do
!    end do
    
    ! we got the face now set up the corresponding
    ! plane and project
    DPoint(:)=p_DvertexCoords(:,iVindex)
    DP2(:)=p_DvertexCoords(:,p_IverticesAtFace(1,iface))
    DR1(:)=p_DvertexCoords(:,p_IverticesAtFace(2,iface)) ! &
    !       -p_DvertexCoords(:,p_IverticesAtFace(1,iface))
    DR1(:)=DR1(:)-DP2(:)
    
    DR2(:)=p_DvertexCoords(:,p_IverticesAtFace(3,iface)) !&
           !-p_DvertexCoords(:,p_IverticesAtFace(1,iface))
    DR2(:)=DR2(:)-DP2(:)
    
    DQ(:)=p_DvertexCoordsNew(:,iVert)
    DQP=(/0.0_dp,0.0_dp,0.0_dp/)

    call gaux_projectPointPlane(DPoint,DR1,DR2,DQ,DQP)

    p_DvertexCoordsNew(:,iVindex)=DQP(:)
    
  end do ! end for all boundary nodes
  
  contains
  
  subroutine findFace(IverticesAtFace,InodalProperty,IfacesAtElement,nvt,nmt,iV,iface)
  
  integer, dimension(:,:),intent(in) :: IverticesAtFace
  integer, dimension(:),intent(in)   :: InodalProperty
  integer, dimension(:,:),intent(in) :: IfacesAtElement
  integer, intent(in) :: nvt
  integer, intent(in) :: nmt
  integer, intent(in) :: iV
  integer, intent(inout) :: iface
  integer :: iae,ive
  
    ! get the plane from the hexahedron and project...
    ! look for the face that is a boundary face
    ! AND contains the vertex iVindex
    do iae=1,6
      iface=IfacesAtElement(iae,iel)
      ! found inner face... skip
      if(InodalProperty(nvt+nmt+iface) .eq. 0)cycle
      do ive=1,4
        if(IverticesAtFace(ive,iface) .eq. iV)return
      end do
    end do
  
  end subroutine findFace

  end subroutine ! end griddef_project2Boundary
 
end module

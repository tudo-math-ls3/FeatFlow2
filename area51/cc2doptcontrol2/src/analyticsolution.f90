!##############################################################################
!# ****************************************************************************
!# <name> analyticsolution </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module encapsules an analytic solution. The solution may be given
!# as stationary FE solution with an arbitrary underlying mesh, as a
!# nonstationary solution or as analytical expression.
!#
!# The following routines can be found here:
!#
!# 1.) ansol_init
!#     -> Basic initialisation
!#
!# 2.) ansol_done
!#     -> Clean up, release memory
!#
!# 3.) ansol_configAnalytical
!#     -> Configures a solution to be analytically given
!#
!#     ansol_configExpressions
!#     -> Configures a solution to be given as expressions or
!#        configures the expressions in a blended flow.
!#
!# 4.) ansol_configStationaryFile
!#     -> Reads a solution from a file, configures the function to be stationary
!#
!# 5.) ansol_configNonstationaryFile
!#     -> Reads a solution from a sequence of files, configures the function
!#        to be nonstationary
!#
!# 6.) ansol_configNonstatPrecalc
!#     -> Configures a nonstationary function based on a precalculated solution
!#
!# 7.) ansol_prepareEval =
!#       ansol_prepareEvalCollection / ansol_prepareEvalDirect
!#     -> Prepares the evaluation of an solution using a collection or
!#     -> Prepares the direct evaluation of an solution
!#
!# 8.) ansol_evaluate =
!#       ansol_evaluateByCollection / ansol_evaluateDirect
!#     -> Evaluates a solution using a collection or
!#     -> Evaluates a solution using an analytic solution structure
!#     -> Does not include any boundary conditions
!#
!# 9.) ansol_doneEval
!#       ansol_doneEvalCollection / ansol_doneEvalDirect
!#     -> Cleans up the evaluation of a solution using a collection or
!#     -> Cleans up a direct evaluation of a FE solution
!#
!# 10.) ansol_prjToVector
!#      -> Project a solution to a vector
!#
!# 11.) ansol_prjToSpaceTimeVector
!#      -> Project a solution to a space-time vector
!#
!# 12.) ansol_getSpaceTimeDiscr
!#      -> Obtains the space/time discretisation of an underlying
!#         space-time solution if possible.
!#
!# Notes
!# -----
!# * The so called "blended" functions allow to create a function by a blending
!#   of two sub-functions. Such functions are defined by an expression for each
!#   component. The expression must return a value in the range [0..1] where
!#   0 represents the 1st solution and 1 the 2nd.
!# </purpose>
!##############################################################################

module analyticsolution

  use fsystem
  use genoutput
  use fparser
  use boundary
  use triangulation
  use basicgeometry
  use derivatives
  use linearsystemscalar
  use linearsystemblock
  use collection
  use feevaluation
  use analyticprojection
  
  use spatialdiscretisation
  use timediscretisation
  use timeevaluation

  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy

  use spacetimevectors
  
  implicit none
  
  private
  
  public :: t_anSolution
  public :: ansol_init
  public :: ansol_done

  public :: ansol_configExpressions
  public :: ansol_configStationaryFile
  public :: ansol_configNonstationaryFile
  public :: ansol_configNonstatPrecalc

  public :: ansol_prepareEval
  public :: ansol_evaluate
  public :: ansol_doneEval

  public :: ansol_prepareEvalCollection
  public :: ansol_evaluateByCollection
  public :: ansol_doneEvalCollection

  public :: ansol_prepareEvalDirect
  public :: ansol_evaluateDirect
  public :: ansol_doneEvalDirect
  
  public :: ansol_prjToVector
  public :: ansol_prjToSpaceTimeVector
  
  public :: ansol_getSpaceTimeDiscr
  
!<constants>

!<constantblock description="Variables that can be used in the expressions of the solution.">
    character(LEN=10), dimension(4), parameter :: ANSOL_VARIABLES = &
      (/'TIME ','X    ','Y    ','Z    '/)
!</constantblock>

!<constantblock description="Type of the function solution">

  ! Analytically given as expressions in p_Sexpressions.
  integer, parameter, public :: ANSOL_TP_EXPRESSIONS = -3
  
  ! Analytically given.
  ! Such a function cannot be evaluated. The caller must check if the function type
  ! is ANSOL_TP_ANALYTICAL and may evaluate the solution using t_anSolution%iid
  ! as a hint what to evaluate. The evaluation routines will return immediately
  ! with an appropriate error flag set in case the function is tried to be evaluated!
  integer, parameter, public :: ANSOL_TP_ANALYTICAL = -2

  ! Zero function.
  integer, parameter, public :: ANSOL_TP_ZERO = -1
  
  ! Undefined.
  integer, parameter, public :: ANSOL_TP_UNDEFINED = 0
  
  ! Mesh-based, but solution is still undefined.
  integer, parameter, public :: ANSOL_TP_MBUNDEFINED = 1
  
  ! Stationary solution, read from file, mesh-based. Arbitrary level.
  integer, parameter, public :: ANSOL_TP_MBSTATIONARYFILE = 2

  ! Stationary solution, specified by the application, mesh-based. Arbitrary level.
  integer, parameter, public :: ANSOL_TP_MBSTATIONARY = 3

  ! Nonstationary solution, read from file, mesh-based. Arbitrary level.
  integer, parameter, public :: ANSOL_TP_MBNONSTATIONARYFILE = 4

  ! Nonstationary solution, specified by the application, mesh-based. Arbitrary level.
  integer, parameter, public :: ANSOL_TP_MBNONSTATIONARY = 5
  
  ! Blended function. The expressions for each component return a parameter
  ! value in [0,1] that blends between two sub-solutions.
  integer, parameter, public :: ANSOL_TP_BLENDED = 6

!</constantblock>

!<constants>

!<types>

!<typeblock>

  ! This type encapsules an analytic solution.
  type t_anSolution
  
    ! Type of the solution.
    ! One of the ANSOL_TP_xxxx constants.
    integer :: ctype = ANSOL_TP_UNDEFINED
    
    ! This is a user defined ID for the solution.
    integer :: iid = 0
    
    ! Number of components in the solution.
    integer :: ncomponents = 0
    
    ! Analytical expression that defines the solution for the all coordinates
    ! in case ctype = ANSOL_TP_MBNONSTATIONARY/ANSOL_TP_MBNONSTATIONARYFILE.
    character(SYS_STRLEN), dimension(:), pointer :: p_Sexpressions => null()
    
    ! Parser object for the expressions.
    type(t_fparser) :: rparserExpression
    
    ! Name of the file or file sequence with the solution.
    character(SYS_STRLEN) :: sfilename = ''
    
    ! Element type used for the discretisation of the solution.
    ! This is a user-defined ID and is used in the case where the
    ! solution is represented as a vector.
    integer :: ielementType = 0
    
    ! Name of the TRI file with the basic mesh corresponding to the solution.
    ! ='': Use the same mesh as for the computation of the solution.
    character(SYS_STRLEN) :: smesh = ''

    ! Underlying FE space
    type(t_feSpaceLevel) :: rfeSpace
    
    ! Solution vector containing a stationary solution.
    ! Only used if ctype = ANSOL_TP_MBSTATIONARY/ANSOL_TP_MBSTATIONARYFILE/
    ! ANSOL_TP_MBNONSTATIONARY/ANSOL_TP_MBNONSTATIONARYFILE.
    type(t_vectorBlock) :: rstationary
    
    ! Time discretisation for a nonstationary solution.
    ! Only used if ctype = ANSOL_TP_MBNONSTATIONARY/ANSOL_TP_MBNONSTATIONARYFILE.
    type(t_timeDiscretisation) :: rtimeDiscr
    
    ! Solution vector containing a nonstationary solution.
    ! Only used if ctype = ANSOL_TP_MBNONSTATIONARY/ANSOL_TP_MBNONSTATIONARYFILE.
    type(t_spaceTimeVector) :: rnonstationary

    ! Auxiliary variable: Current point in time.
    real(DP) :: dtime
    
    ! 1st sub-solution. Used for blended functions.
    type(t_anSolution), pointer :: p_rsubsolution1 => null()
    
    ! 2nd sub-solution. Used for blended functions.
    type(t_anSolution), pointer :: p_rsubsolution2 => null()
    
  end type

!</typeblock>

!</types>

  interface ansol_init
    module procedure ansol_init_zero
    module procedure ansol_init_analytical
    module procedure ansol_init_discr
    module procedure ansol_init_fespace
    module procedure ansol_init_tria
    module procedure ansol_init_file
    module procedure ansol_init_blended
  end interface

  interface ansol_prepareEval
    module procedure ansol_prepareEvalCollection
    module procedure ansol_prepareEvalDirect
  end interface
  
  interface ansol_evaluate
    module procedure ansol_evaluateByCollection
    module procedure ansol_evaluateDirect
  end interface

  interface ansol_doneEval
    module procedure ansol_doneEvalCollection
    module procedure ansol_doneEvalDirect
  end interface


contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_zero (rsolution,ncomponents)

!<description>
  ! Initialises a solution as zero solution, used without a mesh.
  !
  ! The function is initialised as zero solution by default. This can be
  ! changed using ansol_configExpressions.
!</description>

!<input>
  ! Number of components
  integer, intent(in) :: ncomponents
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_ZERO
    rsolution%ncomponents = ncomponents
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_blended (rsolution,ncomponents,rsolution1,rsolution2)

!<description>
  ! Initialises a solution as blended solution between two other solutions.
  ! The solution depends on an expression for each component. The expression
  ! must return a value in the range [0,1] for each point in time, where
  ! 0 represents the values of rsolution1 and 1 the values of rsolution2.
!</description>

!<input>
  ! Number of components
  integer, intent(in) :: ncomponents

  ! 1st solution.
  type(t_anSolution), intent(in), target :: rsolution1

  ! 2nd solution
  type(t_anSolution), intent(in), target :: rsolution2
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_BLENDED
    rsolution%ncomponents = ncomponents
    rsolution%p_rsubsolution1 => rsolution1
    rsolution%p_rsubsolution2 => rsolution2
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_analytical (rsolution,ncomponents,iid)

!<description>
  ! Initialises a solution structure as analytical solution with id iid.
  !
  ! The function is initialised as zero solution by default. This can be
  ! changed using ansol_configAnalytical.
!</description>

!<input>
  ! Number of components
  integer, intent(in) :: ncomponents

  ! User defined ID that is attached to the function.
  integer, intent(in) :: iid
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_ANALYTICAL
    rsolution%ncomponents = ncomponents
    rsolution%iid = iid
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_discr (rsolution,ilevel,rdiscr,ilevelDiscr,&
      ielementType,fgetDiscr,rcollection)

!<description>
  ! Initialises a solution structure based on a space discretisation.
  ! rdiscr describes the underlying discretisation on refinement level ilevel
  ! (mesh + domain). ilevelDiscr defines the refinement level of the mesh
  ! in rsolution; if necessary, the mesh in rdiscr is automatically refined.
  ! Only applicable for solutions given as analytic expressions.
!</description>

!<input>
  ! Basic discretisation structure of the solution.
  type(t_blockDiscretisation), intent(in), target :: rdiscr
  
  ! Refinement level of the discretisation structure rdiscr.
  integer, intent(in) :: ilevelDiscr
  
  ! Refinement level used for the solution. If this is > ilevelDiscr,
  ! the underlying mesh is automatically refined using regular refinement.
  integer, intent(in) :: ilevel
  
  ! User defined ID that identifies the FE space.
  integer, intent(in) :: ielementType

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(ilevel,rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    integer, intent(in) :: ilevel
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_MBUNDEFINED
    rsolution%ielementType = ielementType
    
    ! Create the underlying FE space
    call fesph_createFEspace (rsolution%rfeSpace,ilevel,&
      rdiscr,ilevelDiscr,fgetDiscr,rcollection)
      
    ! Determine the number of components
    rsolution%ncomponents = rsolution%rfeSpace%p_rdiscretisation%ncomponents

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_fespace (rsolution,rfespace,ielementType)

!<description>
  ! Initialises a solution structure based on a FEM-space structure.
!</description>

!<input>
  ! The underlying FEM space structure.
  type(t_feSpaceLevel), intent(in) :: rfeSpace

  ! OPTIONAL: Element id.
  integer, intent(in), optional :: ielementType
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_MBUNDEFINED
    
    ! Take the FE space.
    rsolution%rfeSpace = rfeSpace
    if (present(ielementType)) then
      rsolution%ielementType = ielementType
    end if
      
    ! Determine the number of components
    rsolution%ncomponents = rsolution%rfeSpace%p_rdiscretisation%ncomponents

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_tria (rsolution,ilevel,rtriangulation,ilevelTri,&
      ielementType,fgetDiscr,rcollection,rboundary)

!<description>
  ! Initialises a solution structure to be used with a mesh.
  ! rboundary/rtriangulation describes the underlying discretisation on refinement
  ! level ilevel (mesh + domain). ilevelDiscr defines the refinement level of the
  ! mesh in rsolution; if necessary, the mesh in rdiscr is automatically refined.
  ! Only applicable for solutions given as analytic expressions.
!</description>

!<input>
  ! Underlying mesh, coarse level.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! Refinement level of the triangulation rtriangulation.
  integer, intent(in) :: ilevelTri
  
  ! Refinement level used for the solution. If this is > ilevelDiscr,
  ! the underlying mesh is automatically refined using regular refinement.
  integer, intent(in) :: ilevel
  
  ! User defined ID that identifies the FE space.
  integer, intent(in) :: ielementType

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(ilevel,rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    integer, intent(in) :: ilevel
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: Boundary definition of the domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_MBUNDEFINED
    rsolution%ielementType = ielementType
    
    ! Create the underlying FE space
    call fesph_createFEspace (rsolution%rfeSpace,ilevelTri,&
      rtriangulation,ilevel,fgetDiscr,rcollection,rboundary)

    ! Determine the number of components
    rsolution%ncomponents = rsolution%rfeSpace%p_rdiscretisation%ncomponents

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_init_file (rsolution,ilevel,ndim,smesh,&
      ielementType,fgetDiscr,rcollection,rboundary)

!<description>
  ! Initialises a solution structure to be used with a mesh.
  ! smesh defines the underlying mesh and ilevel the destination refinement
  ! level. ndim defines the dimension of the mesh.
  ! Only applicable for solutions given as analytic expressions.
!</description>

!<input>
  ! Refinement level used for the solution.
  ! The underlying mesh is automatically refined using regular refinement
  ! up to this level. Level 1 identrifies the coarse mesh without refinement.
  integer, intent(in) :: ilevel

  ! Dimension of the underlying mesh.
  integer, intent(in) :: ndim
  
  ! Name of the file containing the mesh.
  character(len=*), intent(in) :: smesh
  
  ! User defined ID that identifies the FE space.
  integer, intent(in) :: ielementType

  ! A callback routine that creates a discretisation based on
  ! a triangulation.
  interface
    subroutine fgetDiscr(ilevel,rtriangulation,rdiscr,rboundary,rcollection)
    use triangulation
    use boundary
    use spatialdiscretisation
    use collection
    integer, intent(in) :: ilevel
    type(t_triangulation), intent(in) :: rtriangulation
    type(t_blockDiscretisation), intent(out) :: rdiscr
    type(t_collection), intent(inout), optional :: rcollection
    type(t_boundary), intent(in), optional :: rboundary
    end subroutine
  end interface
  
  ! OPTIONAL: Collection structure which is passed to fgetDiscr.
  type(t_collection), intent(inout), optional :: rcollection
  
  ! OPTIONAL: Structure defining the domain.
  type(t_boundary), intent(in), optional :: rboundary
!</input>

!<output>
  ! The solution to initialise.
  type(t_anSolution), intent(out) :: rsolution
!</output>

!</subroutine>

    rsolution%ctype = ANSOL_TP_MBUNDEFINED
    rsolution%ielementType = ielementType
    
    ! Create the underlying FE space
    call fesph_createFEspace (rsolution%rfeSpace,ilevel,smesh,ndim,&
      fgetDiscr,rcollection,rboundary)

    ! Determine the number of components
    rsolution%ncomponents = rsolution%rfeSpace%p_rdiscretisation%ncomponents

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_done (rsolution)

!<description>
  ! Cleans up a solution structure. Releases all allocated memory.
!</description>

!<inputoutput>
  ! The solution to clean up.
  type(t_anSolution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    select case (rsolution%ctype)
    case (ANSOL_TP_EXPRESSIONS)
      ! Release the parser and the expressions, that is all.
      if (associated(rsolution%p_Sexpressions)) then
        call fparser_release(rsolution%rparserExpression)
        deallocate(rsolution%p_Sexpressions)
      end if
      
    case (ANSOL_TP_BLENDED)
      ! Release the parser and the expressions, that is all.
      if (associated(rsolution%p_Sexpressions)) then
        call fparser_release(rsolution%rparserExpression)
        deallocate(rsolution%p_Sexpressions)
      end if

    case (ANSOL_TP_MBUNDEFINED)
      ! Release the mesh, that is all.
      call fesph_releaseFEspace(rsolution%rfeSpace)
      rsolution%smesh = ""
      
    case (ANSOL_TP_MBSTATIONARYFILE,ANSOL_TP_MBSTATIONARY)
      ! Release the mesh
      call fesph_releaseFEspace(rsolution%rfeSpace)
      rsolution%smesh = ""
      
      ! Release the stationary solution
      call lsysbl_releaseVector(rsolution%rstationary)
    
    case (ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)
      ! Release the mesh
      call fesph_releaseFEspace(rsolution%rfeSpace)
      rsolution%smesh = ""

      ! Release the solution
      call sptivec_releaseVector(rsolution%rnonstationary)
      
    end select
    
    rsolution%ctype = ANSOL_TP_UNDEFINED
    rsolution%iid = 0
    rsolution%ncomponents = 0
    rsolution%ielementType = 0

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_configExpressions (rsolution,Sexpressions)

!<description>
  ! Configures an analytically solution by specifying expressions for the
  ! variables.
  ! The function must have been initialised by ansol_init_zero.
!</description>

!<input>
  ! OPTIONAL: List of expressions for all the components in the solution.
  ! If not present, the function is set to be a zero function.
  character(len=*), dimension(:), intent(inout), optional :: Sexpressions
!</input>

!<inputoutput>
  ! The solution to set up.
  type(t_anSolution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Basic checks
    if ((rsolution%ctype .ne. ANSOL_TP_ANALYTICAL) .and. &
        (rsolution%ctype .ne. ANSOL_TP_EXPRESSIONS) .and. &
        (rsolution%ctype .ne. ANSOL_TP_BLENDED) .and. &
        (rsolution%ctype .ne. ANSOL_TP_ZERO)) then
      call output_line('Incorrect solution initialisation!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_configAnalytical')
      call sys_halt()
    end if

    ! Release old expressions if they exist.
    if (associated(rsolution%p_Sexpressions)) then
      call fparser_release(rsolution%rparserExpression)
      deallocate(rsolution%p_Sexpressions)
    end if

    if (present(Sexpressions)) then
      ! Allocate new memory and initialise the expressions
      rsolution%ncomponents = ubound(Sexpressions,1)
      allocate(rsolution%p_Sexpressions(rsolution%ncomponents))
      rsolution%p_Sexpressions(:) = Sexpressions
      
      ! Pipe all expressions through a newly created parser
      call fparser_create(rsolution%rparserExpression,rsolution%ncomponents)
      
      do i=1,rsolution%ncomponents
        if (trim(Sexpressions(i)) .ne. "") then
          call fparser_parseFunction (rsolution%rparserExpression, i, &
              Sexpressions(i),ANSOL_VARIABLES)
        else
          ! Empty string; use zero function
          call fparser_parseFunction (rsolution%rparserExpression, i, &
              "0",ANSOL_VARIABLES)
        end if
      end do
      
      if (rsolution%ctype .ne. ANSOL_TP_BLENDED) then
        rsolution%ctype = ANSOL_TP_EXPRESSIONS
      end if
      
    else
    
      ! Zero function
      rsolution%ctype = ANSOL_TP_ZERO
    
    end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_configStationaryFile (rsolution,sfilename,bformatted)

!<description>
  ! Configures a stationary solution which is read from a file.
  ! The function must have been initialised by ansol_init_file.
!</description>

!<input>
  ! Name of the file containing the solution.
  character(len=*), intent(in) :: sfilename

  ! Whether to read formatted or unformatted data from disc.
  logical, intent(in) :: bformatted
!</input>

!<inputoutput>
  ! The solution to set up.
  type(t_anSolution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: svector
    integer :: i
    type(t_vectorblock) :: rtempVector

    ! Basic checks
    if ((rsolution%ctype .ne. ANSOL_TP_MBUNDEFINED) .and. &
        (rsolution%ctype .ne. ANSOL_TP_MBSTATIONARYFILE)) then
      call output_line('Incorrect solution initialisation!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_configStationaryFile')
      call sys_halt()
    end if
    
    if (sfilename .eq. "") then
      call output_line('NO filename specified!!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_configStationaryFile')
      call sys_halt()
    end if

    ! Save the filename
    rsolution%sfilename = sfilename
    
    ! Create the vector holding the solution if not yet created.
    if (rsolution%ctype .eq. ANSOL_TP_MBUNDEFINED) then
      call lsysbl_createVectorBlock(rsolution%rfeSpace%p_rdiscretisation,&
          rsolution%rstationary)
    endif
    
    ! Read a vector from the file.
    call vecio_readBlockVectorHR (&
        rtempVector, svector, .true., 0, sfilename, bformatted)
        
    ! Copy the subvectors and release the memory.
    do i=1,min(rsolution%rstationary%nblocks,rtempVector%nblocks)
      call lsyssc_duplicateVector (rtempVector%RvectorBlock(i),&
          rsolution%rstationary%RvectorBlock(i),&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
    end do
    
    call lsysbl_releaseVector (rtempVector)

    rsolution%ctype = ANSOL_TP_MBSTATIONARYFILE

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_configNonstationaryFile (rsolution, dstartTime,dendTime,ntimesteps,&
      sfilename,istart,idelta,bformatted)

!<description>
  ! Configures a nonstationary solution which is read from a file.
  ! The function must have been initialised by ansol_init_file.
  !
  ! sfilename is a directory/file name pattern in the format of
  ! a FORMAT statement that forms the filename; this pattern must contain
  ! exactly one integer format specifier, which is replaced by the
  ! file number in this routine (e.g. ' (''vector.txt.'',I5.5) ' will
  ! load a file sequence 'vector.txt.00001','vector.txt.00002',
  ! 'vector.txt.00003', ...).
  ! istart and iend prescribe the minimum/maximum file number that is
  ! inserted into the filename: The routine loads in all all files
  ! from istart to iend and forms a space-time vector from that.
!</description>

!<input>
  ! Start time of the nonstationary solution.
  real(dp), intent(in) :: dstartTime
  
  ! End time of the nonstationary solution.
  real(dp), intent(in) :: dendTime
  
  ! Filename pattern + path where to form a filename from.
  character(LEN=*), intent(IN) :: sfilename

  ! Number of the first file to be read in
  integer, intent(IN) :: istart
  
  ! Number of timesteps in the nonstationary solution.
  integer, intent(in) :: ntimesteps
  
  ! Delta parameter that specifies how to increase the filename suffix.
  ! Standard is =1.
  integer, intent(IN) :: idelta
  
  ! Whether to read formatted or unformatted data from disc.
  logical, intent(IN) :: bformatted

!</input>

!<inputoutput>
  ! The solution to set up.
  type(t_anSolution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    ! Basic checks
    if ((rsolution%ctype .ne. ANSOL_TP_MBUNDEFINED) .and. &
        (rsolution%ctype .ne. ANSOL_TP_MBSTATIONARYFILE)) then
      call output_line('Incorrect solution initialisation!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_configNonstationaryFile')
      call sys_halt()
    end if
    
    ! Save the filename
    rsolution%sfilename = sfilename
    
    ! Create the vector holding the solution if not yet created.
    if (rsolution%ctype .eq. ANSOL_TP_MBUNDEFINED) then
      call tdiscr_initOneStepTheta (rsolution%rtimeDiscr,&
          dstartTime, dendTime - mod(ntimesteps,idelta)*(dendTime-dstartTime)/idelta, &
          ntimesteps/idelta, 1.0_DP)

      call sptivec_initVector (rsolution%rnonstationary,&
          rsolution%rtimeDiscr,&
          rsolution%rfeSpace%p_rdiscretisation)

    endif
    
    ! Read the vector from the file.
    call sptivec_loadFromFileSequence (&
        rsolution%rnonstationary,&
        sfilename,istart,istart+ntimesteps,idelta,bformatted,.true.)

    rsolution%ctype = ANSOL_TP_MBNONSTATIONARYFILE

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_configNonstatPrecalc (rsolution, rvector)

!<description>
  ! Configures a nonstationary solution which is taken from a precalculated
  ! space-time solution vector.
  ! The function must have been initialised by ansol_init_file.
!</description>

!<input>
  ! Space-time vector that defines the precalculated solution vector.
  type(t_spaceTimeVector), intent(in) :: rvector
!</input>

!<inputoutput>
  ! The solution to set up.
  type(t_anSolution), intent(inout) :: rsolution
!</inputoutput>

!</subroutine>

    ! Basic checks
    if ((rsolution%ctype .ne. ANSOL_TP_MBUNDEFINED) .and. &
        (rsolution%ctype .ne. ANSOL_TP_MBSTATIONARYFILE)) then
      call output_line('Incorrect solution initialisation!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_configNonstationaryPrecalc')
      call sys_halt()
    end if
    
    ! Create the vector holding the solution if not yet created.
    if (rsolution%ctype .eq. ANSOL_TP_MBUNDEFINED) then
    
      call tdiscr_copy (rvector%p_rtimeDiscr,rsolution%rtimeDiscr)

      call sptivec_initVector (rsolution%rnonstationary,&
          rsolution%rtimeDiscr,&
          rsolution%rfeSpace%p_rdiscretisation)

    endif
    
    ! Copy the vector.
    call sptivec_copyVector (rvector,rsolution%rnonstationary)

    ! This is a solution specified by the application. No file is behind.
    rsolution%ctype = ANSOL_TP_MBNONSTATIONARY

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  recursive subroutine ansol_prepareEvalCollection (rsolution,rcollection,sname,dtime)

!<description>
  ! Adds a solution structure at a definite time as variable to a collection.
  ! The routine evaluates the solution at time dtime and adds this snapshot
  ! to the collection rcollection. (For stationary solutions, dtime is ignored.)
!</description>

!<input>
  ! Name of the solution.
  character(len=*), intent(in) :: sname
  
  ! The solution to be added to the collection.
  type(t_anSolution), intent(inout), target :: rsolution
  
  ! Point in time of the solution which is put to the collection.
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Collection structure where to add the solution to.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    type(t_vectorBlock), pointer :: p_rvector,p_rvector2

    ! Put the type and id to the collection.
    call collct_setvalue_int (rcollection, trim(sname)//"_CTYPE", &
        rsolution%ctype, .true.)

    call collct_setvalue_int (rcollection, trim(sname)//"_ID", &
        rsolution%iid, .true.)

    ! What to do depends on the type.
    select case (rsolution%ctype)
    case (ANSOL_TP_ANALYTICAL)
      call output_line('Solution undefined! Caller must evaluate the solution! Use iid as hint!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_prepareEvalCollection')
      call sys_halt()

    case (ANSOL_TP_ZERO)
      ! Nothing to do

    case (ANSOL_TP_BLENDED)

      ! Add a reference to the parser to the collection
      call collct_setvalue_pars (rcollection, trim(sname)//"_PARS", &
          rsolution%rparserExpression, .true.)
          
      ! Add the current time to the collection
      call collct_setvalue_real (rcollection, trim(sname)//"_TIME", &
          dtime, .true.)
          
      ! Prepare the evaluation of the sub-functions.
      call ansol_prepareEvalCollection (rsolution%p_rsubsolution1,rcollection,trim(sname)//"1",dtime)
      call ansol_prepareEvalCollection (rsolution%p_rsubsolution1,rcollection,trim(sname)//"2",dtime)

    case (ANSOL_TP_EXPRESSIONS)
    
      ! Add a reference to the parser to the collection
      call collct_setvalue_pars (rcollection, trim(sname)//"_PARS", &
          rsolution%rparserExpression, .true.)
          
      ! Add the current time to the collection
      call collct_setvalue_real (rcollection, trim(sname)//"_TIME", &
          dtime, .true.)
    
    case (ANSOL_TP_MBSTATIONARYFILE,ANSOL_TP_MBSTATIONARY)
    
      ! Add a reference to the stationary solution vector to the collection.
      call collct_setvalue_vec (rcollection, trim(sname)//"_VEC", &
          rsolution%rstationary, .true.)
    
    case (ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)

      ! Create a new vector for that timestep.
      allocate(p_rvector)
      call lsysbl_createVectorBlock(rsolution%rfeSpace%p_rdiscretisation,p_rvector,.false.)

      ! Evaluate the function at that time.
      call tmevl_evaluate(rsolution%rnonstationary,dtime,p_rvector)

      ! Add a reference to the stationary solution vector to the collection.
      call collct_setvalue_vec (rcollection, trim(sname)//"_VEC", p_rvector, .true.)

    case default
      call output_line('Case not implemented!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_doneEvalCollection')
      call sys_halt()
    end select


  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  recursive subroutine ansol_doneEvalCollection (rcollection,sname)

!<description>
  ! Removes a solution structure from a collection.
!</description>

!<input>
  ! Name of the solution.
  character(len=*), intent(in) :: sname
!</input>

!<inputoutput>
  ! Collection structure where to remove the solution from.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    type(t_vectorBlock), pointer :: p_rvector
    integer :: ctype
    
    ! Type of the solution?
    ctype = collct_getvalue_int (rcollection, trim(sname)//"_CTYPE")
    call collct_deletevalue (rcollection, trim(sname)//"_CTYPE")

    call collct_deletevalue (rcollection, trim(sname)//"_ID")

    ! What to do depends on the type.
    select case (ctype)
    case (ANSOL_TP_ANALYTICAL)
      call output_line('Solution undefined! Caller must evaluate the solution! Use iid as hint!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_doneEvalCollection')
      call sys_halt()

    case (ANSOL_TP_ZERO)
      ! Nothing to do

    case (ANSOL_TP_EXPRESSIONS)
    
      ! Remove our parser from the collection
      call collct_deletevalue (rcollection, trim(sname)//"_PARS")
      
      call collct_deletevalue (rcollection, trim(sname)//"_TIME")

    case (ANSOL_TP_BLENDED)
    
      ! Remove our parser from the collection
      call collct_deletevalue (rcollection, trim(sname)//"_PARS")
      
      call collct_deletevalue (rcollection, trim(sname)//"_TIME")
      
      ! Clean up the sub-solutions
      call ansol_doneEvalCollection(rcollection,trim(sname)//"1")
      call ansol_doneEvalCollection(rcollection,trim(sname)//"2")

    case (ANSOL_TP_MBSTATIONARYFILE,ANSOL_TP_MBSTATIONARY)
    
      ! Delete the reference to the stationary solution.
      call collct_deletevalue (rcollection, trim(sname)//"_VEC")
    
    case (ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)

      ! Release the vector from the heap.
      p_rvector => collct_getvalue_vec (rcollection, trim(sname)//"_VEC")
      call lsysbl_releaseVector(p_rvector)
      deallocate(p_rvector)

      ! Delete the reference to the stationary solution.
      call collct_deletevalue (rcollection, trim(sname)//"_VEC")

    case default
      call output_line('Case not implemented!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_doneEvalCollection')
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  recursive subroutine ansol_evaluateByCollection (rcollection,sname,idim,Dvalues,&
      npoints,nelements,Dpoints,Ielements,ierror,iid)

!<description>
  ! Evaluate component idim of a solution in a set of points on a set of
  ! elements. The solution is given as part of a collection structure
  ! rcollection. sname identifies the name of the solution.
  ! For nonstationary functions, the solution is given as a snapshot
  ! of a point in time. The solution must have been added to the collection
  ! using ansol_prepareEvalCollection.
  !
  ! Note: The routine will not implement any bondary conditions into the function!
!</description>

!<input>
  ! Collection structure holding the solution.
  type(t_collection), intent(inout) :: rcollection
  
  ! Name of the solution in the collection.
  character(len=*), intent(in) :: sname
  
  ! Number of the comonent to evaluate.
  integer, intent(in) :: idim
  
  ! Number of points per element
  integer, intent(in) :: npoints
  
  ! Number of elements where to evaluate
  integer, intent(in) :: nelements
  
  ! List of coordinates where to evaluate. This is a list of x/y
  ! coordinates for a couple of points on a couple of elements.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
  ! OPTIONAL: List of elements where to evaluate.
  integer, dimension(:), intent(in), optional :: Ielements
!</input>

!<output>
  ! Receives the computed values in all points on all elements.
  real(dp), dimension(:,:), intent(out) :: Dvalues
  
  ! OPTIONAL: Returns whether the evaluation was successfull.
  ! =0: success. =1: error. =-1: Function is analytical and cannot be evaluated.
  integer, intent(out), optional :: ierror
  
  ! OPTIONAL: Returns the user defined ID of the function.
  integer, intent(out), optional :: iid
!</output>
  
!</subroutine>

    real(DP) :: dtime
    integer :: i,j
    type(t_vectorBlock), pointer :: p_rvector
    
    integer :: iel,NEL,ctype
    real(DP), dimension(:), allocatable :: DvaluesAct
    real(DP), dimension(:,:), allocatable :: DpointsAct
    integer, dimension(:), allocatable :: IelementsAct
    real(dp), dimension(:,:), allocatable :: p_Dval
    
    type(t_fparser), pointer :: p_rparser
    real(DP) :: dweight
    real(DP), dimension(:,:), allocatable :: Dvalues2
    real(DP), dimension(4) :: DvaluesBlending
    
    ! DEBUG!!!
    !real(DP), dimension(:), pointer :: p_Ddata

    if (present(iid)) &
      iid = collct_getvalue_int (rcollection, trim(sname)//"_ID")
      
    if (present(ierror)) ierror = 0

    ! Type of the solution?
    ctype = collct_getvalue_int (rcollection, trim(sname)//"_CTYPE")

    select case (ctype)

    case (ANSOL_TP_BLENDED)
    
      ! Get the parser object with the RHS expressions from the collection
      p_rparser => collct_getvalue_pars (rcollection, trim(sname)//"_PARS")
      dtime = collct_getvalue_real (rcollection, trim(sname)//"_TIME")
      
      ! Call the parser, evaluare the weighting parameter.
      DvaluesBlending(:) = (/dtime,0.0_DP,0.0_DP,0.0_DP/)
      call fparser_evalFunction (p_rparser, idim, DvaluesBlending, dweight)
      dweight = min(max(dweight,0.0_DP),1.0_DP)

      ! Allocate a 2nd array that receives the values from the 2nd subfunction.
      allocate (Dvalues2(ubound(Dvalues,1),ubound(Dvalues,2)))
      
      ! Evaluate the functions
      call ansol_evaluateByCollection (rcollection,trim(sname)//"1",idim,Dvalues,&
          npoints,nelements,Dpoints,Ielements,ierror)
          
      if (ierror .ne. 0) then
        ! Error.
        Dvalues(:,:) = 0.0_DP
        deallocate(Dvalues2)
        return
      end if
          
      call ansol_evaluateByCollection (rcollection,trim(sname)//"2",idim,Dvalues2,&
          npoints,nelements,Dpoints,Ielements,ierror)

      if (ierror .ne. 0) then
        ! Error.
        Dvalues(:,:) = 0.0_DP
        deallocate(Dvalues2)
        return
      end if
          
      ! Calculate the blended value, release, finish.

      ! Reshape the data, that's it.
      do i=1,nelements
        do j=1,npoints
          Dvalues(j,i) = (1-dweight) * Dvalues(j,i) + dweight * Dvalues2(j,i)
        end do
      end do
      
      deallocate(Dvalues2)

    case (ANSOL_TP_EXPRESSIONS)
    
      ! Get the parser object with the RHS expressions from the collection
      p_rparser => collct_getvalue_pars (rcollection, trim(sname)//"_PARS")
      dtime = collct_getvalue_real (rcollection, trim(sname)//"_TIME")
      
      ! Prepare the array with the values for the function.
      allocate(p_Dval(4,npoints*nelements))
      select case (ubound(Dpoints,1))
      case (NDIM2D)
        do i=1,nelements
          do j=1,npoints
            p_Dval (1,(i-1)*npoints+j) = dtime
            p_Dval (2,(i-1)*npoints+j) = Dpoints(1,j,i)
            p_Dval (3,(i-1)*npoints+j) = Dpoints(2,j,i)
            p_Dval (4,(i-1)*npoints+j) = 0.0_DP
          end do
        end do

      case (NDIM3D)
        do i=1,nelements
          do j=1,npoints
            p_Dval (1,(i-1)*npoints+j) = dtime
            p_Dval (2,(i-1)*npoints+j) = Dpoints(1,j,i)
            p_Dval (3,(i-1)*npoints+j) = Dpoints(2,j,i)
            p_Dval (4,(i-1)*npoints+j) = Dpoints(3,j,i)
          end do
        end do
        
      case default
        call output_line('Case not implemented!',&
            OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateByCollection')
        call sys_halt()
      end select
      
      ! Evaluate the 1st expression for the X-rhs
      allocate(DvaluesAct(npoints*nelements))
      call fparser_evalFunction (p_rparser, idim, 2, p_Dval, DvaluesAct)

      ! Reshape the data, that's it.
      do i=0,nelements-1
        do j=1,npoints
          Dvalues(j,i+1) = DvaluesAct(i*npoints+j)
        end do
      end do
      
      deallocate(DvaluesAct)
      deallocate(p_Dval)
      
    case (ANSOL_TP_ANALYTICAL)
    
      ! The caller must evaluate this, we cannot do it here!
      if (present(ierror)) ierror = -1
      return
      
      !call output_line('Solution undefined! Caller must evaluate the solution! Use iid as hint!',&
      !    OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateByCollection')
      !call sys_halt()
      
    case (ANSOL_TP_ZERO)
      ! Zero target function
      Dvalues(:,:) = 0.0_DP
      
    case(ANSOL_TP_MBSTATIONARYFILE,ANSOL_TP_MBSTATIONARY,&
          ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)
      ! Stationary/nonstationary function, specified by a block vector.
      !
      ! For every point, find the element of an element nearby the point.
      ! The evaluation routine uses this as hint to speed up the evaluation.
      allocate(DpointsAct(NDIM2D,npoints*nelements))
      allocate(DvaluesAct(npoints*nelements))

      ! Get the function
      p_rvector => collct_getvalue_vec (rcollection, trim(sname)//"_VEC")
      
      if (present(Ielements)) then
      
        ! Evaluate with Ielements giving us a hint where in which elements
        ! the points can be found.
      
        allocate(IelementsAct(npoints*nelements))
        NEL = p_rvector%p_rblockDiscr%p_rtriangulation%NEL
        do i=0,nelements-1
          do j=1,npoints
            iel = Ielements(i+1)
            
            ! If the FE function is evaluated on a level higher than the
            ! discretisation allows, the element number may be out of bounds.
            ! Reduce it by dividing by 4 which simulates coarsening for a mesh
            ! which is refined by 2-level ordering.
            do while (iel .gt. NEL)
              iel = iel / 4
            end do
            
            IelementsAct(i*npoints+j) = iel
            
            DpointsAct(1,i*npoints+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npoints+j) = Dpoints(2,j,i+1)
          end do
        end do
        
        ! Evaluate at the given points.
        call fevl_evaluate (DER_FUNC, DvaluesAct, p_rvector%RvectorBlock(idim), &
          DpointsAct, IelementsHint=IelementsAct,cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
          
      else
      
        ! Evaluate with searching for the points.
        call fevl_evaluate (DER_FUNC, DvaluesAct, p_rvector%RvectorBlock(idim), &
          DpointsAct, cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
      
      end if
        
      deallocate(IelementsAct)
      
      do i=0,nelements-1
        do j=1,npoints
          Dvalues(j,i+1) = DvaluesAct(i*npoints+j)
        end do
      end do
      
      deallocate(DpointsAct)
      deallocate(DvaluesAct)
      
    case default
      call output_line('Function not valid!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateByCollection')
      call sys_halt()
    end select
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  recursive subroutine ansol_prepareEvalDirect (rsolution,dtime)

!<description>
  ! Must be called in advance of a call to ansol_evaluateDirect.
  ! Prepares the evaluation of the solution
!</description>

!<input>
  ! Time where the solution should be evaluated.
  ! For a stationary solution, this is ignored.
  real(dp), intent(in) :: dtime
!</input>

!<input>
  ! The solution to prepare for evaluation.
  type(t_anSolution), intent(inout), target :: rsolution
!</input>

!</subroutine>

    select case (rsolution%ctype)
    case (ANSOL_TP_BLENDED)
      ! Pass to the subfunctions
      call ansol_prepareEvalDirect (rsolution%p_rsubsolution1,dtime)
      call ansol_prepareEvalDirect (rsolution%p_rsubsolution2,dtime)
    
    case (ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)
      
      ! Create a solution vector as temporary space for evaluations.
      call lsysbl_createVectorBlock(rsolution%rfeSpace%p_rdiscretisation,&
          rsolution%rstationary)

      ! Evaluate the function at the specific time.
      ! Write the result to our temp vector.
      call tmevl_evaluate(rsolution%rnonstationary,dtime,rsolution%rstationary)
      
    end select

    ! Remember the current time.
    rsolution%dtime = dtime
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  recursive subroutine ansol_doneEvalDirect (rsolution)

!<description>
  ! Must be called when the evaluation of a solution is finished.
  ! Cleans up temporary data.
!</description>

!<input>
  ! The solution to prepare for evaluation.
  type(t_anSolution), intent(inout), target :: rsolution
!</input>

!</subroutine>

    select case (rsolution%ctype)
    case (ANSOL_TP_BLENDED)
      ! Pass to the subfunctions
      call ansol_doneEvalDirect (rsolution%p_rsubsolution1)
      call ansol_doneEvalDirect (rsolution%p_rsubsolution2)

    case (ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)
      
      ! Release the temp vector
      call lsysbl_releaseVector(rsolution%rstationary)
      
    end select
      
    rsolution%dtime = 0.0_DP
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  recursive subroutine ansol_evaluateDirect (rsolution,idim,Dvalues,&
      npoints,nelements,Dpoints,Ielements,ierror,iid)

!<description>
  ! Evaluate component idim of a solution in a set of points on a set of
  ! elements. The solution is given as part of a collection structure
  ! rcollection. sname identifies the name of the solution.
  ! For nonstationary functions, the solution is given as a snapshot
  ! of a point in time.
  !
  ! Before calling this routine, the evaluation must have been prepared
  ! with ansol_prepareEvalDirect!
  !
  ! Note: The routine will not implement any bondary conditions into the function!
!</description>

!<input>
  ! The solution to evaluate.
  type(t_anSolution), intent(inout), target :: rsolution
  
  ! Number of the component to evaluate.
  integer, intent(in) :: idim
  
  ! Number of points per element
  integer, intent(in) :: npoints
  
  ! Number of elements where to evaluate
  integer, intent(in) :: nelements
  
  ! List of coordinates where to evaluate. This is a list of x/y
  ! coordinates for a couple of points on a couple of elements.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints
  
  ! OPTIONAL: List of elements where to evaluate.
  integer, dimension(:), intent(in), optional :: Ielements
!</input>

!<output>
  ! Receives the computed values in all points on all elements.
  real(dp), dimension(:,:), intent(out) :: Dvalues

  ! OPTIONAL: Returns whether the evaluation was successfull.
  ! =0: success. =1: error. =-1: Function is analytical and cannot be evaluated.
  integer, intent(out), optional :: ierror
  
  ! OPTIONAL: Returns the user defined ID of the function.
  integer, intent(out), optional :: iid
!</output>
  
!</subroutine>

    integer :: i,j
    
    integer :: iel,NEL
    real(DP), dimension(:), allocatable :: DvaluesAct
    real(DP), dimension(:,:), allocatable :: DpointsAct
    integer, dimension(:), allocatable :: IelementsAct
    real(DP) :: dweight
    real(DP), dimension(:,:), allocatable :: Dvalues2
    real(DP), dimension(4) :: DvaluesBlending
    
    real(dp), dimension(:,:), allocatable :: p_Dval
    
    ! DEBUG!!!
    !real(DP), dimension(:), pointer :: p_Ddata

    if (present(iid)) &
      iid = rsolution%iid
      
    if (present(ierror)) ierror = 0

    ! Evaluate, depending on the type.
    
    select case (rsolution%ctype)
    case (ANSOL_TP_BLENDED)
    
      ! Call the parser, evaluare the weighting parameter.
      DvaluesBlending(:) = (/rsolution%dtime,0.0_DP,0.0_DP,0.0_DP/)
      call fparser_evalFunction (rsolution%rparserExpression, idim, DvaluesBlending, dweight)
      dweight = min(max(dweight,0.0_DP),1.0_DP)

      ! Allocate a 2nd array that receives the values from the 2nd subfunction.
      allocate (Dvalues2(ubound(Dvalues,1),ubound(Dvalues,2)))
      
      ! Evaluate the functions
      call ansol_evaluateDirect (rsolution%p_rsubsolution1,idim,Dvalues,&
          npoints,nelements,Dpoints,Ielements,ierror)
          
      if (ierror .ne. 0) then
        ! Error.
        Dvalues(:,:) = 0.0_DP
        deallocate(Dvalues2)
        return
      end if
          
      call ansol_evaluateDirect (rsolution%p_rsubsolution2,idim,Dvalues,&
          npoints,nelements,Dpoints,Ielements,ierror)

      if (ierror .ne. 0) then
        ! Error.
        Dvalues(:,:) = 0.0_DP
        deallocate(Dvalues2)
        return
      end if
          
      ! Calculate the blended value, release, finish.

      ! Reshape the data, that's it.
      do i=1,nelements
        do j=1,npoints
          Dvalues(j,i) = (1-dweight) * Dvalues(j,i) + dweight * Dvalues2(j,i)
        end do
      end do
      
      deallocate(Dvalues2)

    case (ANSOL_TP_EXPRESSIONS)
      ! Prepare the array with the values for the function.
      allocate(p_Dval(3,npoints*nelements))
      select case (ubound(Dpoints,1))
      case (NDIM2D)
        do i=1,nelements
          do j=1,npoints
            p_Dval (1,(i-1)*npoints+j) = rsolution%dtime
            p_Dval (2,(i-1)*npoints+j) = Dpoints(1,j,i)
            p_Dval (3,(i-1)*npoints+j) = Dpoints(2,j,i)
            p_Dval (4,(i-1)*npoints+j) = 0.0_DP
          end do
        end do

      case (NDIM3D)
        do i=1,nelements
          do j=1,npoints
            p_Dval (1,(i-1)*npoints+j) = rsolution%dtime
            p_Dval (2,(i-1)*npoints+j) = Dpoints(1,j,i)
            p_Dval (3,(i-1)*npoints+j) = Dpoints(2,j,i)
            p_Dval (4,(i-1)*npoints+j) = Dpoints(3,j,i)
          end do
        end do

      case default
        call output_line('Case not implemented!',&
            OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateByCollection')
        call sys_halt()
      end select
      
      ! Evaluate the 1st expression for the X-rhs
      allocate(DvaluesAct(npoints*nelements))
      call fparser_evalFunction (rsolution%rparserExpression, idim, 2, p_Dval, DvaluesAct)

      ! Reshape the data, that's it.
      do i=0,nelements-1
        do j=1,npoints
          Dvalues(j,i+1) = DvaluesAct(i*npoints+j)
        end do
      end do
      
      deallocate(DvaluesAct)
      deallocate(p_Dval)
      
    case (ANSOL_TP_ANALYTICAL)
      if (present(ierror)) ierror = -1
      !call output_line('Solution undefined! Caller must evaluate the solution! Use iid as hint!',&
      !    OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateDirect')
      !call sys_halt()
      
    case (ANSOL_TP_ZERO)
      ! Zero target function
      Dvalues(:,:) = 0.0_DP
      
    case(ANSOL_TP_MBSTATIONARYFILE,ANSOL_TP_MBSTATIONARY,&
          ANSOL_TP_MBNONSTATIONARYFILE,ANSOL_TP_MBNONSTATIONARY)
      ! Stationary/nonstationary function, specified by a block vector.
      !
      ! For every point, find the element of an element nearby the point.
      ! The evaluation routine uses this as hint to speed up the evaluation.
      allocate(DvaluesAct(npoints*nelements))
      NEL = rsolution%rstationary%p_rblockDiscr%p_rtriangulation%NEL

      ! Set up an array with point coordinates.
      allocate(DpointsAct(ubound(Dpoints,1),npoints*nelements))

      select case (ubound(Dpoints,1))
      case (NDIM2D)
      
        do i=0,nelements-1
          do j=1,npoints
            DpointsAct(1,i*npoints+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npoints+j) = Dpoints(2,j,i+1)
          end do
        end do
        
      case (NDIM3D)
            
        do i=0,nelements-1
          do j=1,npoints
            DpointsAct(1,i*npoints+j) = Dpoints(1,j,i+1)
            DpointsAct(2,i*npoints+j) = Dpoints(2,j,i+1)
            DpointsAct(3,i*npoints+j) = Dpoints(3,j,i+1)
          end do
        end do

      case default
        call output_line('Case not implemented!',&
            OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateDirect')
        call sys_halt()

      end select
      
      if (present(Ielements)) then
      
        ! Prepare an array that gives us a hint where to find the points
        ! in the mesh.
        allocate(IelementsAct(npoints*nelements))
        
        do i=0,nelements-1
          do j=1,npoints
            iel = Ielements(i+1)
            
            ! If the FE function is evaluated on a level higher than the
            ! discretisation allows, the element number may be out of bounds.
            ! Reduce it by dividing by 4 which simulates coarsening for a mesh
            ! which is refined by 2-level ordering.
            do while (iel .gt. NEL)
              iel = iel / 4
            end do
            
            IelementsAct(i*npoints+j) = iel
          end do
        end do
      
        ! Evaluate our temp vector at the given points.
        call fevl_evaluate (DER_FUNC, DvaluesAct, rsolution%rstationary%RvectorBlock(idim), &
          DpointsAct, IelementsHint=IelementsAct,cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
          
        deallocate(IelementsAct)
  
      else
  
        ! Evaluate our temp vector at the given points.
        call fevl_evaluate (DER_FUNC, DvaluesAct, rsolution%rstationary%RvectorBlock(idim), &
          DpointsAct, cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
  
      end if
        
      do i=0,nelements-1
        do j=1,npoints
          Dvalues(j,i+1) = DvaluesAct(i*npoints+j)
        end do
      end do
      
      deallocate(DpointsAct)
      deallocate(DvaluesAct)
      
    case default
      call output_line('Function not valid!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_evaluateByCollection')
      call sys_halt()
    end select
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunctionPrj (cderivative, rdiscretisation, &
      nelements, npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset, &
      Dvalues, rcollection)
  
  use fsystem
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! Local callback function for projection of analytical solutions
  ! into vectors.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN) :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>
  
!</subroutine>

    ! Evaluate the solution. rcollection%IquickAccess(1) is the component.
    call ansol_evaluateByCollection (rcollection,"SOL",&
        rcollection%IquickAccess(1),Dvalues,&
        npointsPerElement,nelements,Dpoints,rdomainIntSubset%p_Ielements)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine ffunctionPrj2 (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, &
                IdofsTest, rdomainIntSubset, &
                Dcoefficients, rcollection)
  
  use fsystem
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection
  
!<description>
  ! This subroutine is called during the vector assembly. It has to compute
  ! the coefficients in front of the terms of the linear form.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points and all the terms in the linear form
  ! the corresponding coefficients in front of the terms.
!</description>
  
!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
  
  ! The linear form which is currently to be evaluated:
  type(t_linearForm), intent(IN) :: rform
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN) :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN) :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  ! DIMENSION(dimension,npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(IN) :: Dpoints

  ! An array accepting the DOF`s on all elements test in the test space.
  ! DIMENSION(\#local DOF`s in test space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional :: rcollection
  
!</input>

!<output>
  ! A list of all coefficients in front of all terms in the linear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the linear form.
  real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
!</output>
  
!</subroutine>

    real(dp), dimension(:,:), allocatable :: Dvalues
    integer, dimension(3) :: Ibounds
    
    Ibounds = ubound(Dcoefficients)
    allocate(Dvalues(Ibounds(2),Ibounds(3)))
    
    call ffunctionPrj (DER_FUNC, rdiscretisation, &
        nelements, npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset, &
        Dvalues, rcollection)
    
    Dcoefficients(1,:,:) = Dvalues(:,:)
    
    deallocate(Dvalues)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_prjToVector (rsolution,dtime,rvector,idimstart,idimend,&
      rmassMatrix,idimdest,imethod)

!<description>
  ! Projects an analytically given solution to a discrete vector.
!</description>

!<input>
  ! The solution to evaluate.
  type(t_anSolution), intent(inout), target :: rsolution
  
  ! Time where to project.
  ! Can be =0 for a stationary solution.
  real(dp), intent(in) :: dtime
  
  ! First dimension to evaluate.
  ! =0: First dimension
  integer, intent(in) :: idimstart

  ! Last dimension to evaluate.
  ! =0: Last dimension
  integer, intent(in) :: idimend
  
  ! Block matrix with mass matrices on the diagonal for all components,
  ! corresponding to the destination vector.
  type(t_matrixBlock), intent(in) :: rmassMatrix
  
  ! Target dimension in rvector where to write the evaluated solution
  ! to. Only valid for idim > 0. If not specified, this defaults to idim.
  integer, intent(in), optional :: idimdest
  
  ! OPTIONAL: Method how to project.
  ! =0: simple evaluation (default)
  ! =1: L2-projection
  integer, intent(in), optional :: imethod
!</input>

!<inputoutput>
  ! Target vector that receives the FEM representation of rsolution.
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>
  
!</subroutine>
    
    ! local variables
    type(t_collection) :: rcollection
    integer :: imeth,istart,iend,idim2, i
    type(t_configL2ProjectionByMass) :: rL2ProjectionConfig
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Prepare a collection
    call collct_init (rcollection)
    
    ! Prepare the collection for evaluation.
    call ansol_prepareEvalCollection (rsolution,rcollection,"SOL",dtime)
    istart = 1
    iend = max(rsolution%ncomponents,rvector%nblocks)
    if (idimstart .ne. 0) istart = idimstart
    if (idimend .ne. 0) iend = idimend

    imeth = 0
    if (present(imethod)) imeth = imethod
    
    select case (imeth)
    case (0)
      ! Direct evaluation.
      !
      ! Loop through the dimensions to evaluate
      do i = istart,iend
        
        ! Set up destination dimension relative to idest
        idim2 = i
        if (present(idimdest)) &
          idim2 = idimdest + i - istart
        
        ! Prepare the collection
        rcollection%IquickAccess(1) = i
        
        ! DEBUG!!!
        call lsyssc_getbase_double (rvector%RvectorBlock(idim2),p_Ddata)
        
        ! Evaluate
        !call anprj_discrDirect (rvector%RvectorBlock(idim2),ffunctionPrj,rcollection)
        rL2ProjectionConfig%depsRel = 1E-20_DP
        call anprj_analytL2projectionByMass (rvector%RvectorBlock(idim2),&
            rmassMatrix%RmatrixBlock(idim2,idim2),ffunctionPrj2,rcollection,&
            rL2ProjectionConfig)
        
      end do
    
    case (1)
      call output_line('L2 projection currently not supported!',&
          OU_CLASS_ERROR, OU_MODE_STD,'ansol_prjToVector')
      call sys_halt()
    end select
    
    ! Release data
    call ansol_doneEvalCollection (rcollection,"SOL")
    
    call collct_done (rcollection)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_prjToSpaceTimeVector (rsolution,rvector,idimstart,idimend,rmassmatrix,imethod)

!<description>
  ! Projects an analytically given solution to a discrete spac-time vector.
!</description>

!<input>
  ! The solution to evaluate.
  type(t_anSolution), intent(inout), target :: rsolution
  
  ! First spatial dimension to evaluate.
  ! =0: First dimension
  integer, intent(in) :: idimstart

  ! Last spatial dimension to evaluate.
  ! =0: Last dimension
  integer, intent(in) :: idimend

  ! Block matrix with mass matrices on the diagonal for all components,
  ! corresponding to the destination vector.
  type(t_matrixBlock), intent(in) :: rmassMatrix
  
  ! OPTIONAL: Method how to project.
  ! =0: simple evaluation (default)
  ! =1: L2-projection
  integer, intent(in), optional :: imethod
!</input>

!<inputoutput>
  ! Target vector that receives the FEM representation of rsolution.
  type(t_spacetimeVector), intent(inout) :: rvector
!</inputoutput>
  
!</subroutine>
    
    integer :: iinterval
    real(dp) :: dtime
    type(t_vectorBlock) :: rvectemp

    ! Create a temp vector
    call lsysbl_createVectorBlock (rvector%p_rspaceDiscr,rvectemp,.true.)

    ! Loop through all timesteps and calculate the projection there.
    do iinterval = 1,rvector%NEQtime
      ! For an interval value out of bounds, getTimeStep returns the start or end
      ! time, so the following is ok:
      call tdiscr_getTimestep(rvector%p_rtimediscr,iinterval,dtime)
      call ansol_prjToVector (rsolution,dtime,rvectemp,idimstart,idimend,&
          rmassMatrix,imethod=imethod)
      call sptivec_setTimestepData (rvector, iinterval, rvectemp)
    end do
    
    call lsysbl_releaseVector (rvectemp)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine ansol_getSpaceTimeDiscr (rsolution,p_rspaceDiscr,p_rtimeDiscr)

!<description>
  ! If the analytical solution is connected to a space and time discretisation,
  ! this routine returns pointers to these underlying discretisation structures.
  ! If not, NULL() is returned.
!</description>

!<input>
  ! The solution.
  type(t_anSolution), intent(inout), target :: rsolution
!</input>

!<output>
  ! Pointer to the underlying space discretisation or NULL()
  type(t_blockDiscretisation), pointer :: p_rspaceDiscr
  
  ! Pointer to the underlying time discretisation or NULL()
  type(t_timeDiscretisation), pointer :: p_rtimeDiscr
  
!</output>
  
!</subroutine>
    
    nullify(p_rtimeDiscr)
    nullify(p_rspaceDiscr)

    if ((rsolution%ctype .eq. ANSOL_TP_MBNONSTATIONARY) .or. &
        (rsolution%ctype .eq. ANSOL_TP_MBNONSTATIONARYFILE)) then
      p_rtimeDiscr => rsolution%rtimeDiscr
      p_rspaceDiscr => rsolution%rfeSpace%p_rdiscretisation
    end if

  end subroutine

end module

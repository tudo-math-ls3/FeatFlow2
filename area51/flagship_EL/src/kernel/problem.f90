!##############################################################################
!# ****************************************************************************
!# <name> problem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures and subroutines for handling
!#  the complete problem configuration both in space and time.
!#
!#
!# The following routines are available:
!#
!# 1.) problem_initProblem
!#     -> Initializes a new problem structure
!#
!# 2.) problem_createProblem
!#     -> Creates a new problem structure
!#
!# 3.) problem_releaseProblem
!#     -> Releases an existing problem structure
!#
!# 4.) problem_done
!#     -> Release an existing sequence of problem structures
!#
!# 5.) problem_appendProblem
!#     -> Appends a problem structure to the linked list of problem structures
!#
!# 6.) problem_prependProblem
!#     -> Prepends a problem structure to the linked list of problem structures
!#
!# 7.) problem_removeProblem
!#     -> Removes a problem structure from the linked list of problem structures
!#
!# 8.) problem_infoProblem
!#     -> Outputs information about the problem structure
!#
!# 9.) problem_createLevel
!#     -> Creates a new problem level structure
!#
!# 10.) problem_releaseLevel
!#      -> Releases an existing problem level structure
!#
!# 11.) problem_appendLevel
!#     -> Appends a problem level structure into the linked list of problem
!#        level structures
!#
!# 12.) problem_prependLevel
!#     -> Prepends a problem level structure into the linked list of problem
!#        level structures
!#
!# 13.) problem_removeLevel
!#      -> Removes an existing problem level structure from the linked list of
!#         problem level structures
!#
!# 14.) problem_infoLevel
!#      -> Outputs information about the problem level structure
!#
!# 15.) problem_getLevel
!#      -> Returns a pointer to the problem level specified by the level number
!#
!# </purpose>
!##############################################################################

module problem

  use afcstabilisation
  use basicgeometry
  use boundary
  use boundaryfilter
  use fparser
  use fsystem
  use genoutput
  use io
  use linearsystemblock
  use linearsystemscalar
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private
  public :: t_problem
  public :: t_problemLevel
  public :: t_problemDescriptor
  public :: problem_initProblem
  public :: problem_createProblem
  public :: problem_releaseProblem
  public :: problem_done
  public :: problem_appendProblem
  public :: problem_prependProblem
  public :: problem_removeProblem
  public :: problem_infoProblem
  public :: problem_createLevel
  public :: problem_releaseLevel
  public :: problem_appendLevel
  public :: problem_prependLevel
  public :: problem_removeLevel
  public :: problem_infoLevel
  public :: problem_getLevel

  !*****************************************************************************

!<constants>

!<constantblock description="Flags for the problem descriptor specification bitfield">

  ! Standard problem
  integer(I32), parameter, public :: PROBDESC_MSPEC_STANDARD         = 0

  ! Generate extended raw mesh
  integer(I32), parameter, public :: PROBDESC_MSPEC_EXTENDEDRAW      = 2**0

  ! Convert quadrilaterals to  triangles
  integer(I32), parameter, public :: PROBDESC_MSPEC_CONVTRIANGLES    = 2**1

  ! Convert hexahedrals to tetrahedrals
  integer(I32), parameter, public :: PROBDESC_MSPEC_CONVTETRAHEDRALS = 2**2

!</constantblock>


!<constantblock description="Flags for the problem level specification bitfield">

  ! Standard problem level
  integer(I32), parameter, public :: PROBLEV_MSPEC_STANDARD         = 0

  ! Problem level requires update
  integer(I32), parameter, public :: PROBLEV_MSPEC_UPDATE           = 2**0

  ! Problem level requires initialisation
  integer(I32), parameter, public :: PROBLEV_MSPEC_INITIALIZE       = 2**1

!</constantblock>

!</constants>

  !*****************************************************************************

!<types>

!<typeblock>

  ! This data structure contains the abstract problem description

  type t_problemDescriptor

    ! Problem descriptor specification tag. This is a bitfield coming
    ! from an OR combination of different PROBDESC_MSPEC_xxxx
    ! constants and specifies various details of the problem
    ! descriptor. If it is =PROBDESC_MSPEC_STANDARD, the problem
    ! descriptor is a usual descriptor that needs no special handling.
    integer(I32) :: iproblemSpec = PROBDESC_MSPEC_STANDARD

    ! Spatial dimension
    integer :: ndimension

    ! Maximum problem level
    integer :: nlmax

    ! Minimum problem level
    integer :: nlmin

    ! Number of discretisations
    integer :: ndiscretisation

    ! Number of AFC stabilisations
    integer :: nafcstab

    ! Number of scalar matrices
    integer :: nmatrixScalar

    ! Number of block matrices
    integer :: nmatrixBlock

    ! Number of scalar vectors
    integer :: nvectorScalar

    ! Number of block vectors
    integer :: nvectorBlock

    ! Name of the triangulation file
    character(LEN=SYS_STRLEN) :: trifile

    ! Name of the boundary parametrisation file
    character(LEN=SYS_STRLEN) :: prmfile

  end type t_problemDescriptor

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure contains the complete problem configuration

  type t_problem

    ! Boundary parametrisation
    type(t_boundary) :: rboundary

    ! Pointer to the previous problem instance
    type(t_problem), pointer :: p_rproblemPrev => null()

    ! Pointer to the next problem instance
    type(t_problem), pointer :: p_rproblemNext => null()

    ! Pointers to the minimum problem level
    type(t_problemLevel), pointer :: p_rproblemLevelMin => null()

    ! Pointers to the maximum problem level
    type(t_problemLevel), pointer :: p_rproblemLevelMax => null()

  end type t_problem

!</typeblock>

  !*****************************************************************************

!<typeblock>

  ! This data structure contains all required problem data on one level

  type t_problemLevel

    ! Problem level specification tag. This is a bitfield coming from
    ! an OR combination of different PROBLEV_MSPEC_xxxx constants and
    ! specifies variour details of the problem level. If it is
    ! =PROBLEV_MSPEC_INITIALIZE, the problem level is a usual problem
    ! level taht needs no special handling.
    integer(I32) :: iproblemSpec = PROBLEV_MSPEC_INITIALIZE

    ! Number of the problem level
    integer :: ilev

    ! Triangulation structure
    type(t_triangulation) :: rtriangulation

    ! Array of discretisation structure
    type(t_blockDiscretisation), dimension(:), pointer :: Rdiscretisation => null()

    ! Array of AFC stabilisations
    type(t_afcstab), dimension(:), pointer :: Rafcstab => null()

    ! Array of scalar matrices
    type(t_matrixScalar), dimension(:), pointer :: Rmatrix => null()

    ! Array of block matrices
    type(t_matrixBlock), dimension(:), pointer :: RmatrixBlock => null()

    ! Array of scalar vectors
    type(t_vectorScalar), dimension(:), pointer :: Rvector => null()

    ! Array of block vectors
    type(t_vectorBlock), dimension(:), pointer :: RvectorBlock => null()

    ! Pointer to the global problem structure
    type(t_problem), pointer :: p_rproblem => null()

    ! Pointer to next coarse problem level
    type(t_problemLevel), pointer :: p_rproblemLevelCoarse => null()

    ! Pointer to next finer problem level
    type(t_problemLevel), pointer :: p_rproblemLevelFine => null()

  end type t_problemLevel

!</typeblock>

!</types>

contains

  !*****************************************************************************

!<subroutine>

  subroutine problem_initProblem(rproblemDescriptor, rproblem)

!<description>
    ! This subroutine initializes an abstract problem structure
!</description>

!<input>
    ! abstract problem descriptor
    type(t_problemDescriptor), intent(in) :: rproblemDescriptor
!</input>

!<output>
    ! problem data structure
    type(t_problem), intent(out) :: rproblem
!</output>
!</subroutine>

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel, p_rproblemLevel
    integer :: ilev
    logical :: bnoExtendedRaw

    ! Initialize global problem structure
    call problem_createProblem(rproblem)

    ! Initialize coarse level
    nullify(rproblemLevel); allocate(rproblemLevel)
    call problem_createLevel(rproblemLevel, rproblemDescriptor%nlmin)

    bnoExtendedRaw = (iand(rproblemDescriptor%iproblemSpec,&
                           PROBDESC_MSPEC_EXTENDEDRAW) .eq. 0)

    select case(rproblemDescriptor%ndimension)
    case (NDIM1D)
      ! Read coarse mesh from TRI-file
      call tria_readTriFile1D(rproblemLevel%rtriangulation,&
          rproblemDescriptor%trifile, bnoExtendedRaw)

    case (NDIM2D)
      ! Create new boundary and read from PRM-file
      call boundary_read_prm(rproblem%rboundary,&
          rproblemDescriptor%prmfile)

      ! Read coarse mesh from TRI-file
      call tria_readTriFile2D(rproblemLevel%rtriangulation,&
          rproblemDescriptor%trifile, rproblem%rboundary,&
          bnoextendedRaw)

      ! Convert quadrilaterals to triangules if required
      if (iand(rproblemDescriptor%iproblemSpec,&
               PROBDESC_MSPEC_CONVTRIANGLES) .ne. 0) then
        call tria_rawGridToTri(rproblemLevel%rtriangulation)
      end if

    case (NDIM3D)
      ! Read coarse mesh from TRI-file
      call tria_readTriFile3D(rproblemLevel%rtriangulation,&
          rproblemDescriptor%trifile, rproblem%rboundary,&
          bnoextendedRaw)

      ! Convert hexahedrals to tetrahedrals if required
      if (iand(rproblemDescriptor%iproblemSpec,&
               PROBDESC_MSPEC_CONVTETRAHEDRALS) .ne. 0) then
        call output_line('Conversion to tetrahedrals is not available yet!',&
            OU_CLASS_WARNING,OU_MODE_STD,'problem_initProblem')
        call sys_halt()
      end if

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_WARNING,OU_MODE_STD,'problem_initProblem')
      call sys_halt()
    end select

    ! Refine coarse mesh to minimum problem level
    call tria_quickRefine2LevelOrdering(rproblemDescriptor%nlmin-1,&
        rproblemLevel%rtriangulation, rproblem%rboundary)

    ! Create standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
        rproblem%rboundary)

    ! Allocate matrices, vectors, stabilisations, etc.
    if (rproblemDescriptor%ndiscretisation .gt. 0)&
        allocate(rproblemLevel%Rdiscretisation(&
        rproblemDescriptor%ndiscretisation))
    if (rproblemDescriptor%nmatrixScalar .gt. 0)&
        allocate(rproblemLevel%Rmatrix(&
        rproblemDescriptor%nmatrixScalar))
    if (rproblemDescriptor%nmatrixBlock .gt. 0)&
        allocate(rproblemLevel%RmatrixBlock(&
        rproblemDescriptor%nmatrixBlock))
    if (rproblemDescriptor%nvectorScalar .gt. 0)&
        allocate(rproblemLevel%Rvector(&
        rproblemDescriptor%nvectorScalar))
    if (rproblemDescriptor%nvectorBlock .gt. 0)&
        allocate(rproblemLevel%RvectorBlock(&
        rproblemDescriptor%nvectorBlock))
    if (rproblemDescriptor%nafcstab .gt. 0)&
        allocate(rproblemLevel%Rafcstab(&
        rproblemDescriptor%nafcstab))

    ! Append level to global problem
    call problem_appendLevel(rproblem, rproblemLevel)
    p_rproblemLevel => rproblemLevel

    ! Generate fine levels
    do ilev = rproblemDescriptor%nlmin+1,&
              rproblemDescriptor%nlmax

      ! Initialize current level
      nullify(rproblemLevel); allocate(rproblemLevel)
      call problem_createLevel(rproblemLevel, ilev)

      ! Generate regularly refined mesh
      call tria_refine2LevelOrdering(p_rproblemLevel%rtriangulation,&
          rproblemLevel%rtriangulation, rproblem%rboundary)

      ! Create standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
          rproblem%rboundary)

      ! Allocate matrices, vectors and stabilisations
      if (rproblemDescriptor%ndiscretisation .gt. 0)&
        allocate(rproblemLevel%Rdiscretisation(&
        rproblemDescriptor%ndiscretisation))
      if (rproblemDescriptor%nmatrixScalar .gt. 0)&
          allocate(rproblemLevel%Rmatrix(&
          rproblemDescriptor%nmatrixScalar))
      if (rproblemDescriptor%nmatrixBlock .gt. 0)&
          allocate(rproblemLevel%RmatrixBlock(&
          rproblemDescriptor%nmatrixBlock))
      if (rproblemDescriptor%nvectorScalar .gt. 0)&
          allocate(rproblemLevel%Rvector(&
          rproblemDescriptor%nvectorScalar))
      if (rproblemDescriptor%nvectorBlock .gt. 0)&
          allocate(rproblemLevel%RvectorBlock(&
          rproblemDescriptor%nvectorBlock))
      if (rproblemDescriptor%nafcstab .gt. 0)&
          allocate(rproblemLevel%Rafcstab(&
          rproblemDescriptor%nafcstab))

      ! Append current level to global problem
      call problem_appendLevel(rproblem, rproblemLevel)
      p_rproblemLevel => rproblemLevel
    end do

    ! Compress triangulation structure
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while (associated(p_rproblemLevel) .and.&
              associated(p_rproblemLevel%p_rproblemLevelCoarse))
      call tria_compress2LevelOrdHierarchy(&
          p_rproblemLevel%rtriangulation,&
          p_rproblemLevel%p_rproblemLevelCoarse%rtriangulation)
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
    end do

  end subroutine problem_initProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_createProblem(rproblem)

!<description>
    ! This subroutine creates a new problem structure
!</description>

!<intputoutput>
    ! problem data structure
    type(t_problem), intent(out) :: rproblem
!</inputoutput>
!</subroutine>

    ! Reset data
    nullify(rproblem%p_rproblemPrev, rproblem%p_rproblemNext)
    nullify(rproblem%p_rproblemLevelMin, rproblem%p_rproblemLevelMax)

  end subroutine problem_createProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_releaseProblem(rproblem)

!<description>
    ! This subroutine releases an existing problem structure
!</description>

!<inputoutput>
    ! problem data structure
    type(t_problem), intent(inout) :: rproblem
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! Initialisation
    p_rproblemLevel => rproblem%p_rproblemLevelMax

    ! Loop over all problem levels and destroy them
    do while(associated(p_rproblemLevel))
      call problem_removeLevel(rproblem, p_rproblemLevel)
      call problem_releaseLevel(p_rproblemLevel)

      deallocate(p_rproblemLevel)
      p_rproblemLevel => rproblem%p_rproblemLevelMax
    end do

    ! Release boundary
    call boundary_release(rproblem%rboundary)

    ! Reset data
    nullify(rproblem%p_rproblemPrev, rproblem%p_rproblemNext)
    nullify(rproblem%p_rproblemLevelMin, rproblem%p_rproblemLevelMax)

  end subroutine problem_releaseProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_done(p_rproblemFirst, p_rproblemLast)

!<description>
    ! This subroutine releases a sequence of existing problem structures
!</description>

!<inputoutput>
    ! first problem structure
    type(t_problem), pointer :: p_rproblemFirst

    ! last problem structure
    type(t_problem), pointer :: p_rproblemLast
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_problem), pointer :: p_rproblem


    p_rproblem => p_rproblemFirst
    do while(associated(p_rproblem))

      call problem_removeProblem(p_rproblemFirst,&
          p_rproblemLast, p_rproblem)
      call problem_releaseProblem(p_rproblem)
      deallocate(p_rproblem)
      p_rproblem => p_rproblemFirst
    end do

  end subroutine problem_done

  !*****************************************************************************

!<subroutine>

  subroutine problem_appendProblem(p_rproblemFirst, p_rproblemLast,&
      rproblem, rproblemRef)

!<description>
    ! This subroutine appends a problem structure to the linked list
    ! of problem structures. If the optional reference problem
    ! rproblemRef is given, the problem structure is appended to
    ! rproblemRef. Otherwise, the last problem is used as reference
    ! structure.
!</description>

!<input>
    ! first problem structure
    type(t_problem), pointer :: p_rproblemFirst

    ! last problem structure
    type(t_problem), pointer :: p_rproblemLast
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! OPTIONAL: reference problem structure
    type(t_problem), intent(inout), target, optional :: rproblemRef
!</inputoutput>
!</subroutine>

    if (associated(p_rproblemFirst) .and.&
        associated(p_rproblemLast)) then

      if (present(rproblemRef)) then

        ! Insert rproblem after rproblemRef
        rproblem%p_rproblemPrev => rproblemRef
        rproblem%p_rproblemNext => rproblemRef%p_rproblemNext

        if (associated(rproblemRef%p_rproblemNext)) then
          rproblemRef%p_rproblemNext%p_rproblemPrev => rproblem
        else
          p_rproblemLast => rproblem
        end if

        rproblemRef%p_rproblemNext => rproblem

      else

        ! Set pointer to last problem structure
        rproblem%p_rproblemPrev => p_rproblemLast
        nullify(rproblem%p_rproblemNext)

        p_rproblemLast%p_rproblemNext => rproblem
        p_rproblemLast => rproblem

      end if

    else

      ! List of problem structures is completely empty
      p_rproblemFirst => rproblem
      p_rproblemLast  => rproblem

      nullify(rproblem%p_rproblemPrev, rproblem%p_rproblemNext)

    end if

  end subroutine problem_appendProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_prependProblem(p_rproblemFirst, p_rproblemLast,&
      rproblem, rproblemRef)

!<description>
    ! This subroutine prepends a problem structure to the linked list
    ! of problem structures. If the optional reference problem
    ! rproblemRef is given, the problem structure is prepended to
    ! rproblemRef. Otherwise, the first problem is used as reference
    ! structure.
!</description>

!<input>
    ! first problem structure
    type(t_problem), pointer :: p_rproblemFirst

    ! last problem structure
    type(t_problem), pointer :: p_rproblemLast
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! OPTIONAL: reference problem structure
    type(t_problem), intent(inout), target, optional :: rproblemRef
!</inputoutput>
!</subroutine>

    if (associated(p_rproblemFirst) .and.&
        associated(p_rproblemLast)) then

      if (present(rproblemRef)) then

        ! Insert rproblem before rproblemRef
        rproblem%p_rproblemPrev => rproblemRef%p_rproblemPrev
        rproblem%p_rproblemNext => rproblemRef

        if (associated(rproblemRef%p_rproblemPrev)) then
          rproblemRef%p_rproblemPrev%p_rproblemNext => rproblem
        else
          p_rproblemFirst => rproblem
        end if

        rproblemRef%p_rproblemPrev => rproblem

      else

        ! Set pointer to first problem structure
        rproblem%p_rproblemNext => p_rproblemFirst
        nullify(rproblem%p_rproblemPrev)

        p_rproblemFirst%p_rproblemPrev => rproblem
        p_rproblemFirst => rproblem

      end if

    else

      ! List of problem structures is completely empty
      p_rproblemFirst => rproblem
      p_rproblemLast  => rproblem

      nullify(rproblem%p_rproblemPrev, rproblem%p_rproblemNext)

    end if

  end subroutine problem_prependProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_removeProblem(p_rproblemFirst, p_rproblemLast,&
      rproblem)

!<description>
    ! This subroutine removes an existing problem structure
    ! from the linked list of problem structures
!</description>

!<input>
    ! first problem structure
    type(t_problem), pointer :: p_rproblemFirst

    ! last problem structure
    type(t_problem), pointer :: p_rproblemLast
!</input>

!<inputoutput>
    ! problem structure
    type(t_problem), intent(inout) :: rproblem
!</inputoutput>
!</subroutine>

    if (.not.associated(rproblem%p_rproblemPrev)) then
      p_rproblemFirst => rproblem%p_rproblemNext
    else
      rproblem%p_rproblemPrev%p_rproblemNext => rproblem%p_rproblemNext
    end if

    if (.not.associated(rproblem%p_rproblemNext)) then
      p_rproblemLast => rproblem%p_rproblemPrev
    else
      rproblem%p_rproblemNext%p_rproblemPrev => rproblem%p_rproblemPrev
    end if

  end subroutine problem_removeProblem

  !*****************************************************************************

!<subroutine>

    subroutine problem_infoProblem(rproblemFirst, rproblemLast)

!<description>
    ! This subroutine outputs information about the problem structure
!</description>

!<input>
    ! first problem structure
    type(t_problem), intent(in), target :: rproblemFirst

    ! OPTIONAL: last problem structure
    type(t_problem), intent(in), target, optional :: rproblemLast
!</input>
!</subroutine>

    ! local variables
    type(t_problem), pointer :: p_rproblem


    if (present(rproblemLast)) then

      p_rproblem => rproblemFirst

      do while(associated(p_rproblem))

        call doInfo(p_rproblem)
        call output_lbrk()

        if (associated(p_rproblem, rproblemLast)) exit
        p_rproblem => p_rproblem%p_rproblemNext

      end do

    else

      call doInfo(rproblemFirst)

    end if

  contains

    ! Here, the real info routines follow.

    !**************************************************************

    subroutine doInfo(rproblem)

      type(t_problem), intent(in), target :: rproblem

      ! local variables
      type(t_problemLevel), pointer :: p_rproblemLevel


      call output_line ('Problem:')
      call output_line ('---------')

      if (associated(rproblem%p_rproblemLevelMin)) then
        call output_line ('minimum level: '//&
            trim(sys_siL(rproblem%p_rproblemLevelMin%ilev,15)))
      else
        call output_line ('minimum level: not associated')
      end if

      if (associated(rproblem%p_rproblemLevelMax)) then
        call output_line ('maximum level: '//&
            trim(sys_siL(rproblem%p_rproblemLevelMax%ilev,15)))
      else
        call output_line ('maximum level: not associated')
      end if

      call output_lbrk()

      ! Initialisation
      p_rproblemLevel => rproblem%p_rproblemLevelMax

      ! Loop over all problem levels
      do while(associated(p_rproblemLevel))
        call problem_infoLevel(p_rproblemLevel)
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

    end subroutine doInfo

  end subroutine problem_infoProblem

  !*****************************************************************************

!<subroutine>

  subroutine problem_createLevel(rproblemLevel, ilev)

!<description>
    ! This subroutine creates a new problem level structure
!</description>

!<input>
    ! level number
    integer, intent(in) :: ilev
!</input>

!<output>
    ! problem level structure
    type(t_problemLevel), intent(out) :: rproblemLevel
!</output>
!</subroutine>

    ! Set problem level
    rproblemLevel%ilev = ilev

  end subroutine problem_createLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_releaseLevel(rproblemLevel)

!<description>
    ! This subroutine releases an existing problem level structure
!</description>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Release triangulation structure
    call tria_done(rproblemLevel%rtriangulation)

    ! Release discretisation structures
    if (associated(rproblemLevel%Rdiscretisation)) then
      do i = lbound(rproblemLevel%Rdiscretisation,1),&
             ubound(rproblemLevel%Rdiscretisation,1)
        call spdiscr_releaseBlockDiscr(rproblemLevel%Rdiscretisation(i))
      end do
      deallocate(rproblemLevel%Rdiscretisation)
    end if

    ! Release all scalar matrices
    if (associated(rproblemLevel%Rmatrix)) then
      do i = lbound(rproblemLevel%Rmatrix,1),&
             ubound(rproblemLevel%Rmatrix,1)
        call lsyssc_releaseMatrix(rproblemLevel%Rmatrix(i))
      end do
      deallocate(rproblemLevel%Rmatrix)
    end if

    ! Release all block matries
    if (associated(rproblemLevel%RmatrixBlock)) then
      do i = lbound(rproblemLevel%RmatrixBlock,1),&
             ubound(rproblemLevel%RmatrixBlock,1)
        call lsysbl_releaseMatrix(rproblemLevel%RmatrixBlock(i))
      end do
      deallocate(rproblemLevel%RmatrixBlock)
    end if

    ! Release all scalar vectors
    if (associated(rproblemLevel%Rvector)) then
      do i = lbound(rproblemLevel%Rvector,1),&
             ubound(rproblemLevel%Rvector,1)
        call lsyssc_releaseVector(rproblemLevel%Rvector(i))
      end do
      deallocate(rproblemLevel%Rvector)
    end if

    ! Release all block vectors
    if (associated(rproblemLevel%RvectorBlock)) then
      do i = lbound(rproblemLevel%RvectorBlock,1),&
             ubound(rproblemLevel%RvectorBlock,1)
        call lsysbl_releaseVector(rproblemLevel%RvectorBlock(i))
      end do
      deallocate(rproblemLevel%RvectorBlock)
    end if

    ! Release stabilisation structure
    if (associated(rproblemLevel%Rafcstab)) then
      do i = lbound(rproblemLevel%Rafcstab,1),&
             ubound(rproblemLevel%Rafcstab,1)
        call afcstab_releaseStabilisation(rproblemLevel%Rafcstab(i))
      end do
      deallocate(rproblemLevel%Rafcstab)
    end if

  end subroutine problem_releaseLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_appendLevel(rproblem, rproblemLevel,&
      rproblemLevelRef)

!<description>
    ! This subroutine appends a problem level structure to the linked
    ! list of problem level structures. If the optional reference
    ! level rproblemLevelRef is given, the problem level structure is
    ! appended to rproblemLevelRef. Otherwise, the maximum problem
    ! level is used as reference structure.
!</description>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! OPTIONAL: reference problem structure
    type(t_problemLevel), intent(inout), target, optional :: rproblemLevelRef
!</inputoutput>
!</subroutine>

    ! Set pointer to global problem structure
    rproblemLevel%p_rproblem => rproblem

    if (associated(rproblem%p_rproblemLevelMin) .and.&
        associated(rproblem%p_rproblemLevelMax)) then

      if (present(rproblemLevelRef)) then

        ! Insert rproblemLevel after rproblemLevelRef
        rproblemLevel%p_rproblemLevelCoarse => rproblemLevelRef
        rproblemLevel%p_rproblemLevelFine   => rproblemLevelRef%p_rproblemLevelFine

        if (associated(rproblemLevelRef%p_rproblemLevelFine)) then
          rproblemLevelRef%p_rproblemLevelFine%p_rproblemLevelCoarse => rproblemLevel
        else
          rproblem%p_rproblemLevelMax => rproblemLevel
        end if

        rproblemLevelRef%p_rproblemLevelFine => rproblemLevel

      else

        ! Set pointer to maximum problem level
        rproblemLevel%p_rproblemLevelCoarse => rproblem%p_rproblemLevelMax
        nullify(rproblemLevel%p_rproblemLevelFine)

        rproblem%p_rproblemLevelMax%p_rproblemLevelFine => rproblemLevel
        rproblem%p_rproblemLevelMax => rproblemLevel

      end if

    else

      ! Problem structure is completely empty
      rproblem%p_rproblemLevelMin => rproblemLevel
      rproblem%p_rproblemLevelMax => rproblemLevel

      nullify(rproblemLevel%p_rproblemLevelCoarse, rproblemLevel%p_rproblemLevelFine)

    end if
  end subroutine problem_appendLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_prependLevel(rproblem, rproblemLevel,&
      rproblemLevelRef)

!<description>
    ! This subroutine prepends a problem level structure to the linked
    ! list of problem level structures. If the optional reference
    ! level rproblemLevelRef is given, the problem level structure is
    ! prepended to rproblemLevelRef. Otherwise, the minimum problem
    ! level is used as reference structure.
!</description>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! OPTIONAL: reference problem structure
    type(t_problemLevel), intent(inout), target, optional :: rproblemLevelRef
!</inputoutput>
!</subroutine>

    ! Set pointer to global problem structure
    rproblemLevel%p_rproblem => rproblem

    if (associated(rproblem%p_rproblemLevelMin) .and.&
        associated(rproblem%p_rproblemLevelMax)) then

      if (present(rproblemLevelRef)) then

        ! Insert rproblemLevel before rproblemLevelRef
        rproblemLevel%p_rproblemLevelCoarse => rproblemLevelRef%p_rproblemLevelCoarse
        rproblemLevel%p_rproblemLevelFine   => rproblemLevelRef

        if (associated(rproblemLevelRef%p_rproblemLevelCoarse)) then
          rproblemLevelRef%p_rproblemLevelCoarse%p_rproblemLevelFine => rproblemLevel
        else
          rproblem%p_rproblemLevelMin => rproblemLevel
        end if

        rproblemLevelRef%p_rproblemLevelCoarse => rproblemLevel

      else

        ! Set pointer to minimum problem level
        rproblemLevel%p_rproblemLevelFine => rproblem%p_rproblemLevelMin
        nullify(rproblemLevel%p_rproblemLevelCoarse)

        rproblem%p_rproblemLevelMin%p_rproblemLevelCoarse => rproblemLevel
        rproblem%p_rproblemLevelMin => rproblemLevel

      end if

    else

      ! Problem structure is completely empty
      rproblem%p_rproblemLevelMin => rproblemLevel
      rproblem%p_rproblemLevelMax => rproblemLevel

      nullify(rproblemLevel%p_rproblemLevelCoarse, rproblemLevel%p_rproblemLevelFine)

    end if

  end subroutine problem_prependLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_removeLevel(rproblem, rproblemLevel)

!<description>
    ! This subroutine removes an existing problem level structure from
    ! an existing problem structure
!</description>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(inout), target :: rproblem

    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel
!</inputoutput>
!</subroutine>

    if (.not.associated(rproblemLevel%p_rproblem, rproblem)) then
      call output_line('Problem level structure does not belong to problem structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'problem_removeLevel')
      call sys_halt()
    end if

    if (.not.associated(rproblemLevel%p_rproblemLevelCoarse)) then
      rproblem%p_rproblemLevelMin => rproblemLevel%p_rproblemLevelFine
    else
      rproblemLevel%p_rproblemLevelCoarse%p_rproblemLevelFine => rproblemLevel%p_rproblemLevelFine
    end if

    if (.not.associated(rproblemLevel%p_rproblemLevelFine)) then
      rproblem%p_rproblemLevelMax => rproblemLevel%p_rproblemLevelCoarse
    else
      rproblemLevel%p_rproblemLevelFine%p_rproblemLevelCoarse => rproblemLevel%p_rproblemLevelCoarse
    end if

  end subroutine problem_removeLevel

  !*****************************************************************************

!<subroutine>

  subroutine problem_infoLevel(rproblemLevel)

!<description>
    ! This subroutine outputs information about the problem level structure
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel
!</input>
!</subroutine>

    call output_line ('ProblemLevel: '//trim(sys_siL(rproblemLevel%ilev,15)))
    call output_line ('-------------')

    if (associated(rproblemLevel%p_rproblemLevelCoarse)) then
      call output_line ('coarser level:       '//&
                        trim(sys_siL(rproblemLevel%p_rproblemLevelCoarse%ilev,15)))
    else
      call output_line ('coarser level:       not associated')
    end if

    if (associated(rproblemLevel%p_rproblemLevelFine)) then
      call output_line ('finer level:         '//&
                        trim(sys_siL(rproblemLevel%p_rproblemLevelFine%ilev,15)))
    else
      call output_line ('finer level:         not associated')
    end if

    if (associated(rproblemLevel%Rafcstab)) then
      call output_line ('Rafcstab:            '//&
                        trim(sys_siL(size(rproblemLevel%Rafcstab),15)))
    else
      call output_line ('Rafcstab:            not associated')
    end if

    if (associated(rproblemLevel%Rmatrix)) then
      call output_line ('Rmatrix:             '//&
                        trim(sys_siL(size(rproblemLevel%Rmatrix),15)))
    else
      call output_line ('Rmatrix:             not associated')
    end if

    if (associated(rproblemLevel%RmatrixBlock)) then
      call output_line ('RmatrixBlock:        '//&
                        trim(sys_siL(size(rproblemLevel%RmatrixBlock),15)))
    else
      call output_line ('RmatrixBlock:        not associated')
    end if

    if (associated(rproblemLevel%Rvector)) then
      call output_line ('Rvector:             '//&
                        trim(sys_siL(size(rproblemLevel%Rvector),15)))
    else
      call output_line ('Rvector:             not associated')
    end if

    if (associated(rproblemLevel%RvectorBlock)) then
      call output_line ('RvectorBlock:        '//&
                        trim(sys_siL(size(rproblemLevel%RvectorBlock),15)))
    else
      call output_line ('RvectorBlock:        not associated')
    end if
    call output_lbrk()

  end subroutine problem_infoLevel

  !*****************************************************************************

!<function>

  function problem_getLevel(rproblem, ilev, btopdown)&
      result(p_rproblemLevel)

!<description>
    ! This subroutine returns a pointer to the problem level structure
    ! with specified level number ilve. If such problem level does not
    ! exist p_rproblemLevel points to null().
!</description>

!<input>
    ! global problem structure
    type(t_problem), intent(in) :: rproblem

    ! level number
    integer, intent(in) :: ilev

    ! OPTIONAL: search direction
    ! If this parameter is .false., the loop starts at the minimum
    ! problem level and proceeds to the next finer problem
    ! level. Otherwise, the loop starts at the maximum problem level
    ! and continues with the next coarser level. The default value is
    ! .true., that is, top-down search is performed
    logical, intent(in), optional :: btopdown
!</input>

!<result
    ! problem level structure
    type(t_problemLevel), pointer :: p_rproblemLevel
!</result>
!</function>

    ! local variable
    logical :: bisTopdown

    bisTopdown = .true.
    if (present(btopdown)) bisTopdown=btopdown

    if (bisTopdown) then

      ! Top-down search
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      do while (associated(p_rproblemLevel))
        if(p_rproblemLevel%ilev .eq. ilev) return
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

    else

      ! Bottom-up search
      p_rproblemLevel => rproblem%p_rproblemLevelMin
      do while (associated(p_rproblemLevel))
        if(p_rproblemLevel%ilev .eq. ilev) return
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelFine
      end do

    end if

    nullify(p_rproblemLevel)

  end function problem_getLevel

end module problem

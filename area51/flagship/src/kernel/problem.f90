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
!# 1.) problem_createProblem
!#     -> Creates a new problem structure
!#
!# 2.) problem_releaseProblem
!#     -> Releases an existing problem structure
!#
!# 3.) problem_done
!#     -> Release an existing sequence of problem structures
!#
!# 4.) problem_appendProblem
!#     -> Appends a problem structure to the linked list of problem structures
!#
!# 5.) problem_prependProblem
!#     -> Prepends a problem structure to the linked list of problem structures
!#
!# 6.) problem_removeProblem
!#     -> Removes a problem structure from the linked list of problem structures
!#
!# 7.) problem_infoProblem
!#     -> Outputs information about the problem structure
!#
!# 8.) problem_createLevel
!#     -> Creates a new multigrid level structure
!#
!# 9.) problem_releaseLevel
!#     -> Releases an existing multigrid level structure
!#
!# 10.) problem_appendLevel
!#     -> Appends a multigrid level structure into the linked list of multigrid
!#        level structures
!#
!# 11.) problem_prependLevel
!#     -> Prepends a multigrid level structure into the linked list of multigrid
!#        level structures
!#
!# 12.) problem_removeLevel
!#      -> Removes an existing multigrid level structure from the linked list of
!#         multigrid level structures
!#
!# 13.) problem_infoLevel
!#      -> Outputs information about the multigrid level structure
!#
!# 14.) problem_createProfile
!#      -> Initializes vector profile from parameter file
!#
!# 15.) problem_info
!#      -> Outputs information about all problem structures in the linked list
!# </purpose>
!##############################################################################

module problem

  use afcstabilisation
  use boundary
  use fparser
  use io
  use linearsystemblock
  use linearsystemscalar
  use spatialdiscretisation
  use triangulation
  use fsystem

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
  public :: problem_createProfile
  
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

    ! Maximum multigrid level
    integer :: nlmax

    ! Minimum multigrid level
    integer :: nlmin

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
  !   
  ! 1. Each "problem" is allowed to exhibit a single boundary
  !    parametrization which is valid for all levels.
  ! 
  ! 2. Each "problem" is made up from multiple "levels". Starting from
  !    the coarsest level. Finer levels are generated by either
  !    regular refinement of the mesh or local mesh adaptivity. Beside
  !    the collection of meshes, each level provides memory access to
  !    matrices, vectors, index vectors, etc.  The levels are stored
  !    as a double-linked list which allows each level to access the
  !    next coarser and finer level directly. At the same time, linked
  !    lists provide the flexibility to easily remove or insert new
  !    levels without reorganizing the complete data structure.
  !
  ! 3. Each "problem" provides pointers to the previous and the next
  !    problem structure. This double-linked list allows to easily
  !    move through the temporal sequence of problem instances which
  !    may be useful if primal and dual problems for time-dependent
  !    flows are computed forward and backward in time.
  type t_problem

    ! Boundary parametrization
    type(t_boundary) :: rboundary

    ! Pointer to the previous problem instance
    type(t_problem), pointer :: p_rproblemPrev => null()

    ! Pointer to the next problem instance
    type(t_problem), pointer :: p_rproblemNext => null()

    ! Pointers to the minimum multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevelMin => null()

    ! Pointers to the maximum multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevelMax => null()
  end type t_problem

!</typeblock>

  !*****************************************************************************
  
!<typeblock>
  
  ! This data structure contains all required problem data on one level
  !
  ! 1. Each "level" provides a set of meshes which represent the 
  !    geometric basement of each finite element discretization.
  ! 2. Each "level" provides a set of matrix structures. Each matrix
  !    structure contains multiple data arrays which are build upon the
  !    underlying matrix structure. 
  type t_problemLevel

    ! Problem level specification tag. This is a bitfield coming from
    ! an OR combination of different PROBLEV_MSPEC_xxxx constants and
    ! specifies variour details of the problem level. If it is
    ! =PROBLEV_MSPEC_STANDARD, the problem level is a usual problem
    ! level taht needs no special handling.
    integer(I32) :: iproblemSpec = PROBLEV_MSPEC_STANDARD

    ! Number of the multigrid level
    integer :: ilev

    ! Triangulation structure
    type(t_triangulation) :: rtriangulation

    ! Discretization structure
    type(t_blockDiscretisation) :: rdiscretisation

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

    ! Pointer to next coarse multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevelCoarse => null()

    ! Pointer to next finer multigrid level
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
    type(t_problemDescriptor), intent(IN) :: rproblemDescriptor
!</input>

!<output>
    ! problem data structure
    type(t_problem), intent(OUT) :: rproblem
!</output>
!</subroutine>

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevel, p_rproblemLevel
    integer :: ilev
    logical :: berror, bnoExtendedRaw

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
      call boundary_read_prm(rproblem%rboundary, rproblemDescriptor%prmfile)
      
      ! Read coarse mesh from TRI-file
      call tria_readTriFile2D(rproblemLevel%rtriangulation, rproblemDescriptor%trifile,&
                              rproblem%rboundary, bnoextendedRaw)

      ! Convert quadrilaterals to triangules if required
      if (iand(rproblemDescriptor%iproblemSpec, PROBDESC_MSPEC_CONVTRIANGLES) .ne. 0)&
          call tria_rawGridToTri(rproblemLevel%rtriangulation)

    case (NDIM3D)
      ! Read coarse mesh from TRI-file
      call tria_readTriFile3D(rproblemLevel%rtriangulation, rproblemDescriptor%trifile,&
                              rproblem%rboundary, bnoextendedRaw)

      ! Convert hexahedrals to tetrahedrals if required
      if (iand(rproblemDescriptor%iproblemSpec, PROBDESC_MSPEC_CONVTETRAHEDRALS) .ne. 0) then
        print *, "This feature is not yet implemented"
        stop
      end if

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'problem_initProblem')
      call sys_halt()
    end select

    ! Refine coarse mesh to minimum multigrid level
    call tria_quickRefine2LevelOrdering(rproblemDescriptor%nlmin-1,&
                                        rproblemLevel%rtriangulation,&
                                        rproblem%rboundary)
    
    ! Create standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
                                      rproblem%rboundary)

    ! Allocate matrices, vectors and stabilisations
    if (rproblemDescriptor%nmatrixScalar .gt. 0)&
        allocate(rproblemLevel%Rmatrix(rproblemDescriptor%nmatrixScalar))
    if (rproblemDescriptor%nmatrixBlock .gt. 0)&
        allocate(rproblemLevel%RmatrixBlock(rproblemDescriptor%nmatrixBlock))
    if (rproblemDescriptor%nvectorScalar .gt. 0)&
        allocate(rproblemLevel%Rvector(rproblemDescriptor%nvectorScalar))
    if (rproblemDescriptor%nvectorBlock .gt. 0)&
        allocate(rproblemLevel%RvectorBlock(rproblemDescriptor%nvectorBlock))
    if (rproblemDescriptor%nafcstab .gt. 0)&
        allocate(rproblemLevel%Rafcstab(rproblemDescriptor%nafcstab))

    ! Append level to global problem
    call problem_appendLevel(rproblem, rproblemLevel)
    p_rproblemLevel => rproblemLevel

    ! Generate fine levels
    do ilev = rproblemDescriptor%nlmin+1, rproblemDescriptor%nlmax
      
      ! Initialize current level
      nullify(rproblemLevel); allocate(rproblemLevel)
      call problem_createLevel(rproblemLevel, ilev)
      
      ! Generate regularly refined mesh
      call tria_refine2LevelOrdering(p_rproblemLevel%rtriangulation,&
                                     rproblemLevel%rtriangulation,&
                                     rproblem%rboundary)

      ! Create standard mesh from raw mesh
      call tria_initStandardMeshFromRaw(rproblemLevel%rtriangulation,&
                                        rproblem%rboundary)

      ! Allocate matrices, vectors and stabilisations
      if (rproblemDescriptor%nmatrixScalar .gt. 0)&
          allocate(rproblemLevel%Rmatrix(rproblemDescriptor%nmatrixScalar))
      if (rproblemDescriptor%nmatrixBlock .gt. 0)&
          allocate(rproblemLevel%RmatrixBlock(rproblemDescriptor%nmatrixBlock))
      if (rproblemDescriptor%nvectorScalar .gt. 0)&
          allocate(rproblemLevel%Rvector(rproblemDescriptor%nvectorScalar))
      if (rproblemDescriptor%nvectorBlock .gt. 0)&
          allocate(rproblemLevel%RvectorBlock(rproblemDescriptor%nvectorBlock))
      if (rproblemDescriptor%nafcstab .gt. 0)&
          allocate(rproblemLevel%Rafcstab(rproblemDescriptor%nafcstab))
      
      ! Append current level to global problem
      call problem_appendLevel(rproblem, rproblemLevel)
      p_rproblemLevel => rproblemLevel
    end do
    
    ! Compress triangulation structure
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    do while (associated(p_rproblemLevel) .and.&
              associated(p_rproblemLevel%p_rproblemLevelCoarse))
      call tria_compress2LevelOrdHierarchy(p_rproblemLevel%rtriangulation,&
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
    type(t_problem), intent(OUT) :: rproblem
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
    type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! Initialization
    p_rproblemLevel => rproblem%p_rproblemLevelMax
    
    ! Loop over all multigrid levels and destroy them
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

      call problem_removeProblem(p_rproblemFirst, p_rproblemLast, p_rproblem)
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
    type(t_problem), intent(INOUT), target :: rproblem

    ! OPTIONAL: reference problem structure
    type(t_problem), intent(INOUT), target, optional :: rproblemRef
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
    type(t_problem), intent(INOUT), target :: rproblem

    ! OPTIONAL: reference problem structure
    type(t_problem), intent(INOUT), target, optional :: rproblemRef
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

  subroutine problem_removeProblem(p_rproblemFirst, p_rproblemLast, rproblem)

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
    type(t_problem), intent(INOUT) :: rproblem
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
    type(t_problem), intent(IN), target :: rproblemFirst

    ! OPTIONAL: last problem structure
    type(t_problem), intent(IN), target, optional :: rproblemLast
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

      type(t_problem), intent(IN), target :: rproblem

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
      
      ! Initialization
      p_rproblemLevel => rproblem%p_rproblemLevelMax
      
      ! Loop over all multigrid levels
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
    ! This subroutine creates a new multigrid level structure
!</description>

!<input>
    ! level number
    integer, intent(IN) :: ilev
!</input>

!<output>
    ! multigrid level structure
    type(t_problemLevel), intent(OUT) :: rproblemLevel
!</output>
!</subroutine>    
    
    ! Set multigrid level
    rproblemLevel%ilev = ilev
    
  end subroutine problem_createLevel
  
  !*****************************************************************************

!<subroutine>

  subroutine problem_releaseLevel(rproblemLevel)

!<description>
    ! This subroutine releases an existing multigrid level structure
!</description>
 
!<inputoutput>  
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i
    
    ! Release triangulation structure
    call tria_done(rproblemLevel%rtriangulation)

    ! Release discretization structure
    call spdiscr_releaseBlockDiscr(rproblemLevel%rdiscretisation)
    
    ! Release all scalar matrices
    if (associated(rproblemLevel%rmatrix)) then
      do i = lbound(rproblemLevel%rmatrix,1), ubound(rproblemLevel%rmatrix,1)
        call lsyssc_releaseMatrix(rproblemLevel%rmatrix(i))
      end do
      deallocate(rproblemLevel%Rmatrix)
    end if
    
    ! Release all block matries
    if (associated(rproblemLevel%rmatrixBlock)) then
      do i = lbound(rproblemLevel%rmatrixBlock,1),&
             ubound(rproblemLevel%rmatrixBlock,1)
        call lsysbl_releaseMatrix(rproblemLevel%rmatrixBlock(i))
      end do
      deallocate(rproblemLevel%RmatrixBlock)
    end if
    
    ! Release all scalar vectors
    if (associated(rproblemLevel%rvector)) then
      do i = lbound(rproblemLevel%rvector,1),&
             ubound(rproblemLevel%rvector,1)
        call lsyssc_releaseVector(rproblemLevel%rvector(i))
      end do
      deallocate(rproblemLevel%Rvector)
    end if
    
    ! Release all block vectors
    if (associated(rproblemLevel%rvectorBlock)) then
      do i = lbound(rproblemLevel%rvectorBlock,1),&
             ubound(rproblemLevel%rvectorBlock,1)
        call lsysbl_releaseVector(rproblemLevel%rvectorBlock(i))
      end do
      deallocate(rproblemLevel%RvectorBlock)
    end if
    
    ! Release stabilization structure
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

  subroutine problem_appendLevel(rproblem, rproblemLevel, rproblemLevelRef)

!<description>
    ! This subroutine appends a multigrid level structure to the
    ! linked list of multigrid level structures. If the optional
    ! reference level rproblemLevelRef is given, the multigrid level structure
    ! is appended to rproblemLevelRef. Otherwise, the maximum multigrid level
    ! is used as reference structure.
!</description>
  
!<inputoutput>
    ! global problem structure
    type(t_problem), intent(INOUT), target :: rproblem

    ! multigrid level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! OPTIONAL: reference multigrid structure
    type(t_problemLevel), intent(INOUT), target, optional :: rproblemLevelRef
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

        ! Set pointer to maximum multigrid level
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

  subroutine problem_prependLevel(rproblem, rproblemLevel, rproblemLevelRef)

!<description>
    ! This subroutine prepends a multigrid level structure to the
    ! linked list of multigrid level structures. If the optional
    ! reference level rproblemLevelRef is given, the multigrid level structure
    ! is prepended to rproblemLevelRef. Otherwise, the minimum multigrid level
    ! is used as reference structure.
!</description>
  
!<inputoutput>
    ! global problem structure
    type(t_problem), intent(INOUT), target :: rproblem

    ! multigrid level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! OPTIONAL: reference multigrid structure
    type(t_problemLevel), intent(INOUT), target, optional :: rproblemLevelRef
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

        ! Set pointer to minimum multigrid level
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
    ! This subroutine removes an existing multigrid level structure
    ! from an existing problem structure
!</description>

!<inputoutput>
    ! global problem structure
    type(t_problem), intent(INOUT), target :: rproblem

    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
!</inputoutput>
!</subroutine>
    
    if (.not.associated(rproblemLevel%p_rproblem, rproblem)) then
      call output_line('Multigrid level structure does not belong to problem structure!',&
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
    ! This subroutine outputs information about the multigrid level structure
!</description>

!<input>
    ! multigrid level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel
!</input>
!</subroutine>

    call output_line ('ProblemLevel: '//trim(sys_siL(rproblemLevel%ilev,15)))
    call output_line ('-------------')

    if (associated(rproblemLevel%p_rproblemLevelCoarse)) then
      call output_line ('coarser level: '//&
                        trim(sys_siL(rproblemLevel%p_rproblemLevelCoarse%ilev,15)))
    else
      call output_line ('coarser level: not associated')
    end if

    if (associated(rproblemLevel%p_rproblemLevelFine)) then
      call output_line ('finer level:   '//&
                        trim(sys_siL(rproblemLevel%p_rproblemLevelFine%ilev,15)))
    else
      call output_line ('finer level:   not associated')
    end if

    if (associated(rproblemLevel%Rafcstab)) then
      call output_line ('Rafcstab:      '//&
                        trim(sys_siL(size(rproblemLevel%Rafcstab),15)))
    else
      call output_line ('Rafcstab:      not associated')
    end if

    if (associated(rproblemLevel%Rmatrix)) then
      call output_line ('Rmatrix:       '//&
                        trim(sys_siL(size(rproblemLevel%Rmatrix),15)))
    else
      call output_line ('Rmatrix:       not associated')
    end if

    if (associated(rproblemLevel%RmatrixBlock)) then
      call output_line ('RmatrixBlock:  '//&
                        trim(sys_siL(size(rproblemLevel%RmatrixBlock),15)))
    else
      call output_line ('RmatrixBlock:  not associated')
    end if

    if (associated(rproblemLevel%Rvector)) then
      call output_line ('Rvector:       '//&
                        trim(sys_siL(size(rproblemLevel%Rvector),15)))
    else
      call output_line ('Rvector:       not associated')
    end if

    if (associated(rproblemLevel%RvectorBlock)) then
      call output_line ('RvectorBlock:  '//&
                        trim(sys_siL(size(rproblemLevel%RvectorBlock),15)))
    else
      call output_line ('RvectorBlock:  not associated')
    end if
    call output_lbrk()

  end subroutine problem_infoLevel

  !*****************************************************************************
  
!<subroutine>

  subroutine problem_createProfile(rproblem, sfilename, ssectionname,&
                                   cvariables, rparser, istatus, time, rvector)

!<description>
    ! This subroutine initializes a vector profile from the parameters
    ! specified in the parameter file by calling a function parser w.r.t.
    ! the specified variable names
!</description>
    
!<input>
    ! global problem structure
    type(t_problem), intent(IN) :: rproblem

    ! name of parameter file
    character(LEN=*), intent(IN) :: sfilename

    ! name of the parameter section
    character(LEN=*), intent(IN) :: ssectionname

    ! symbolic variable names
    character(LEN=*), dimension(:), intent(IN) :: cvariables
    
    ! OPTIONAL: simulation time
    ! if this parameter is not specified, then time=0 is assumed
    real(DP), intent(IN), optional :: time
!</input>

!<inputoutput>
    ! OPTIONAL: profile vector
    ! If this vector is not present then no explicit profile is
    ! generated and only the parser is filled with data
    type(t_vectorBlock), intent(INOUT), optional :: rvector
!</inputoutput>

!<output>
    ! Parser
    type(t_fparser), intent(OUT) :: rparser

    ! Error flag: < 0 if profile was not created
    integer, intent(OUT) :: istatus
!</output>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(size(cvariables)) :: DvariableValues
    real(DP) :: ttime
    character(SYS_STRLEN) :: keyword
    character(LEN=1024) :: sdata
    character(LEN=1024) :: expression
    integer :: ieq,neq,iunit,ncomp,icomp,ipos,ios,idatalen
    
    !----------------------------------------------------------------
    ! Initialize parser
    !----------------------------------------------------------------

    ! Try to open the file
    call io_openFileForReading(sfilename, iunit, .true.)

    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'problem_createProfile')
      istatus = -1
      return
    end if

    ! Read through the input file until the given keyword is reached.
    ! If the key word is not present, then set ISTATUS=-1 and return.
    ios = 0
    do while(ios .eq. 0)
      
      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)
      if (ios .ne. 0) then
        call output_line('Unable to read KEYWORD from input file!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'problem_createProfile')
        istatus = -1
        close(iunit)
        return
      end if

      ! Check for keyword
      call sys_tolower(sdata(1:idatalen), keyword)
      if (trim(adjustl(keyword)) .eq. trim(adjustl(ssectionname))) exit
    end do
    
    ! We found the keyword. String NCOMP must be the next line to read.
    call io_readlinefromfile(iunit, sdata, idatalen, ios)
    if (ios .ne. 0) then
      call output_line('Unable to read data from input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'problem_createProfile')
      istatus = -1
      close(iunit)
      return
    end if

    ! Check for keyword NCOMP
    call sys_tolower(sdata(1:idatalen), keyword)
    if (trim(adjustl(keyword)) .ne. 'ncomp') then
      call output_line('Syntax error in input file!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'problem_createProfile')
      call sys_halt()
    end if
    
    ! Read value for NCOMP
    read(iunit,*,IOSTAT=ios) ncomp
    if (ios .ne. 0) then
      call output_line('Unable to read value of NCOMP from input file!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'problem_createProfile')
      call sys_halt()
    end if
    

    ! Create parser structure
    call fparser_create(rparser, ncomp)
    
    ! Read expressions into parser, relations and concatinations
    do icomp = 1, ncomp
      
      ! Read until expression is finished
      ios  = 0
      ipos = 1
      expression = " "

      do while(ios .eq. 0)

        ! Read next line in file
        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'problem_createProfile')
          call sys_halt()
        end if

        ! Append line to expression
        expression(ipos:) = sdata(1:idatalen)
        ipos = len_trim(expression)
        
        ! Check if expression is continued in the following line
        if (expression(max(1,ipos-2):ipos) .eq. '...') then
          ipos = ipos-2
        else
          exit
        end if
      end do

      call fparser_parseFunction(rparser, icomp, expression(1:ipos), cvariables)
    end do

    ! We are done. Close the unit.
    close(iunit)
    istatus = 0
    
    !----------------------------------------------------------------
    ! Create profile vector (if required)
    !----------------------------------------------------------------
    
    if (.not.present(rvector)) return
    
    ! Set pointers
    call storage_getbase_double2D(&
        rproblem%p_rproblemLevelMax%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    
    if (present(time)) then
      ttime = time
    else
      ttime = 0.0_DP
    end if
    
    if (ncomp .eq. rvector%nblocks) then
      ! Each component is stored separately in one block.
      ! Thus, loop over all components and evaluate the profile blockwise.
      do icomp = 1, ncomp

        ! Get number of equations for current block
        neq = rvector%RvectorBlock(icomp)%NEQ
        call lsyssc_getbase_double(rvector%RvectorBlock(icomp), p_Dx)
        
        ! Loop over all equations in the block
        do ieq = 1, neq
          DvariableValues = (/p_DvertexCoords(:,ieq), ttime/)
          call fparser_evalFunction(rparser, icomp, DvariableValues, p_Dx(ieq))
        end do
      end do
      
    elseif(rvector%nblocks .eq. 1 .and. &
           ncomp .eq. rvector%RvectorBlock(1)%NVAR) then
      ! The entire vector is stored in interleave format.

      ! Get number of equations for first vector block
      neq = rvector%RvectorBlock(1)%NEQ
      call lsyssc_getbase_double(rvector%RvectorBlock(1), p_Dx)

      ! Loop over all equations
      do ieq = 1, neq
        DvariableValues = (/p_DvertexCoords(:,ieq), ttime/)
        
        ! Loop over components
        do icomp = 1, ncomp
          call fparser_evalFunction(rparser, icomp, DvariableValues, p_Dx((ieq-1)*ncomp+icomp))
        end do
      end do

    else

      call output_line('Wrong number of vector components!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'problem_createProfile')
      call sys_halt()

    end if
  end subroutine problem_createProfile
end module problem

!##############################################################################
!# ****************************************************************************
!# <name> hadaptivity </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module contains all routines which are required to perform
!# h-adaptivity, namely, grid refinement and grid coarsening.  In
!# order to apply h-adaptivity, the initial mesh must be conforming,
!# that is, it is allowed to contain a mixture of triangular and
!# quadrilateral elements without 'hanging' nodes. After one grid
!# adaptivity step, the resulting mesh is also conforming without
!# 'hanging' nodes.
!#
!# This module makes extensive use of dynamic data structures such as
!# quadtrees, octrees, maps and so-called arraylists which need to be
!# generated from the static mesh structure. After grid adaptivity,
!# the dynamic data needs to be reconverted to the static structures
!# required for the simulation.
!#
!# Note that this module is implemented as a 'black-box' tool for grid
!# adaptivity.  One of the building blocks is the t_hadapt data
!# structure which provides all required information.
!#
!#
!# The following routines are available:
!#
!#  1.) hadapt_initFromParameterlist
!#      -> Initialises adaptivity structure from parameterlist
!#
!#  2.) hadapt_initFromTriangulation
!#      -> Initialises adaptivity structure from triangulation structure
!#
!#  3.) hadapt_generateRawMesh
!#      ->  Generates the raw mesh from the adaptivity structure
!#
!#  4.) hadapt_releaseAdaptation
!#      -> Releases all internal adaptation structures
!#
!#  5.) hadapt_duplicateAdaptation
!#      -> Creates a duplicate / backup of an adaptivity structure.
!#
!#  6.) hadapt_restoreAdaptation
!#      -> Restores an adaptivity structure previously backed up with
!#         hadapt_duplicateAdapation
!#
!#  7.) hadapt_refreshAdaptation
!#      -> Refreshes pointers of adaptation structure
!#
!#  8.) hadapt_performAdaptation
!#      -> Performes one step of grid adaptation
!#
!#  9.) hadapt_infoStatistics
!#      -> Outputs information about the adaptivity structure
!#
!# 10.) hadapt_writeGridGMV
!#      -> Writes the adapted grid to file in GMV format
!#
!# 12.) hadapt_checkConsistency
!#      -> Checks the internal consistency of dynamic data structures
!#
!# </purpose>
!##############################################################################

module hadaptivity

!$use omp_lib
  use arraylistInt
  use basicgeometry
  use boundary
  use collection
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptaux1d
  use hadaptaux2d
!!$  use hadaptaux3d
  use io
  use linearsystemscalar
  use mapInt_DP
  use octreeDP
  use paramlist
  use quadtreeDP
  use sort
  use storage
  use triangulation

  implicit none

  private
  public :: t_hadapt
  public :: hadapt_initFromParameterlist
  public :: hadapt_initFromTriangulation
  public :: hadapt_generateRawMesh
  public :: hadapt_releaseAdaptation
  public :: hadapt_duplicateAdaptation
  public :: hadapt_restoreAdaptation
  public :: hadapt_refreshAdaptation
  public :: hadapt_performAdaptation
  public :: hadapt_infoStatistics
  public :: hadapt_writeGridGMV
  public :: hadapt_checkConsistency

contains

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_initFromParameterlist(rhadapt, rparlist, ssectionName)

!<description>
    ! This subroutine initialises the adaptivity structure
    ! with the values supplied by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! name of the section
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! adaptivity structure
    type(t_hadapt), intent(out) :: rhadapt
!</output>
!</subroutine>

    ! Get mandatory parameters from list
    call parlst_getvalue_int   (rparlist, ssectionName, "nsubdividemax",&
                                rhadapt%nsubdividemax)
    call parlst_getvalue_int   (rparlist, ssectionName, "iadaptationStrategy",&
                                rhadapt%iadaptationStrategy)
    call parlst_getvalue_double(rparlist, ssectionName, "drefinementTolerance",&
                                rhadapt%drefinementTolerance)
    call parlst_getvalue_double(rparlist, ssectionName , "dcoarseningTolerance",&
                                rhadapt%dcoarseningTolerance)

    ! Initialise data
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_PARAMETERS)

  end subroutine hadapt_initFromParameterlist

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_initFromTriangulation(rhadapt, rtriangulation)

!<description>
    ! This subroutine initialises all required components of the adaptativit
    ! structure from the triangulation structure rtriangulation.
!</description>

!<input>
    ! triangulation structure
    type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Initialise duplication flag
    rhadapt%iduplicationFlag = 0

    ! Set coordinates
    select case(rtriangulation%ndim)
    case (NDIM1D)
      call hadapt_setVertexCoords1D(rhadapt,&
                                    rtriangulation%h_DvertexCoords,&
                                    rtriangulation%NVT)
    case(NDIM2D)
      call hadapt_setVertexCoords2D(rhadapt,&
                                    rtriangulation%h_DvertexCoords,&
                                    rtriangulation%NVT)
    case(NDIM3D)
      call hadapt_setVertexCoords3D(rhadapt,&
                                    rtriangulation%h_DvertexCoords,&
                                    rtriangulation%NVT)
    case default
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_initFromTriangulation')
      call sys_halt()
    end select

    ! Set nodal properties
    call hadapt_setNodalProperty(rhadapt,&
                                 rtriangulation%h_InodalProperty)

    ! Set element numbers
    call hadapt_setNelOfType(rhadapt,&
                             rtriangulation%InelOfType)

    ! Set vertices at element
    call hadapt_setVerticesAtElement(rhadapt,&
                                     rtriangulation%h_IverticesAtElement,&
                                     rtriangulation%NEL)

    ! Set elements adjacent to element
    call hadapt_setNeighboursAtElement(rhadapt,&
                                       rtriangulation%h_IneighboursAtElement)

    ! Set boundary parametrisation for 2D
    if ((rtriangulation%ndim .eq. NDIM2D) .and.&
        (rtriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE)) then
      call hadapt_setBoundary(rhadapt, rtriangulation%h_IboundaryCpIdx,&
                              rtriangulation%h_IverticesAtBoundary,&
                              rtriangulation%h_DvertexParameterValue,&
                              rtriangulation%NBCT, rtriangulation%NVBD)
    end if

    ! Generate "elements-meeting-at-vertex" structure
    call hadapt_genElementsAtVertex(rhadapt)

    ! Create generation array and initialise all nodes with "age" 0
    call storage_new('hadapt_initFromTriangulation','p_IvertexAge',&
                     rhadapt%NVT, ST_INT, rhadapt%h_IvertexAge, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rhadapt%h_IvertexAge, rhadapt%p_IvertexAge)

  end subroutine hadapt_initFromTriangulation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_generateRawMesh(rhadapt, rtriangulation)

!<description>
    ! This subroutine generates an (extended)
    ! raw mesh from the adaptivity structure.
!</description>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! triangulation structure
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
!</subroutine>

    ! Get dimensions
    rtriangulation%ndim = rhadapt%ndim

    ! Get coordinates
    select case(rhadapt%ndim)
    case (NDIM1D)
      call hadapt_getVertexCoords1D(rhadapt,&
                                    rtriangulation%h_DvertexCoords,&
                                    rtriangulation%NVT)
    case(NDIM2D)
      call hadapt_getVertexCoords2D(rhadapt,&
                                    rtriangulation%h_DvertexCoords,&
                                    rtriangulation%NVT)
    case(NDIM3D)
      call hadapt_getVertexCoords3D(rhadapt,&
                                    rtriangulation%h_DvertexCoords,&
                                    rtriangulation%NVT)
    case default
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_generateRawMesh')
      call sys_halt()
    end select

    ! Get number of elements
    call hadapt_getNelOfType(rhadapt,&
                             rtriangulation%InelOfType)

    ! Get vertices at element list
    call hadapt_getVerticesAtElement(rhadapt,&
                                     rtriangulation%h_IverticesAtElement,&
                                     rtriangulation%NEL)

    ! Get element neighbours
    call hadapt_getNeighboursAtElement(rhadapt,&
                                       rtriangulation%h_IneighboursAtElement)

    ! Get nodal property list
    call hadapt_getNodalProperty(rhadapt,&
                                 rtriangulation%h_InodalProperty)

    ! Get boundary for 2D
    if ((rtriangulation%ndim .eq. NDIM2D) .and. &
        (rtriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE)) then
      call hadapt_getBoundary(rhadapt, rtriangulation%h_IboundaryCpIdx,&
                              rtriangulation%h_IverticesAtBoundary,&
                              rtriangulation%h_DvertexParameterValue,&
                              rtriangulation%NVBD, rtriangulation%NBCT)
    end if

    ! Generate extended raw mesh information for edge (and face)
    ! numbering.
    call tria_initExtendedRawMesh (rtriangulation)

  end subroutine hadapt_generateRawMesh

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_releaseAdaptation(rhadapt)

!<description>
    ! This subroutine releases all internal structures of the
    ! adaptivity structure rhadapt.
!</description>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer(I32) :: idupflag
    integer :: ibct

    idupflag = rhadapt%iduplicationFlag

    ! Check if quadtree exists
    if (iand(rhadapt%iSpec,HADAPT_HAS_COORDS) .eq.&
                           HADAPT_HAS_COORDS) then
      select case(rhadapt%ndim)
      case (NDIM1D)
        call storage_free(rhadapt%h_DvertexCoords1D)
        nullify(rhadapt%p_DvertexCoords1D)

      case(NDIM2D)
        call qtree_release(rhadapt%rVertexCoordinates2D)

      case(NDIM3D)
        call otree_release(rhadapt%rVertexCoordinates3D)

      case default
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_releaseAdaptation')
        call sys_halt()
      end select
    end if

    ! Check if boundary structure exists
    if (associated(rhadapt%rBoundary)) then
      do ibct = 1, size(rhadapt%rBoundary, 1)
        call map_release(rhadapt%rBoundary(ibct))
      end do
      deallocate(rhadapt%rBoundary)
      nullify(rhadapt%rBoundary)
    end if

    ! Release elements-meeting-at-vertex arraylist
    call alst_release(rhadapt%relementsAtVertex)

    ! Release storage which is no longer in use
    call checkAndRelease(idupflag, HADAPT_SHARE_IMARKER,&
                         rhadapt%h_Imarker)
    call checkAndRelease(idupflag, HADAPT_SHARE_IVERTEXAGE,&
                         rhadapt%h_IvertexAge)
    call checkAndRelease(idupflag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
                         rhadapt%h_ImidneighboursAtElement)

    ! Nullify "performance-pointers"
    nullify(rhadapt%p_IvertexAge)
    nullify(rhadapt%p_InodalProperty)
    nullify(rhadapt%p_IverticesAtElement)
    nullify(rhadapt%p_IneighboursAtElement)
    nullify(rhadapt%p_ImidneighboursAtElement)

    ! Clear parameters
    rhadapt%iadaptationStrategy  = HADAPT_NOADAPTATION
    rhadapt%drefinementTolerance = 0._DP
    rhadapt%dcoarseningTolerance = 0._DP

    ! Clear data
    rhadapt%iSpec            = HADAPT_UNDEFINED
    rhadapt%nRefinementSteps = 0
    rhadapt%nCoarseningSteps = 0
    rhadapt%nSmoothingSteps  = 0
    rhadapt%nGreenElements   = 0
    rhadapt%ndim             = 0
    rhadapt%NVT              = 0
    rhadapt%NVT0             = 0
    rhadapt%increaseNVT      = 0
    rhadapt%NVBD             = 0
    rhadapt%NVBD0            = 0
    rhadapt%NBCT             = 0
    rhadapt%NEL              = 0
    rhadapt%NEL0             = 0
    rhadapt%NELMAX           = 0
    rhadapt%InelOfType       = 0
    rhadapt%InelOfType0      = 0

  contains

    subroutine checkAndRelease (idupFlag, ibitfield, ihandle)

      ! Checks if idupFlag has all bits ibitfield set.
      ! If yes, the handle ihandle is set removed
      ! Otherwise, the the memory behind ihandle is freed.

      integer(I32), intent(in) :: ibitfield
      integer(I32), intent(in) :: idupFlag
      integer, intent(inout) :: ihandle

      if (iand(idupFlag, ibitfield) .ne. ibitfield) then
        if (ihandle .ne. ST_NOHANDLE) call storage_free(ihandle)
      else
        ihandle = ST_NOHANDLE
      end if

    end subroutine checkAndRelease

  end subroutine hadapt_releaseAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_duplicateAdaptation(rhadapt, rhadaptBackup,&
                                        iduplicationFlag, bupdate)

!<description>
    ! This subroutine makes a copy of an adaptivity structure in memory. The
    ! variable iduplicationFlag decides on which arrays are copied in memory
    ! and which are not.
    !
    ! By setting the corresponding bit in iduplicationFlag to 0, the array is
    ! duplicated in memory, and any change to the new array will not harm
    ! the original one.
    ! By setting a flag HADAPT_SHARE_xxxx in iduplicationFlag, the corresponding
    ! array is not duplicated. The handle of the original structure is simply put
    ! into the new structure to make the information accessable. Then this array
    ! is shared between two adaptivity structures!
!</description>

!<input>
    ! The source structure which provides all information
    type(t_hadapt), intent(in) :: rhadapt

    ! Bitfield that decides which handles are a copy of another structure, thus
    ! which arrays are shared between the new and the old structure.
    ! The bitfield contains a combination of HADAPT_SHARE_xxxx constants. Every
    ! information whose flag is set in iduplicationFlag is shared between rhadapt
    ! and rhadaptBackup.
    ! Therefore e.g., iduplicationFlag=0 copies all arrays from rhadapt in memory,
    ! while HADAPT_SHARE_ALL will copy nothing, but will share everything between
    ! rhadapt and rhadaptBackup.
    integer(I32), intent(in) :: iduplicationFlag

    ! OPTIONAL. Defines how to create the backup.
    ! = .FALSE.: Treat rhadaptBackup as empty destination structure. If necessary,
    !    information in rhadaptBackup is released. rhadaptBackup is rebuild
    !    according to rhadapt and iduplicationFlag.
    !    This is the standard setting if bupdate is not specified.
    ! = .TRUE. : Treat rhadaptBackup as existing copy of rhadapt which has to be
    !    updated. Recover all arrays of rhadaptBackup by those of rhadapt.
    !    (I.e. those arrays which were duplicated by a previous
    !    call to iduplicationFlag with IUPD=0.)
    !    iduplicationFlag can used to specify which data do copy
    !    from rhadapt to rhadaptBackup. It is OR`ed with rhadaptBackup%iduplicationFlag
    !    to get the actual duplication flag. This leads to the following interpretation
    !    of iduplicationFlag:
    !     =0:  Copy all data that was copied previously. This is the usual setting.
    !    <>0:  Copy all arrays where the corresponding flag in iduplicationFlag is not
    !          set and which exist as duplicates. Arrays corresponding to flags which
    !          are set or where the handles in rhadapt and rhadaptBackup coincide
    !          are not touched.
    logical, intent(in), optional :: bupdate
!</input>

!<inputoutput>
    ! The destination structure which receives the information
    type(t_hadapt), intent(inout) :: rhadaptBackup
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: idupFlag
    integer :: ibct
    logical :: bupd

    bupd = .false.
    if (present(bupdate)) bupd = bupdate

    if (.not. bupd) then
      ! Release any old data.
      call hadapt_releaseAdaptation(rhadaptBackup)

      rhadaptBackup%iSpec                = rhadapt%iSpec
      rhadaptBackup%iadaptationStrategy  = rhadapt%iadaptationStrategy
      rhadaptBackup%NSUBDIVIDEMAX        = rhadapt%NSUBDIVIDEMAX
      rhadaptBackup%nRefinementSteps     = rhadapt%nRefinementSteps
      rhadaptBackup%nCoarseningSteps     = rhadapt%nCoarseningSteps
      rhadaptBackup%nSmoothingSteps      = rhadapt%nSmoothingSteps
      rhadaptBackup%drefinementTolerance = rhadapt%dRefinementTolerance
      rhadaptBackup%dCoarseningTolerance = rhadapt%dCoarseningTolerance
      rhadaptBackup%ndim                 = rhadapt%ndim
      rhadaptBackup%NVT0                 = rhadapt%NVT0
      rhadaptBackup%NVT                  = rhadapt%NVT
      rhadaptBackup%increaseNVT          = rhadapt%increaseNVT
      rhadaptBackup%NVBD0                = rhadapt%NVBD0
      rhadaptBackup%NVBD                 = rhadapt%NVBD
      rhadaptBackup%NBCT                 = rhadapt%NBCT
      rhadaptBackup%NEL0                 = rhadapt%NEL0
      rhadaptBackup%NEL                  = rhadapt%NEL
      rhadaptBackup%NELMAX               = rhadapt%NELMAX
      rhadaptBackup%nGreenElements       = rhadapt%nGreenElements
      rhadaptBackup%InelOfType0          = rhadapt%InelOfType0
      rhadaptBackup%InelOfType           = rhadapt%InelOfType

      ! Decide on iduplicationFlag which arrays to copy
      rhadaptBackup%iduplicationFlag = iduplicationFlag
      idupFlag = iduplicationFlag

    else

      ! Create a bitfiled what to copy by ORing iduplicationFlag with what
      ! we have in rhadaptBackup. That way, only arrays that exist as real
      ! duplicates are copied from rhadapt to rhadaptBackup.

      idupFlag = ior(iduplicationFlag, rhadaptBackup%iduplicationFlag)

    end if

    ! Call checkAndCopy for all components. This will either copy the handle
    ! or allocate new memory and copy the content of the component.

    ! Bit   0: Imarker
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMARKER,&
                      rhadapt%h_Imarker, rhadaptBackup%h_Imarker)

    ! Bit   1: IvertexAge
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTEXAGE,&
                      rhadapt%h_IvertexAge, rhadaptBackup%h_IvertexAge)
    call storage_getbase_int(rhadaptBackup%h_IvertexAge,&
                             rhadaptBackup%p_IvertexAge)

    ! Bit   2: InodalProperty
    call checkAndCopy(idupFlag, HADAPT_SHARE_INODALPROPERTY,&
                      rhadapt%h_InodalProperty, rhadaptBackup%h_InodalProperty)
    call storage_getbase_int(rhadaptBackup%h_InodalProperty,&
                             rhadaptBackup%p_InodalProperty)

    ! Bit   3: IverticesAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTICESATELEMENT,&
                      rhadapt%h_IverticesAtElement,&
                      rhadaptBackup%h_IverticesAtElement)
    call storage_getbase_int2D(rhadaptBackup%h_IverticesAtElement,&
                               rhadaptBackup%p_IverticesAtElement)

    ! Bit   4: IneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_INEIGHATELEMENT,&
                      rhadapt%h_IneighboursAtElement,&
                      rhadaptBackup%h_IneighboursAtElement)
    call storage_getbase_int2D(rhadaptBackup%h_IneighboursAtElement,&
                               rhadaptBackup%p_IneighboursAtElement)

    ! Bit   5: ImidneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
                      rhadapt%h_ImidneighboursAtElement,&
                      rhadaptBackup%h_ImidneighboursAtElement)
    call storage_getbase_int2D(rhadaptBackup%h_ImidneighboursAtElement,&
                               rhadaptBackup%p_ImidneighboursAtElement)

    ! Bit   6: rVertexCoordinates
    if (iand(idupFlag, HADAPT_SHARE_RVERTEXCOORDINATES) .ne.&
                       HADAPT_SHARE_RVERTEXCOORDINATES) then

      select case(rhadaptBackup%ndim)
      case (NDIM1D)
        call storage_copy(rhadapt%h_DvertexCoords1D,&
                          rhadaptBackup%h_DvertexCoords1D)
        call storage_getbase_double2D(rhadaptBackup%h_DvertexCoords1D,&
                                      rhadaptBackup%p_DvertexCoords1D)
      case(NDIM2D)
        call qtree_duplicate(rhadapt%rVertexCoordinates2D,&
                             rhadaptBackup%rVertexCoordinates2D)
      case(NDIM3D)
        call otree_duplicate(rhadapt%rVertexCoordinates3D,&
                                   rhadaptBackup%rVertexCoordinates3D)
      case default
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_duplicateAdaptation')
        call sys_halt()
      end select
    end if

    ! Bit   7: rBoundary
    if (iand(idupFlag, HADAPT_SHARE_RBOUNDARY) .ne.&
                       HADAPT_SHARE_RBOUNDARY) then

      if (associated(rhadaptBackup%rBoundary)) then
        do ibct = 1, size(rhadaptBackup%rBoundary,1)
          call map_release(rhadaptBackup%rBoundary(ibct))
        end do
        deallocate(rhadaptBackup%rBoundary)
      end if

      allocate(rhadaptBackup%rBoundary(size(rhadapt%rBoundary,1)))
      do ibct = 1, size(rhadapt%rBoundary,1)
        call map_duplicate(rhadapt%rBoundary(ibct),&
                           rhadaptBackup%rBoundary(ibct))
      end do
    end if

    ! Bit   8: rElementsAtVertex
    if (iand(idupFlag, HADAPT_SHARE_RELEMENTSATVERTEX) .ne.&
                       HADAPT_SHARE_RELEMENTSATVERTEX) then
      call alst_duplicate(rhadapt%rElementsAtVertex,&
                          rhadaptBackup%rElementsAtVertex)
    end if

  contains

    subroutine checkAndCopy (idupFlag, ibitfield, isourcehandle, idesthandle)

      ! Checks if idupFlag has all bits ibitfield set.
      ! If yes, idesthandle is set to isourcehandle.
      ! Otherwise, the memory behind isourcehandle is duplicated in memory
      ! and idesthandle receives the handle to the new memory block.

      integer(I32), intent(in) :: ibitfield
      integer(I32), intent(in) :: idupFlag
      integer, intent(in) :: isourcehandle
      integer, intent(inout) :: idesthandle

      if (iand(idupFlag, ibitfield) .ne. ibitfield) then
        if (isourcehandle .ne. ST_NOHANDLE) then
          call storage_copy(isourcehandle, idesthandle)
        end if
      else
        idesthandle = isourcehandle
      end if

    end subroutine checkAndCopy

  end subroutine hadapt_duplicateAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_restoreAdaptation(rhadaptBackup, rhadapt)

!<description>
    ! This subroutine restares data of an adaptivity structure. All information
    ! arrays which are not shared between rhadaptBackup and another adaptivity
    ! structure is copied (back) to rhadapt.
!</description>

!<input>
    ! Backup of an adaptivity structure
    type(t_hadapt), intent(in) :: rhadaptBackup
!</input>

!<inputoutput>
    ! Destination adaptivity structure
    ! All components where a duplicates exist in rhadaptBackup are copied
    ! to rhadapt, overwriting the old information arrays.
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: idupFlag
    integer :: ibct

    idupFlag = rhadapt%iduplicationFlag

    rhadapt%iSpec                = rhadaptBackup%iSpec
    rhadapt%iadaptationStrategy  = rhadaptBackup%iadaptationStrategy
    rhadapt%NSUBDIVIDEMAX        = rhadaptBackup%NSUBDIVIDEMAX
    rhadapt%nRefinementSteps     = rhadaptBackup%nRefinementSteps
    rhadapt%nCoarseningSteps     = rhadaptBackup%nCoarseningSteps
    rhadapt%nSmoothingSteps      = rhadaptBackup%nSmoothingSteps
    rhadapt%drefinementTolerance = rhadaptBackup%dRefinementTolerance
    rhadapt%dCoarseningTolerance = rhadaptBackup%dCoarseningTolerance
    rhadapt%ndim                 = rhadaptBackup%ndim
    rhadapt%NVT0                 = rhadaptBackup%NVT0
    rhadapt%NVT                  = rhadaptBackup%NVT
    rhadapt%increaseNVT          = rhadaptBackup%increaseNVT
    rhadapt%NVBD0                = rhadaptBackup%NVBD0
    rhadapt%NVBD                 = rhadaptBackup%NVBD
    rhadapt%NBCT                 = rhadaptBackup%NBCT
    rhadapt%NEL0                 = rhadaptBackup%NEL0
    rhadapt%NEL                  = rhadaptBackup%NEL
    rhadapt%NELMAX               = rhadaptBackup%NELMAX
    rhadapt%nGreenElements       = rhadaptBackup%nGreenElements
    rhadapt%InelOfType0          = rhadaptBackup%InelOfType0
    rhadapt%InelOfType           = rhadaptBackup%InelOfType

    ! Call checkAndCopy for all components. This will either copy the handle
    ! or allocate new memory and copy the content of the component.

    ! Bit   0: Imarker
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMARKER,&
                      rhadapt%h_Imarker, rhadaptBackup%h_Imarker)

    ! Bit   1: IvertexAge
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTEXAGE,&
                      rhadapt%h_IvertexAge, rhadaptBackup%h_IvertexAge)
    call storage_getbase_int(rhadapt%h_IvertexAge,&
                             rhadapt%p_IvertexAge)

    ! Bit   2: InodalProperty
    call checkAndCopy(idupFlag, HADAPT_SHARE_INODALPROPERTY,&
                      rhadapt%h_InodalProperty, rhadaptBackup%h_InodalProperty)
    call storage_getbase_int(rhadapt%h_InodalProperty,&
                             rhadapt%p_InodalProperty)

    ! Bit   3: IverticesAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTICESATELEMENT,&
                      rhadapt%h_IverticesAtElement,&
                      rhadaptBackup%h_IverticesAtElement)
    call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                               rhadapt%p_IverticesAtElement)

    ! Bit   4: IneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_INEIGHATELEMENT,&
                      rhadapt%h_IneighboursAtElement,&
                      rhadaptBackup%h_IneighboursAtElement)
    call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                               rhadapt%p_IneighboursAtElement)

    ! Bit   5: ImidneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
                      rhadapt%h_ImidneighboursAtElement,&
                      rhadaptBackup%h_ImidneighboursAtElement)
    call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                               rhadapt%p_ImidneighboursAtElement)

    ! Bit   6: rVertexCoordinates
    if (iand(idupFlag, HADAPT_SHARE_RVERTEXCOORDINATES) .ne.&
                       HADAPT_SHARE_RVERTEXCOORDINATES) then

      select case(rhadaptBackup%ndim)
      case (NDIM1D)
        call storage_copy(rhadaptBackup%h_DvertexCoords1D,&
                          rhadapt%h_DvertexCoords1D)
        call storage_getbase_double2D(rhadapt%h_DvertexCoords1D,&
                                      rhadapt%p_DvertexCoords1D)
      case(NDIM2D)
        call qtree_restore(rhadaptBackup%rVertexCoordinates2D,&
                           rhadapt%rVertexCoordinates2D)
      case(NDIM3D)
        call otree_restore(rhadaptBackup%rVertexCoordinates3D,&
                                 rhadapt%rVertexCoordinates3D)
      case default
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_restoreAdaptation')
        call sys_halt()
      end select
    end if

    ! Bit   7: rBoundary
    if (iand(idupFlag, HADAPT_SHARE_RBOUNDARY) .ne.&
                       HADAPT_SHARE_RBOUNDARY) then
      if (size(rhadaptBackup%rBoundary,1) .ne. &
          size(rhadapt%rBoundary,1)) then
        call output_line('Invalid number of boundary components!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_restoreAdaptation')
        call sys_halt()
      end if
      do ibct = 1, size(rhadapt%rBoundary,1)
        call map_restore(rhadaptBackup%rBoundary(ibct),&
                         rhadapt%rBoundary(ibct))
      end do
    end if


    ! Bit   8: rElementsAtVertex
    if (iand(idupFlag, HADAPT_SHARE_RELEMENTSATVERTEX) .ne.&
                       HADAPT_SHARE_RELEMENTSATVERTEX) then
      call alst_restore(rhadaptBackup%rElementsAtVertex,&
                        rhadapt%rElementsAtVertex)
    end if

  contains

    subroutine checkAndCopy (idupFlag, ibitfield, idesthandle, isourcehandle)

      ! Checks if idupFlag has all bits ibitfield set.
      ! If not, the memory behind isourcehandle is copied to idesthandle
      ! overwriting all previous information.

      integer(I32), intent(in) :: ibitfield
      integer(I32), intent(in) :: idupFlag
      integer, intent(in) :: isourcehandle
      integer, intent(inout) :: idesthandle

      if (iand(idupFlag, ibitfield) .ne. ibitfield) then
        if (isourcehandle .ne. ST_NOHANDLE) then
          call storage_copy(isourcehandle, idesthandle)
        end if
      end if

    end subroutine checkAndCopy

  end subroutine hadapt_restoreAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_refreshAdaptation(rhadapt, rtriangulation)

!<description>
    ! This subroutine refreshes the pointers and internal structure of the
    ! adaptation structure.
    ! NOTE: The triangulation structure must be compatible with the adaptivity
    ! structure. This CANNOT be checked by the routine and must be guaranteed
    ! by the used. The philosophy behind this adaptivity structure is to not
    ! modify the triangulation externally. If external modifications are performed
    ! then all adaptivity structures must be released and rebuild from scratch.
    ! This is quite time-consuming, and moreover, the deletion of vertices is
    ! not possible since all vertices will be considered as initial vertices.
!</description>

!<input>
    ! Triangulation structure
    type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>


    ! Get pointer to InodalProperty
    if (iand(rhadapt%iSpec, HADAPT_HAS_NODALPROP) .eq.&
                            HADAPT_HAS_NODALPROP) then
      if (rhadapt%h_InodalProperty .eq. &
          rtriangulation%h_InodalProperty) then
        call storage_getbase_int(rhadapt%h_InodalProperty,&
                                 rhadapt%p_InodalProperty)
      else
        call output_line('Inconsistent handle h_InodalProperty',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
        call sys_halt()
      end if
    end if

    ! Get pointer to IverticesAtElement
    if (iand(rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .eq.&
                            HADAPT_HAS_VERTATELEM) then
      if (rhadapt%h_IverticesAtElement .eq. &
          rtriangulation%h_IverticesAtElement) then
        call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                                   rhadapt%p_IverticesAtElement)
      else
        call output_line('Inconsistent handle h_IverticesAtElement',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
        call sys_halt()
      end if
    end if

    ! Get pointer to IneighboursAtElement
    if (iand(rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM) .eq.&
                            HADAPT_HAS_NEIGHATELEM) then
      if (rhadapt%h_IneighboursAtElement .eq. &
          rtriangulation%h_IneighboursAtElement) then
        call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                                   rhadapt%p_IneighboursAtElement)
      else
        call output_line('Inconsistent handle h_IneighboursAtElement',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
        call sys_halt()
      end if
    end if

    ! Get pointer to ImidneighboursAtElement
    if (iand(rhadapt%iSpec, HADAPT_HAS_MIDNEIGH) .eq.&
                            HADAPT_HAS_MIDNEIGH) then
      if (rhadapt%h_ImidneighboursAtElement .ne. ST_NOHANDLE) then
        call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                   rhadapt%p_ImidneighboursAtElement)
      else
        call output_line('Inconsistent handle h_ImidneighboursAtElement',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
        call sys_halt()
      end if
    end if

    ! Get pointer to IvertexAge
    if (rhadapt%h_IvertexAge .ne. ST_NOHANDLE) then
      call storage_getbase_int(rhadapt%h_IvertexAge,&
                               rhadapt%p_IvertexAge)
    else
      call output_line('Inconsistent handle rhadapt%h_IvertexAge',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      call sys_halt()
    end if

  end subroutine hadapt_refreshAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_performAdaptation(rhadapt, rindicator, rcollection,&
                                      fcb_hadaptCallback)

!<description>
    ! This subroutine performs the complete adaptation process.
    ! First, the internal data structures are generated and
    ! elements are marked for refinement and coarsening based on
    ! the indicator vector and the tolerances for refinement and
    ! coarsening, respecitvely.
    ! Next, the list of marked elements is updated so as to obtain
    ! a globally conforming triangulation.
    ! Finally, mesh refinement and coarsening is performed
!</description>

!<input>
    ! Indicator vector for refinement
    type(t_vectorScalar), intent(in) :: rindicator

    ! callback routines
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: nvt,nel

    ! Check if dynamic data structures are available
    select case(rhadapt%ndim)
    case (NDIM1D)
      if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA1D) .ne.&
                              HADAPT_HAS_DYNAMICDATA1D) then
        call output_line('Dynamic data structures in 1D are not generated!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        call sys_halt()
      end if

    case (NDIM2D)
      if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                              HADAPT_HAS_DYNAMICDATA2D) then
        call output_line('Dynamic data structures in 2D are not generated!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        call sys_halt()
      end if

    case (NDIM3D)
      if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA3D) .ne.&
                              HADAPT_HAS_DYNAMICDATA3D) then
        call output_line('Dynamic data structures in 3D are not generated!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        call sys_halt()
      end if

    case default
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
      call sys_halt()
    end select


    ! Initialise initial dimensions
    call storage_getsize(rhadapt%h_IverticesAtElement, Isize)
    rhadapt%NELMAX = Isize(2)
    call storage_getsize(rhadapt%h_IneighboursAtElement, Isize)
    rhadapt%NELMAX = min(rhadapt%NELMAX, Isize(2))

    rhadapt%InelOfType0 = rhadapt%InelOfType
    rhadapt%NVT0        = rhadapt%NVT
    rhadapt%NEL0        = rhadapt%NEL
    rhadapt%NVBD0       = rhadapt%NVBD
    rhadapt%increaseNVT = 0

    ! What kind of grid refinement should be performed
    select case(rhadapt%iadaptationStrategy)

    ! No grid refinement
    case (HADAPT_NOADAPTATION)

    ! Red-green grid refinement
    case (HADAPT_REDGREEN)

      ! Which spatial dimensions are we?
      select case(rhadapt%ndim)
      case (NDIM1D)
        ! Mark elements for refinement based on indicator function
        call hadapt_markRefinement1D(rhadapt, rindicator)

        ! Mark element for recoarsening based on indicator function
        call hadapt_markCoarsening1D(rhadapt, rindicator)

        ! Compute new dimensions
        nvt = rhadapt%NVT+rhadapt%increaseNVT
        nel = hadapt_CalcNumberOfElements1D(rhadapt)

        ! Adjust array DvertexCoords1D
        call storage_realloc('hadapt_performAdaptation', nvt,&
                             rhadapt%h_DvertexCoords1D, ST_NEWBLOCK_NOINIT, .true.)
        call storage_getbase_double2D(rhadapt%h_DvertexCoords1D, rhadapt%p_DvertexCoords1D)

      case (NDIM2D)
        ! Mark elements for refinement based on indicator function
        call hadapt_markRefinement2D(rhadapt, rindicator)

        ! Mark additional elements to restore conformity
        call hadapt_markRedgreenRefinement2D(rhadapt, rcollection,&
                                             fcb_hadaptCallback)

        ! Mark element for recoarsening based on indicator function
        call hadapt_markRedgreenCoarsening2D(rhadapt, rindicator)

        ! Compute new dimensions
        nvt = rhadapt%NVT+rhadapt%increaseNVT
        nel = hadapt_CalcNumberOfElements2D(rhadapt)

      case (NDIM3D)
        call output_line('h-adaptivity not yet implemented in 3D!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        call sys_halt()
      end select


      ! Adjust array IvertexAge
      call storage_realloc('hadapt_performAdaptation', nvt,&
                           rhadapt%h_IvertexAge, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int(rhadapt%h_IvertexAge, rhadapt%p_IvertexAge)

      ! Adjust array InodalProperty
      if (iand(rhadapt%iSpec, HADAPT_HAS_NODALPROP) .eq.&
                              HADAPT_HAS_NODALPROP) then
        call storage_realloc('hadapt_performAdaptation', nvt,&
                             rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)
        call storage_getbase_int(rhadapt%h_InodalProperty, rhadapt%p_InodalProperty)
      end if

      ! Adjust array IverticesAtElement
      if (iand(rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .eq.&
                              HADAPT_HAS_VERTATELEM) then
        call storage_realloc('hadapt_performAdaptation', nel,&
                             rhadapt%h_IverticesAtElement, ST_NEWBLOCK_NOINIT, .true.)
        call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                                   rhadapt%p_IverticesAtElement)
      end if

      ! Adjust array IneighboursAtElement
      if (iand(rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM) .eq.&
                              HADAPT_HAS_NEIGHATELEM) then
        call storage_realloc('hadapt_performAdaptation', nel,&
                             rhadapt%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
        call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                                   rhadapt%p_IneighboursAtElement)
      end if

      ! Adjust array ImidneighboursAtElement
      if (iand(rhadapt%iSpec, HADAPT_HAS_MIDNEIGH) .eq.&
                              HADAPT_HAS_MIDNEIGH) then
        call storage_realloc('hadapt_performAdaptation', nel,&
                             rhadapt%h_ImidneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
        call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                   rhadapt%p_ImidneighboursAtElement)
      end if

      ! Adjust dimension of solution vector
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        rcollection%IquickAccess(1) = nvt
        call fcb_hadaptCallback(HADAPT_OPR_ADJUSTVERTEXDIM, rcollection)
      end if


      ! Which spatial dimensions are we?
      select case(rhadapt%ndim)
      case (NDIM1D)
        ! Perform element refinement in 1D
        call hadapt_refine1D(rhadapt, rcollection, fcb_hadaptCallback)

        ! Perform element coarsening in 1D
        call hadapt_coarsen1D(rhadapt, rcollection, fcb_hadaptCallback)

        ! Adjust nodal property array
        call storage_realloc('hadapt_performAdaptation', rhadapt%NVT,&
                             rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)

      case (NDIM2D)
        ! Perform element refinement in 2D
        call hadapt_refine2D(rhadapt, rcollection, fcb_hadaptCallback)

        ! Perform element coarsening in 2D
        call hadapt_coarsen2D(rhadapt, rcollection, fcb_hadaptCallback)

        ! Adjust nodal property array
        call storage_realloc('hadapt_performAdaptation', rhadapt%NVT,&
                             rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)

      case default
        call output_line('Unsupported spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported refinement strategy!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
      call sys_halt()
    end select

  end subroutine hadapt_performAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_infoStatistics(rhadapt)

!<description>
    ! This subroutine outputs statistical info about the adaptivity structure
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(in) :: rhadapt
!</input>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: ibct

    call output_line('Adaptivity statistics:')
    call output_line('----------------------')
    call output_line('Total number of grid refinement steps:           '//&
        trim(sys_siL(rhadapt%nRefinementSteps,3)))
    call output_line('Total number of grid coarsening steps:           '//&
        trim(sys_siL(rhadapt%nCoarseningSteps,3)))
    call output_line('Total number of grid smoothing  steps:           '//&
        trim(sys_siL(rhadapt%nSmoothingSteps,3)))
    call output_line('Strategy for grid refinement and coarsening:     '//&
        trim(sys_siL(rhadapt%iadaptationStrategy,3)))
    call output_line('Total number of elements (initially):            '//&
        trim(sys_siL(rhadapt%NEL,15))//"("//trim(sys_siL(rhadapt%NEL0,15))//")")
    call output_line('Total number of vertices (initially):            '//&
        trim(sys_siL(rhadapt%NVT,15))//"("//trim(sys_siL(rhadapt%NVT0,15))//")")
    call output_line('Total number of vertices at boundary(initially): '//&
        trim(sys_siL(rhadapt%NVBD,15))//"("//trim(sys_siL(rhadapt%NVBD0,15))//")")
    call output_lbrk

    call output_line('Handles:')
    call output_line('--------')
    call output_line('h_Imarker:                 '//&
        trim(sys_siL(rhadapt%h_Imarker,15)))
    call output_line('h_IvertexAge:              '//&
        trim(sys_siL(rhadapt%h_IvertexAge,15)))
    call output_line('h_InodalProperty:          '//&
        trim(sys_siL(rhadapt%h_InodalProperty,15)))
    call output_line('h_IverticesAtElement:      '//&
        trim(sys_siL(rhadapt%h_IverticesAtElement,15)))
    call output_line('h_IneighboursAtElement:    '//&
        trim(sys_siL(rhadapt%h_IneighboursAtElement,15)))
    call output_line('h_ImidneighboursAtElement: '//&
        trim(sys_siL(rhadapt%h_ImidneighboursAtElement,15)))
    call output_lbrk

    call output_line('Coordinates:')
    call output_line('------------')
    select case(rhadapt%ndim)
    case (NDIM1D)
      call storage_getsize(rhadapt%h_DvertexCoords1D, Isize)
      call output_line('DvertexCoords1D')
      call output_line('---------------')
      call output_line('NVT:               '//trim(sys_siL(rhadapt%NVT,15)))
      call output_line('NNVT:              '//trim(sys_siL(Isize(2),15)))
      call output_line('h_DvertexCoords1D: '//trim(sys_siL(rhadapt%h_DvertexCoords1D,15)))
    case(NDIM2D)
      call qtree_info(rhadapt%rVertexCoordinates2D)
    case(NDIM3D)
      call otree_info(rhadapt%rVertexCoordinates3D)
    end select
    call output_lbrk

    call output_line('Boundary:')
    call output_line('---------')
    do ibct=1,rhadapt%NBCT
      call output_line('Boundary component '//trim(sys_siL(ibct,3))//":")
      call output_line('------------------------')
      call map_info(rhadapt%rBoundary(ibct))
    end do
    call output_lbrk

    call output_line('Element at vertex structure:')
    call output_line('----------------------------')
    call alst_info(rhadapt%relementsAtVertex)

  end subroutine hadapt_infoStatistics

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_writeGridGMV(rhadapt, coutputFile)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in GMV format.
!</description>

!<input>
    ! output file name w/o suffix .gmv
    character(LEN=*), intent(in) :: coutputFile
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local parameters
    integer,  dimension(:), pointer :: p_Imarker
    integer :: ivt,iel,iunit,nve
    integer, save :: iout=0

    ! Check if dynamic data structures generated
    select case(rhadapt%ndim)
    case (NDIM1D)
      if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA1D) .ne.&
                              HADAPT_HAS_DYNAMICDATA1D) then
        call output_line('Dynamic data structures in 1D are not generated!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
        call sys_halt()
      end if

    case (NDIM2D)
      if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA2D) .ne.&
                              HADAPT_HAS_DYNAMICDATA2D) then
        call output_line('Dynamic data structures in 2D are not generated!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
        call sys_halt()
      end if

    case (NDIM3D)
      if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA3D) .ne.&
                              HADAPT_HAS_DYNAMICDATA3D) then
        call output_line('Dynamic data structures in 3D are not generated!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
        call sys_halt()
      end if

    case default
      call output_line('Unsupported spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
      call sys_halt()
    end select

    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Increment the sample number
    iout=iout+1

    ! Open output file for writing
    call io_openFileForWriting(trim(adjustl(coutputFile))//'.'//&
        trim(sys_siL(iout,5))//'.gmv', iunit, SYS_REPLACE, bformatted=.true.)
    write(UNIT=iunit,FMT='(A)') 'gmvinput ascii'

    ! Write vertices to output file
    select case(rhadapt%ndim)
    case (NDIM1D)
      write(UNIT=iunit,FMT=*) 'nodes ', rhadapt%NVT
      do ivt = 1, rhadapt%NVT
        write(UNIT=iunit,FMT=10) rhadapt%p_DvertexCoords1D(1,ivt)
      end do
      do ivt = 1, qtree_getSize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) 0._DP
      end do
      do ivt = 1, qtree_getSize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) 0._DP
      end do

    case (NDIM2D)
      write(UNIT=iunit,FMT=*) 'nodes ', rhadapt%NVT
      do ivt = 1, qtree_getSize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) qtree_getX(rhadapt%rVertexCoordinates2D, ivt)
      end do
      do ivt = 1, qtree_getSize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) qtree_getY(rhadapt%rVertexCoordinates2D, ivt)
      end do
      do ivt = 1, qtree_getSize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) 0._DP
      end do

    case (NDIM3D)
      write(UNIT=iunit,FMT=*) 'nodes ', rhadapt%NVT
      do ivt = 1, otree_getSize(rhadapt%rVertexCoordinates3D)
        write(UNIT=iunit,FMT=10) otree_getX(rhadapt%rVertexCoordinates3D, ivt)
      end do
      do ivt = 1, otree_getSize(rhadapt%rVertexCoordinates3D)
        write(UNIT=iunit,FMT=10) otree_getY(rhadapt%rVertexCoordinates3D, ivt)
      end do
      do ivt = 1, otree_getSize(rhadapt%rVertexCoordinates3D)
        write(UNIT=iunit,FMT=10) otree_getZ(rhadapt%rVertexCoordinates3D, ivt)
      end do

    case default
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
      call sys_halt()
    end select


    ! Write cells to output file
    write(UNIT=iunit,FMT=*) 'cells ', rhadapt%NEL
    do iel = 1, rhadapt%NEL
      nve = hadapt_getNVE(rhadapt,iel)
      select case(nve)
      case(TRIA_NVELINE1D)
        write(UNIT=iunit,FMT=*) 'line 2'
        write(UNIT=iunit,FMT=15) rhadapt%p_IverticesAtElement(1:TRIA_NVELINE1D, iel)

      case(TRIA_NVETRI2D)
        write(UNIT=iunit,FMT=*) 'tri 3'
        write(UNIT=iunit,FMT=20) rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D, iel)

      case(TRIA_NVEQUAD2D)
        write(UNIT=iunit,FMT=*) 'quad 4'
        write(UNIT=iunit,FMT=30) rhadapt%p_IverticesAtElement(1:TRIA_NVEQUAD2D, iel)

      case default
        call output_line('Invalid element type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
        call sys_halt()
      end select
    end do

    ! Write velocity to output file
    write(UNIT=iunit,FMT=*) 'velocity 1'
    do ivt = 1, rhadapt%NVT
      write(UNIT=iunit,FMT=10) 0._DP
      write(UNIT=iunit,FMT=10) 0._DP
      write(UNIT=iunit,FMT=10) 0._DP
    end do

    ! Write variable for vertex age
    write(UNIT=iunit,FMT=*) 'variable'
    write(UNIT=iunit,FMT=*) 'vert_age 1'

    do ivt = 1, min(rhadapt%NVT, size(rhadapt%p_IvertexAge, 1))
      write(UNIT=iunit,FMT=40) rhadapt%p_IvertexAge(ivt)
    end do
    do ivt = min(rhadapt%NVT, size(rhadapt%p_IvertexAge, 1))+1, rhadapt%NVT
      write(UNIT=iunit,FMT=40) -99
    end do

    ! Write variable for element marker
    write(UNIT=iunit,FMT=*) 'elem_mark 0'
    do iel = 1, min(rhadapt%NEL, size(p_Imarker, 1))
      write(UNIT=iunit,FMT=40) p_Imarker(iel)
    end do
    do iel = min(rhadapt%NEL, size(p_Imarker, 1))+1, rhadapt%NEL
      write(UNIT=iunit,FMT=40) -99
    end do

    write(UNIT=iunit,FMT=*) 'endvars'
    write(UNIT=iunit,FMT=*) 'probtime ', 0._DP
    write(UNIT=iunit,FMT=*) 'endgmv'

    ! Close output file
    close(iunit)

10  format(E15.6E3)
15  format(2(1X,I8))
20  format(3(1X,I8))
30  format(4(1X,I8))
40  format(I8)

  end subroutine hadapt_writeGridGMV

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_checkConsistency(rhadapt)

!<description>
    ! This subroutine checks the internal consistency of the dynamic data structures.
    ! Note that this routine performs brute-force search, and hence, should not
    ! be called in a productive environment. In is meant for debugging purposes
    ! only. If an error occurs, it stop without possibility to resume.
!</description>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(inout) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_arraylistInt) :: ralstIter
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer :: ivt,idx,iel,jel,jelmid,ive,jve,nve,mve
    integer :: h_IelementsAtVertexIdx,h_IelementsAtVertex
    logical :: btest,bfound

    ! Test #1: Consistency of element numbers
    btest = (rhadapt%NEL .eq. sum(rhadapt%InelOfType))
    call output_line('Test #1: Checking consistency of element numbers '//&
                     merge('PASSED','FAILED',btest))

    ! Test #2: Vertex age must not exceed maximum refinement level
    btest = .true.
    do ivt = 1, rhadapt%NVT
      btest = btest .or. (rhadapt%p_IvertexAge(ivt) .gt. rhadapt%NSUBDIVIDEMAX)
    end do
    call output_line('Test #2: Checking maximum vertex age '//&
                     merge('PASSED','FAILED',btest))

    ! Test #3: Check consistency of element neighbours
    btest = .true.
    do iel = 1, rhadapt%NEL

      ! Get number of vertices per element
      nve = hadapt_getNVE(rhadapt, iel)

      ! Loop over all adjacent elements
      do ive = 1, nve
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)

        ! Check that adjacent element number is not larger than the
        ! total number of elements present in the triangulation
        if (jel > rhadapt%NEL .or. jelmid > rhadapt%NEL) then
          btest = .false.; cycle
        end if

        ! Do nothing if we are adjacent to the boundary
        if (jel .eq. 0 .or. jelmid .eq. 0) cycle

        ! Is neighboring element subdivided?
        if (jel .eq. jelmid) then

          ! Get number of vertices per element
          mve = hadapt_getNVE(rhadapt,jel)

          ! Find element IEL in adjacency list of JEL
          bfound = .false.
          do jve = 1, mve
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. &
                rhadapt%p_ImidneighboursAtElement(jve, jel)) then
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            else
              if (rhadapt%p_IneighboursAtElement(jve, jel)    .eq. iel .or.&
                  rhadapt%p_ImidneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            end if
          end do
          btest = btest .and. bfound

        else

          ! Get number of vertices per element
          mve = hadapt_getNVE(rhadapt, jel)

          ! Find element IEL in adjacency list of JEL
          bfound = .false.
          do jve = 1, mve
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. &
                rhadapt%p_ImidneighboursAtElement(jve, jel)) then
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            else
              if (rhadapt%p_IneighboursAtElement(jve, jel)    .eq. iel .or.&
                  rhadapt%p_ImidneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            end if
          end do
          btest = btest .and. bfound

          ! Get number of vertices per element
          mve = hadapt_getNVE(rhadapt, jelmid)

          ! Find element IEL in adjacency list of JELMID
          bfound = .false.
          do jve = 1, mve
            if (rhadapt%p_IneighboursAtElement(jve, jelmid) .eq. &
                rhadapt%p_ImidneighboursAtElement(jve, jelmid)) then
              if (rhadapt%p_IneighboursAtElement(jve, jelmid) .eq. iel) then
                bfound = .true.; exit
              end if
            else
              if (rhadapt%p_IneighboursAtElement(jve, jelmid)    .eq. iel .or.&
                  rhadapt%p_ImidneighboursAtElement(jve, jelmid) .eq. iel) then
                bfound = .true.; exit
              end if
            end if
          end do
          btest = btest .and. bfound

        end if
      end do
    end do
    call output_line('Test #3: Checking consistency of element neighbours '//&
                     merge('PASSED','FAILED',btest))

    ! Test #4: Check consistency of common vertices between two edges
    btest = .true.
    do iel = 1, rhadapt%NEL

      ! Get number of vertices per element
      nve = hadapt_getNVE(rhadapt, iel)

      ! Loop over all adjacent elements
      do ive = 1, nve
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)

        ! Check that adjacent element number is not larger than the
        ! total number of elements present in the triangulation
        if (jel > rhadapt%NEL .or. jelmid > rhadapt%NEL) then
          btest = .false.; cycle
        end if

        ! Do nothing if we are adjacent to the boundary
        if (jel .eq. 0 .or. jelmid .eq. 0) cycle

        ! Do nothing if there exists a temporal hanging node
        if (jel .ne. jelmid) cycle

        ! Get number of vertices per element
        mve = hadapt_getNVE(rhadapt, jel)

        ! Find element IEL in adjacency list of JEL
        bfound = .false.
        do jve = 1, mve
          if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. &
              rhadapt%p_ImidneighboursAtElement(jve, jel)) then
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
              bfound = .true.; exit
            end if
          else
            ! Do nothing if there exists a temporal hanging node
            cycle
          end if
        end do

        ! If the common edge has been found, check the two endpoints
        if (bfound) then
          bfound = ((rhadapt%p_IverticesAtElement(ive, iel) .eq. &
                     rhadapt%p_IverticesAtElement(modulo(jve,mve)+1, jel)) .and. &
                     rhadapt%p_IverticesAtElement(modulo(ive,nve)+1, iel) .eq. &
                     rhadapt%p_IverticesAtElement(jve, jel))
        end if
        btest = btest .and. bfound
      end do
    end do
    call output_line('Test #4: Checking consistency of common vertices along edges '//&
                     merge('PASSED','FAILED',btest))

    ! Test #5: Check consistency of element-meeting-at-vertex lists
    btest = (alst_ntable(rhadapt%rElementsAtVertex) .eq. rhadapt%NVT)
    if (btest) then

      ! Create index array
      call storage_new('hadapt_checkConsistency', 'IelementAtVertexIdx', rhadapt%NVT+1,&
                       ST_INT, h_IelementsAtVertexIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)

      ! Count number of elements meeting at vertex
      do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)

        ! Loop over corner vertices
        do ive = 1, nve
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          p_IelementsAtVertexIdx(ivt+1) = p_IelementsAtVertexIdx(ivt+1)+1
        end do
      end do

      ! Convert element couter into absolute position
      p_IelementsAtVertexIdx(1) = 1
      do ivt = 2, rhadapt%NVT+1
        p_IelementsAtVertexIdx(ivt) = p_IelementsAtVertexIdx(ivt)+&
                                      p_IelementsAtVertexIdx(ivt-1)
      end do

      ! Create working array
      call storage_new('hadapt_checkConsistency', 'IelementAtVertex',&
                       p_IelementsAtVertexIdx(rhadapt%NVT+1)-1,&
                       ST_INT, h_IelementsAtVertex, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(h_IelementsAtVertex, p_IelementsAtVertex)

      ! Retrieve the element numbers
      do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)

        ! Loop over corner vertices
        do ive = 1, nve
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          idx = p_IelementsAtVertexIdx(ivt)
          p_IelementsAtVertexIdx(ivt) = idx+1
          p_IelementsAtVertex(idx)    = iel
        end do
      end do

      ! Restore index array
      do ivt = rhadapt%NVT+1, 2, -1
        p_IelementsAtVertexIdx(ivt) = p_IelementsAtVertexIdx(ivt-1)
      end do
      p_IelementsAtVertexIdx(1) = 1

      ! Start to compare the temporal elements-meeting-at-vertex list
      ! and the dynamic data structure from the adaptivity structure
      do ivt = 1, rhadapt%NVT

        ! Get first entry in array list
        ralstIter = alst_begin(rhadapt%rElementsAtVertex, ivt)

        ! Repeat until there is no entry left in the array list
        do while (.not.alst_isNull(ralstIter))

          ! Get element number IEL
          iel = alst_get(rhadapt%rElementsAtVertex, ralstIter)

          ! Proceed to next entry in array list
          call alst_next(ralstIter)

          ! Look for element IEL in temporal elements-meeting-at-vertex list
          ! If it does exist, multiply its value by minus one so that it cannot
          ! be found twice. At the end, all positive entries in the temporal
          ! list are not present in the dynamic structure
          bfound = .false.
          do idx = p_IelementsAtVertexIdx(ivt), p_IelementsAtVertexIdx(ivt+1)-1
            if (p_IelementsAtVertex(idx) .eq. iel) then
              p_IelementsAtVertex(idx) = -p_IelementsAtVertex(idx)
              bfound = .true.
              exit
            end if
          end do
          btest = btest .and. bfound
        end do
      end do
      btest = btest .and. all(p_IelementsAtVertex < 0)
      call output_line('Test #5: Checking consistency of elements meeting at vertices '//&
                       merge('PASSED','FAILED',btest))

      ! Release auxiliary storage
      call storage_free(h_IelementsAtVertexIdx)
      call storage_free(h_IelementsAtVertex)
    else

      call output_line('Test #5: Checking consistency of element-meeting-at-vertex list '//&
                       merge('PASSED','FAILED',btest))
    end if

  end subroutine hadapt_checkConsistency

end module hadaptivity

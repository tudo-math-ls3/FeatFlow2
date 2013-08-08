!##############################################################################
!# ****************************************************************************
!# <name> meshregion </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains structures and routines for managing a specific region
!# inside a given triangulation. These region may be used e.g. for boundary
!# condition assembly.
!#
!# The following routines can be found here:
!#
!#  1.) mshreg_createFromNodalProp
!#      -> Creates a mesh region based on the triangulation`s
!#         nodal property array.
!#
!#  2.) mshreg_createFromHitTest
!#      -> Creates a mesh region based on a hit-test callback function.
!#
!#  3.) mshreg_createFromExpression
!#      -> Creates a mesh region based on an expression using the parser
!#
!#  4.) mshreg_done
!#      -> Releases a mesh region structure.
!#
!#  5.) mshreg_recalcVerticesFromEdges
!#      -> Recalculates the vertex index array of the mesh region from the
!#         edge index array of the mesh region.
!#
!#  6.) mshreg_recalcVerticesFromFaces
!#      -> Recalculates the vertex index array of the mesh region from the
!#         face index array of the mesh region.
!#
!#  7.) mshreg_recalcEdgesFromFaces
!#      -> Recalculates the edge index array of the mesh region from the
!#         face index array of the mesh region.
!#
!#  8.) mshreg_calcBoundaryNormals2D
!#      -> Calculates the outer normal vectors for the boundary edges of a
!#         2D mesh region.
!#
!#  9.) mshreg_calcBoundaryNormals3D
!#      -> Calculates the outer normal vectors for the boundary faces of a
!#         3D mesh region.
!#
!# </purpose>
!##############################################################################

module meshregion

!$use omp_lib
  use fsystem
  use genoutput
  use geometryaux
  use basicgeometry
  use storage
  use triangulation
  use fparser
  use numbersets
  use collection, only: t_collection

  implicit none

  private

  public :: t_meshRegion
  public :: mshreg_createFromNodalProp
  public :: mshreg_createFromHitTest
  public :: mshreg_createFromExpression
  public :: mshreg_done
  public :: mshreg_recalcVerticesFromEdges
  public :: mshreg_recalcVerticesFromFaces
  public :: mshreg_recalcEdgesFromFaces
  public :: mshreg_calcBoundaryNormals2D
  public :: mshreg_calcBoundaryNormals3D

!<constants>

!<constantblock description="Cell type identificator flags">

  ! Identification flag for the vertice index array
  integer(I32), parameter, public :: MSHREG_IDX_VERTEX  = 2**0

  ! Identification flag for the edge index array
  integer(I32), parameter, public :: MSHREG_IDX_EDGE    = 2**1

  ! Identification flag for the face index array
  integer(I32), parameter, public :: MSHREG_IDX_FACE    = 2**2

  ! Identification flag for the element index array
  integer(I32), parameter, public :: MSHREG_IDX_ELEMENT = 2**3

  ! Identification flag for no index arrays
  integer(I32), parameter, public :: MSHREG_IDX_NONE    = 0

  ! Identification flag for all index arrays
  integer(I32), parameter, public :: MSHREG_IDX_ALL     = MSHREG_IDX_VERTEX +&
                      MSHREG_IDX_EDGE + MSHREG_IDX_FACE + MSHREG_IDX_ELEMENT

!</constantblock>

!<constantblock description="Mask operator types">

  ! Kick the old index array
  integer(I32), parameter, public :: MSHREG_MASK_KICK = 0

  ! Apply the OR-operator
  integer(I32), parameter, public :: MSHREG_MASK_OR   = 1

  ! Apply the AND-operator
  integer(I32), parameter, public :: MSHREG_MASK_AND  = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Structure defining a mesh region.
  type t_meshRegion

    ! A pointer to the triangulation on which this mesh region is defined.
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! Share flags of the mesh region. A combination of MSHREG_IDX_XXXX
    ! constants defined above. If a MSHREG_IDX_XXXX flag is defined,
    ! then the mesh region shares the corresponding index array with
    ! another structure.
    integer(I32) :: cshareFlags = 0

    ! Number of vertices in the mesh region.
    integer :: NVT = 0

    ! Number of edges in the mesh region (only 2D/3D).
    integer :: NMT = 0

    ! Number of faces in the mesh region (only 3D).
    integer :: NAT = 0

    ! Number of elements in the mesh region.
    integer :: NEL = 0

    ! Handle to integer array holding the indices of all vertices which
    ! belong to this mesh region.
    integer :: h_IvertexIdx = ST_NOHANDLE

    ! Handle to integer array holding the indices of all edges which
    ! belong to this mesh region (only 2D/3D).
    integer :: h_IedgeIdx = ST_NOHANDLE

    ! Handle to integer array holding the indices of all faces which
    ! belong to this mesh region (only 3D).
    integer :: h_IfaceIdx = ST_NOHANDLE

    ! Handle to integer array holding the indces of all cells which
    ! belong to this mesh region.
    integer :: h_IelementIdx = ST_NOHANDLE

  end type

!</typeblock>

!</types>

  contains

!************************************************************************

!<subroutine>

  subroutine mshreg_aux_calcIdxArray1(IpropArray,ioff,ilen,h_IidxArray,&
                                     inumEntries,iidx,IallowedProp)

!<description>
  ! This auxiliary routine creates an index array from a nodal property
  ! array.
  !
  ! This subroutine is used by some other routines in this module and
  ! is not meant to be called from outside this module, therefore, it is
  ! declared as private.
!</description>

!<input>
  ! The nodal property array from which the index array is to be
  ! calculated.
  integer, dimension(:), intent(in)                :: IpropArray

  ! The offset of the first entry in the nodal property array.
  integer, intent(in)                              :: ioff

  ! The number of entries in the nodal property array.
  integer, intent(in)                              :: ilen

  ! ...
  integer, intent(in)                              :: iidx

  ! OPTIONAL: An array holding all allowed nodal properties that
  ! are to be added into the index array. If not given, all entries
  ! with nodal property > 0 are added.
  integer, dimension(:), optional, intent(in)      :: IallowedProp
!</input>

!<output>
  ! A storage handle to the created index array.
  integer, intent(out)                             :: h_IidxArray

  ! The number of entries in the created index array.
  integer, intent(out)                             :: inumEntries
!</output>

!</subroutine>

  ! Three nice local variables
  integer :: i,j,ientry
  integer, dimension(:), pointer :: p_Iarray

    ! Reset the handle
    h_IidxArray = ST_NOHANDLE

    ! First of all, count the entries
    inumEntries = 0
    if (present(IallowedProp)) then
      do i = 1, ilen
        do j = lbound(IallowedProp,1), ubound(IallowedProp,1)
          if (IpropArray(ioff+i) .eq. IallowedProp(j)) then
            inumEntries = inumEntries + 1
            exit
          end if
        end do
      end do

    else
      do i = 1, ilen
        if (IpropArray(ioff+i) .gt. 0) inumEntries = inumEntries + 1
      end do
    end if

    ! If there are no entries, we can leave here
    if (inumEntries .eq. 0) return

    ! Otherwise allocate an index array
    call storage_new('mshreg_aux_calcIdxArray1', 'p_IidxArray', inumEntries, &
                     ST_INT, h_IidxArray, ST_NEWBLOCK_NOINIT)

    ! Get the index array
    call storage_getbase_int(h_IidxArray, p_Iarray)

    ! Go and calculate the indices
    ientry = 1
    if (present(IallowedProp)) then

      do i = 1, ilen
        do j = lbound(IallowedProp,1), ubound(IallowedProp,1)
          if (IpropArray(ioff+i) .eq. IallowedProp(j)) then
            p_Iarray(ientry) = iidx + i
            ientry = ientry + 1
          end if
        end do
      end do

    else

      do i = 1, ilen
        if (IpropArray(ioff+i) .gt. 0) then
          p_Iarray(ientry) = iidx + i
          ientry = ientry + 1
        end if
      end do

    end if

    ! That is it

  end subroutine mshreg_aux_calcIdxArray1

!************************************************************************

!<subroutine>

  subroutine mshreg_aux_calcIdxArray2(IpropArray,ioff,ilen,h_IidxArray,&
                                     inumEntries,Iindex,IallowedProp)

!<description>
  ! This auxiliary routine creates an index array from a nodal property
  ! array.
  !
  ! This subroutine is used by some other routines in this module and
  ! is not meant to be called from outside this module, therefore, it is
  ! declared as private.
!</description>

!<input>
  ! The nodal property array from which the index array is to be
  ! calculated.
  integer, dimension(:), intent(in)                :: IpropArray

  ! The offset of the first entry in the nodal property array.
  integer, intent(in)                              :: ioff

  ! The number of entries in the nodal property array.
  integer, intent(in)                              :: ilen

  ! An array holding the indices of the entries in the
  ! nodal property array.
  integer, dimension(:), intent(in)                :: Iindex

  ! OPTIONAL: An array holding all allowed nodal properties that
  ! are to be added into the index array. If not given, all entries
  ! with nodal property > 0 are added.
  integer, dimension(:), optional, intent(in)      :: IallowedProp
!</input>

!<output>
  ! A storage handle to the created index array.
  integer, intent(out)                             :: h_IidxArray

  ! The number of entries in the created index array.
  integer, intent(out)                             :: inumEntries
!</output>

!</subroutine>

  ! Three nice local variables
  integer :: i,j,ientry
  integer, dimension(:), pointer :: p_Iarray

    ! Reset the handle
    h_IidxArray = ST_NOHANDLE

    ! First of all, count the entries
    inumEntries = 0
    if (present(IallowedProp)) then
      do i = 1, ilen
        do j = lbound(IallowedProp,1), ubound(IallowedProp,1)
          if (IpropArray(ioff+i) .eq. IallowedProp(j)) then
            inumEntries = inumEntries + 1
            exit
          end if
        end do
      end do

    else
      do i = 1, ilen
        if (IpropArray(ioff+i) .gt. 0) inumEntries = inumEntries + 1
      end do
    end if

    ! If there are no entries, we can leave here
    if (inumEntries .eq. 0) return

    ! Otherwise allocate an index array
    call storage_new('mshreg_aux_calcIdxArray2', 'p_IidxArray', inumEntries, &
                     ST_INT, h_IidxArray, ST_NEWBLOCK_NOINIT)

    ! Get the index array
    call storage_getbase_int(h_IidxArray, p_Iarray)

    ! Go and calculate the indices
    ientry = 1
    if (present(IallowedProp)) then

      do i = 1, ilen
        do j = lbound(IallowedProp,1), ubound(IallowedProp,1)
          if (IpropArray(ioff+i) .eq. IallowedProp(j)) then
            p_Iarray(ientry) = Iindex(i)
            ientry = ientry + 1
          end if
        end do
      end do

    else

      do i = 1, ilen
        if (IpropArray(ioff+i) .gt. 0) then
          p_Iarray(ientry) = Iindex(i)
          ientry = ientry + 1
        end if
      end do

    end if

    ! That is it

  end subroutine mshreg_aux_calcIdxArray2

!************************************************************************

!<subroutine>

  subroutine mshreg_createFromNodalProp(rmeshRegion, rtriangulation, &
                                        cidxCalc, IallowedProp)

!<description>
  ! Creates a new mesh region based on the nodal property array of
  ! the triangulation.
!</description>

!<input>
  ! The triangulation that is to be used.
  type(t_triangulation), target, intent(in)      :: rtriangulation

  ! A combination of MSHREG_IDX_XXXX constants defined above which
  ! specifies which index arrays should be calculated from the nodal
  ! property array. May also be MSHREG_IDX_NONE or MSHREG_IDX_ALL.
  integer(I32), intent(in)                       :: cidxCalc

  ! OPTIONAL: A nodal property array for which the region is to be created.
  ! If present, all vertices, edges and faces which have one of the nodal
  ! properties specifies in the array will be added to this mesh region.
  ! If not present, all vertices, edges and faces which have a positive
  ! nodal property will be added to this mesh region.
  integer, dimension(:), optional, intent(in)    :: IallowedProp
!</input>

!<output>
  ! The mesh region of the triangulation.
  type(t_meshRegion), intent(out)                :: rmeshRegion
!</output>

!</subroutine>

  ! Some local variables
  integer :: idim
  integer, dimension(:), pointer :: p_InodalProperty

    ! Hang in the triangulation
    rmeshRegion%p_rtriangulation => rtriangulation

    ! Get the nodal property array of the triangulation
    call storage_getbase_int(rtriangulation%h_InodalProperty, p_InodalProperty)

    ! Get the dimension of the triangulation
    idim = rtriangulation%ndim

    ! Should we calculate the vertice index array?
    if(iand(cidxCalc, MSHREG_IDX_VERTEX) .ne. 0) then

      ! Calculate the vertice index array
      if (present(IallowedProp)) then

        call mshreg_aux_calcIdxArray1(p_InodalProperty, 0, rtriangulation%NVT, &
                                      rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, &
                                      0, IallowedProp)
      else

        call mshreg_aux_calcIdxArray1(p_InodalProperty, 0, rtriangulation%NVT, &
                                     rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0)
      end if

    end if

    ! Should we calculate the edge index array?
    if((idim .ge. 2) .and. (iand(cidxCalc, MSHREG_IDX_EDGE) .ne. 0)) then

      ! Calculate the edge index array
      if (present(IallowedProp)) then

        call mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT, &
                                      rtriangulation%NMT, rmeshRegion%h_IedgeIdx, &
                                      rmeshRegion%NMT, 0, IallowedProp)
      else

        call mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT, &
                                      rtriangulation%NMT, rmeshRegion%h_IedgeIdx, &
                                      rmeshRegion%NMT, 0)
      end if

    end if

    ! Should we calculate the face index array?
    if((idim .ge. 3) .and. (iand(cidxCalc, MSHREG_IDX_FACE) .ne. 0)) then

      ! Calculate the face index array
      if (present(IallowedProp)) then

        call mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT + &
                                      rtriangulation%NMT, rtriangulation%NAT, &
                                      rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, &
                                      0, IallowedProp)
      else

        call mshreg_aux_calcIdxArray1(p_InodalProperty, rtriangulation%NVT + &
                                      rtriangulation%NMT, rtriangulation%NAT, &
                                      rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, 0)
      end if

    end if

    ! That is it

  end subroutine mshreg_createFromNodalProp

!************************************************************************

!<subroutine>

  subroutine mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
          cidxCalc, bonlyBndry, fmshregHitTest, IallowedHit, rcollection)

!<description>
  ! Creates a new mesh region by calling a hit-test callback function
  ! to decide whether a cell belongs to the mesh region or not.
  !
  ! fmshregHitTest must decide on point coordinates if a point belongs
  ! to a mesh region or not; the returned array Ihit(.) must return a
  ! value > 0 for all mesh entities that should belong to the mesh region.
  !
  ! If IallowedHit is specified, the value in Ihit(.) also counts.
  ! In this case, only those mesh entities are collected in the mesh
  ! region whose Ihit(.) identifier also appears in IallowedHit(.).
  ! This allows to use one subroutine fmshregHitTest for the whole mesh
  ! that returns identifiers for different parts of the mesh.
  ! Using IallowedHit(.) then allows to filter for one or multiple parts.
!</description>

!<input>
  ! The triangulation that is to be used.
  type(t_triangulation), target, intent(in)      :: rtriangulation

  ! A combination of MSHREG_IDX_XXXX constants defined above which
  ! specifies which index arrays should be calculated from the nodal
  ! property array. May also be MSHREG_IDX_NONE or MSHREG_IDX_ALL.
  integer(I32), intent(in)                       :: cidxCalc

  ! If .TRUE., then only the cells which belong to the triangulation`s
  ! boundary are hit-tested, otherwise all cells are hit-tested.
  logical, intent(in)                            :: bonlyBndry

  ! A callback to the hit-test function.
  include 'intf_mshreghittest.inc'

  ! OPTIONAL: Hit-test filter.
  ! If this is not specified, the routine generates a mesh region of
  !   all mesh entities where the callback function fmshregHitTest
  !   returns a value > 0 in Ihit(.).
  ! If this is specified, the routine generates a mesh region of
  !   all mesh entities where the callback function fmshregHitTest
  !   returns a value which also appears in IallowedHit(*).
  integer, dimension(:), optional, intent(in)    :: IallowedHit

!</input>

!<output>
  ! The mesh region of the triangulation.
  type(t_meshRegion), intent(out)                :: rmeshRegion
!</output>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  type(t_collection), optional, intent(inout)    :: rcollection
!</inputoutput>

!</subroutine>

  ! Some local variables
  integer :: idim,i,idx
  integer, dimension(:), pointer :: p_Ihit, p_Iindex
  integer, dimension(:,:), pointer :: p_IvertsAtCell
  real(DP), dimension(:,:), pointer :: p_Dcoords, p_Dverts

    ! Hang in the triangulation
    rmeshRegion%p_rtriangulation => rtriangulation

    ! Get the dimension of the triangulation
    idim = rtriangulation%ndim

    ! Get the vertice coordinates
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Dverts)

    ! Do we calculate vertices?
    if (iand(cidxCalc, MSHREG_IDX_VERTEX) .ne. 0) then

      ! Do we process only boundary vertices?
      if (bonlyBndry) then

        ! Allocate a map
        allocate(p_Ihit(rtriangulation%NVBD))
        
        ! Get the vertices-at-boundary array
        call storage_getbase_int(rtriangulation%h_IverticesAtBoundary, p_Iindex)

        ! Create an array with all vertices on the boundary.
        allocate(p_Dcoords(idim,rtriangulation%NVBD))

        do i=1, rtriangulation%NVBD
          p_Dcoords(:,i) = p_Dverts(:, p_Iindex(i))
        end do

        ! And call the hit-test routine
        call fmshregHitTest(rtriangulation%NVBD, p_Dcoords, p_Ihit, rcollection)
        
        ! Calculate the index array
        call mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NVBD, &
            rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, p_Iindex, IallowedHit)

        ! Deallocate vertex array
        deallocate(p_Dcoords)

        ! Free the hit-map
        deallocate(p_Ihit)

      else

        ! Allocate a map
        allocate(p_Ihit(rtriangulation%NVT))

        ! And call the hit-test routine
        call fmshregHitTest(rtriangulation%NVT, p_Dverts, p_Ihit, rcollection)

        ! Calculate the index array
        call mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NVT, &
            rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0, IallowedHit)

        ! Free the hit-map
        deallocate(p_Ihit)

      end if

    end if

    ! Do we calculate edges?
    if ((idim .ge. 2) .and. (iand(cidxCalc, MSHREG_IDX_EDGE) .ne. 0)) then

      ! Do we process only boundary edgess?
      if (bonlyBndry) then

        ! Allocate a map
        allocate(p_Ihit(rtriangulation%NMBD))

        ! Allocate a coordinate array
        allocate(p_Dcoords(idim,rtriangulation%NMBD))

        ! Get the edges-at-boundary array
        call storage_getbase_int(rtriangulation%h_IedgesAtBoundary, p_Iindex)

        ! Get the vertices-at-edge array
        call storage_getbase_int2D(rtriangulation%h_IverticesAtEdge, &
                                   p_IvertsAtCell)

        ! Calculate the edge mid-points
        do i=1, rtriangulation%NMBD

          ! Get the index of the edge
          idx = p_Iindex(i)

          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.5_DP * (p_Dverts(:, p_IvertsAtCell(1,idx)) &
                                  +  p_Dverts(:, p_IvertsAtCell(2,idx)))
        end do

        ! And call the hit-test routine
        call fmshregHitTest(rtriangulation%NMBD, p_Dcoords, p_Ihit, rcollection)

        ! Calculate the index array
        call mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NMBD, &
            rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, p_Iindex, IallowedHit)

        ! Deallocate midpoint array
        deallocate(p_Dcoords)

        ! Free the hit-map
        deallocate(p_Ihit)

      else

        ! Allocate a coordinate array
        allocate(p_Dcoords(idim,rtriangulation%NMT))

        ! Get the vertices-at-edge array
        call storage_getbase_int2D(rtriangulation%h_IverticesAtEdge, &
                                   p_IvertsAtCell)

        ! Calculate the edge mid-points
        do i=1, rtriangulation%NMBD

          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.5_DP * (p_Dverts(:, p_IvertsAtCell(1,i)) &
                                  +  p_Dverts(:, p_IvertsAtCell(2,i)))
        end do

        ! And call the hit-test routine
        call fmshregHitTest(rtriangulation%NMT, p_Dcoords, p_Ihit, rcollection)

        ! Calculate the index array
        call mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NMT, &
            rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, 0, IallowedHit)

        ! Deallocate midpoint array
        deallocate(p_Dcoords)

        ! Free the hit-map
        deallocate(p_Ihit)

      end if

    end if

    ! Do we calculate face?
    if ((idim .ge. 3) .and. (iand(cidxCalc, MSHREG_IDX_FACE) .ne. 0)) then

      ! Do we process only boundary edgess?
      if (bonlyBndry) then

        ! Allocate a map
        allocate(p_Ihit(rtriangulation%NABD))

        ! Allocate a coordinate array
        allocate(p_Dcoords(idim,rtriangulation%NABD))

        ! Get the faces-at-boundary array
        call storage_getbase_int(rtriangulation%h_IfacesAtBoundary, p_Iindex)

        ! Get the vertices-at-face array
        call storage_getbase_int2D(rtriangulation%h_IverticesAtFace, &
                                   p_IvertsAtCell)

        ! Calculate the face mid-points
        do i=1, rtriangulation%NABD

          ! Get the index of the face
          idx = p_Iindex(i)

          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.25_DP * (p_Dverts(:, p_IvertsAtCell(1,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(2,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(3,idx)) &
                                   +  p_Dverts(:, p_IvertsAtCell(4,idx)))
        end do

        ! And call the hit-test routine
        call fmshregHitTest(rtriangulation%NABD, p_Dcoords, p_Ihit, rcollection)

        ! Calculate the index array
        call mshreg_aux_calcIdxArray2(p_Ihit, 0, rtriangulation%NABD, &
            rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, p_Iindex, IallowedHit)

        ! Deallocate midpoint array
        deallocate(p_Dcoords)

        ! Free the hit-map
        deallocate(p_Ihit)

      else

        ! Allocate a coordinate array
        allocate(p_Dcoords(idim,rtriangulation%NMT))

        ! Get the vertices-at-face array
        call storage_getbase_int2D(rtriangulation%h_IverticesAtEdge, &
                                   p_IvertsAtCell)

        ! Calculate the face mid-points
        do i=1, rtriangulation%NMBD

          ! Calculate the mid-point
          p_Dcoords(:,i) = 0.25_DP * (p_Dverts(:, p_IvertsAtCell(1,i)) &
                                   +  p_Dverts(:, p_IvertsAtCell(2,i)) &
                                   +  p_Dverts(:, p_IvertsAtCell(3,i)) &
                                   +  p_Dverts(:, p_IvertsAtCell(4,i)))
        end do

        ! And call the hit-test routine
        call fmshregHitTest(rtriangulation%NAT, p_Dcoords, p_Ihit, rcollection)

        ! Calculate the index array
        call mshreg_aux_calcIdxArray1(p_Ihit, 0, rtriangulation%NAT, &
            rmeshRegion%h_IfaceIdx, rmeshRegion%NAT, 0, IallowedHit)

        ! Deallocate midpoint array
        deallocate(p_Dcoords)

        ! Free the hit-map
        deallocate(p_Ihit)

      end if

    end if

    ! That is it

  end subroutine mshreg_createFromHitTest

!************************************************************************

!<subroutine>

  subroutine mshreg_createFromExpression(rmeshRegion, rtriangulation, &
          cidxCalc, bonlyBndry, scondition, SadditionalVars, DadditionalVars,&
          IallowedHit)

!<description>
  ! Creates a new mesh region by evaluating a condition by a parser
  ! to decide whether a cell belongs to the mesh region or not.
  !
  ! The string scondition represents a parser condition which returns a value
  ! >0 (integer based, =1, e.g.) if a point belongs to a mesh region
  ! or =0, if not.
  !
  ! The routine provides the variables "X", "Y" and "Z" which represent
  ! the X/Y/Z coordinate of the point to check. Via SadditionalVars (names) and
  ! DadditionalVars (values), the user can pass additional variables to the 
  ! condition. The variables have a limit of 10 characters. The number of
  ! possible additional variables is limited to 125.
  !
  ! If IallowedHit is specified, the value returned by scondition also counts.
  ! In this case, only those mesh entities are collected in the mesh
  ! region where the value returned by scondition also appears in IallowedHit(.).
  ! This allows to use one condition for the whole mesh
  ! that returns identifiers for different parts of the mesh.
  ! Using IallowedHit(.) then allows to filter for one or multiple parts.
!</description>

!<remarks>
  ! fparser_init() must have been called prior to this routine!
!</remarks>

!<input>
  ! The triangulation that is to be used.
  type(t_triangulation), target, intent(in)      :: rtriangulation

  ! A combination of MSHREG_IDX_XXXX constants defined above which
  ! specifies which index arrays should be calculated from the nodal
  ! property array. May also be MSHREG_IDX_NONE or MSHREG_IDX_ALL.
  integer(I32), intent(in)                       :: cidxCalc

  ! If .TRUE., then only the cells which belong to the triangulation`s
  ! boundary are hit-tested, otherwise all cells are hit-tested.
  logical, intent(in)                            :: bonlyBndry

  ! Condition
  character(len=*), intent(in) :: scondition
  
  ! OPTINAL: Additional variables passed to the condition.
  character(len=*), dimension(:), intent(in), optional :: SadditionalVars
  
  ! OPTIONAL: Values of the variables. Must be present if SadditionalVars
  ! is specified.
  real(DP), dimension(:), intent(in), optional :: DadditionalVars

  ! OPTIONAL: Hit-test filter.
  ! If this is not specified, the routine generates a mesh region of
  !   all mesh entities where the callback function fmshregHitTest
  !   returns a value > 0 in Ihit(.).
  ! If this is specified, the routine generates a mesh region of
  !   all mesh entities where the callback function fmshregHitTest
  !   returns a value which also appears in IallowedHit(*).
  integer, dimension(:), optional, intent(in)    :: IallowedHit

!</input>

!<output>
  ! The mesh region of the triangulation.
  type(t_meshRegion), intent(out)                :: rmeshRegion
!</output>

!</subroutine>

  ! Some local variables
  character(len=10), dimension(:), allocatable :: Svars
  type(t_collection) :: rcollection
  integer :: i
  type(t_fparser), target :: rparser
  
    ! Prepare a parser.
    call fparser_create (rparser, 1)
    
    ! Collect the variables.
    if (.not. present(SadditionalVars)) then
      allocate (Svars(3))
    else
      allocate (Svars(3+ubound(SadditionalVars,1)))
      do i=1,ubound(SadditionalVars,1)
        Svars(3+i-1) = SadditionalVars(i)
      end do
    end if

    ! X/Y/Z variables are the first
    Svars(1) = "X"
    Svars(2) = "Y"
    Svars(3) = "Z"
    
    ! Parse the condition
    call fparser_parseFunction (rparser, 1, scondition, Svars)
    
    ! Release variable names.
    deallocate(Svars)
    
    ! Copy provided variables to the collection
    rcollection%DquickAccess(:) = 0
    if (present(DadditionalVars)) then
      rcollection%DquickAccess(4:4+ubound(DadditionalVars,1)-1) = DadditionalVars(:)
    end if
    
    ! Provide the parser in the collection.
    rcollection%p_rfparserQuickAccess1 => rparser
    
    ! Call the previous hit-test routine, provide a specialised expression
    ! evaluation function.
    call mshreg_createFromHitTest(rmeshRegion, rtriangulation, &
        cidxCalc, bonlyBndry, mshreg_fEvaluateExpression, IallowedHit, rcollection)
        
    ! Release the parser
    call fparser_release(rparser)

  end subroutine mshreg_createFromExpression
  
!************************************************************************

!<subroutine>

  subroutine mshreg_fEvaluateExpression(inumCells,Dcoords,Ihit,rcollection)

  use collection
  use fsystem

!<description>
  ! This subroutine is called for hit-testing cells to decide whether a
  ! cell belongs to a specific mesh region or not.
!</description>

!<input>
  ! Number of cells (i.e. vertices,edges,faces,etc.) for which the
  ! hit-test is to be performed.
  integer, intent(in)                          :: inumCells

  ! Coordinates of the points for which the hit-test is to be performed.
  ! The dimension of the array is at least (1:idim,1:inumCells), where:
  ! -> idim is the dimension of the mesh, i.e. 1 for 1D, 2 for 2D, etc.
  ! -> inumCells is the parameter passed to this routine.
  real(DP), dimension(:,:), intent(in)         :: Dcoords
!</input>

!<output>
  ! An array that recieves the result of the hit-test.
  ! The dimension of the array is at least (1:inumCells).
  integer, dimension(:), intent(out)           :: Ihit
!</output>

!</inputoutput>
  ! OPTIONAL: A collection structure to provide additional information
  ! to the hit-test routine.
  type(t_collection), intent(inout), optional  :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dresult
    integer :: i,n

    n = min(NDIM3D,ubound(Dcoords,1))

    do i=1,inumCells
      ! Copy the coordinates to the beginning of DquickAccess
      rcollection%DquickAccess(1:n) = Dcoords(1:n,i)
    
      ! For all points, evaluate the parser.
      call fparser_evalFunction (&
          rcollection%p_rfparserQuickAccess1, 1, rcollection%DquickAccess, dresult)
          
      ! Convert to integer, return as hit test ID.
      Ihit(i) = nint(dresult)
    end do

  end subroutine

!************************************************************************

!<subroutine>

  subroutine mshreg_done(rmeshRegion)

!<description>
  ! Releases a mesh region structure and frees all allocated memory from
  ! the heap.
!</description>

!<inputoutput>
  type(t_meshRegion), intent(inout) :: rmeshRegion
!</inputoutput>

!</subroutine>

    ! If the region does not share its vertex index array, release it
    if ((rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) .and. &
        (iand(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .eq. 0)) then

        call storage_free(rmeshRegion%h_IvertexIdx)

    end if

    ! If the region does not share its edge index array, release it
    if ((rmeshRegion%h_IedgeIdx .ne. ST_NOHANDLE) .and. &
        (iand(rmeshRegion%cshareFlags, MSHREG_IDX_EDGE) .eq. 0)) then

        call storage_free(rmeshRegion%h_IedgeIdx)

    end if

    ! If the region does not share its face index array, release it
    if ((rmeshRegion%h_IfaceIdx .ne. ST_NOHANDLE) .and. &
        (iand(rmeshRegion%cshareFlags, MSHREG_IDX_FACE) .eq. 0)) then

        call storage_free(rmeshRegion%h_IfaceIdx)

    end if

    ! If the region does not share its element index array, release it
    if ((rmeshRegion%h_IelementIdx .ne. ST_NOHANDLE) .and. &
        (iand(rmeshRegion%cshareFlags, MSHREG_IDX_ELEMENT) .eq. 0)) then

        call storage_free(rmeshRegion%h_IelementIdx)

    end if

    ! Last but not least: reset the structure
    rmeshRegion%p_rtriangulation => null()
    rmeshRegion%cshareFlags = 0
    rmeshRegion%NVT = 0
    rmeshRegion%NMT = 0
    rmeshRegion%NAT = 0
    rmeshRegion%NEL = 0
    ! Although usually the handles are reset by the storage_free() routine,
    ! we need to reset the handles by hand, too, because if a handle was given,
    ! and the array represented by the handle was marked as shared in the flags,
    ! the storage_free() routine is not called.
    rmeshRegion%h_IvertexIdx = ST_NOHANDLE
    rmeshRegion%h_IedgeIdx = ST_NOHANDLE
    rmeshRegion%h_IfaceIdx = ST_NOHANDLE
    rmeshRegion%h_IelementIdx = ST_NOHANDLE

    ! That is it

  end subroutine mshreg_done

!************************************************************************

!<subroutine>

  subroutine mshreg_recalcVerticesFromEdges(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the vertice index array from the edge index array.
  ! The operator mask parameter decides which vertices will be added into
  ! the mesh region`s vertice index array:
  ! 1. coperatorMask = MSHREG_MASK_KICK:
  !    The old vertice index array is deleted (if it exists), and every
  !    vertice that belongs to an edge in the mesh region will be added to
  !    the vertice index array.
  !
  ! 2. coperatorMask = MSHREG_MASK_OR:
  ! -> Every vertice that already exists in the mesh region or that belongs
  !    to an edge in the mesh region will be added to vertice index array.
  !
  ! 3. coperatorMask = MSHREG_MASK_AND:
  ! -> Every vertice that already exists in the mesh region and that belongs
  !    to an edge in the mesh region will be added to vertice index array.
  !    This implies that if the mesh region`s vertice index array is empty
  !    on entry, it will also be empty on exit.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old vertice index array.
  ! One of the MSHREG_MASK_XXXX constants defined above.
  ! If not given, MSHREG_MASK_KICK is used.
  integer(I32), intent(in), optional   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  type(t_meshRegion), intent(inout)    :: rmeshRegion
!</inputoutput>

!</subroutine>

  integer, dimension(:), allocatable :: Imap
  integer, dimension(:), pointer :: p_IedgeIdx, p_IvertIdx
  integer, dimension(:,:), pointer :: p_IvertsAtEdge
  integer, dimension(1) :: IallowedProp
  integer(I32) :: copMask
  integer :: i, iedge
  type(t_triangulation), pointer :: p_rtria

    ! Decide what operator to use
    copMask = MSHREG_MASK_KICK
    if (present(coperatorMask)) copMask = coperatorMask

    ! If we use the AND-operator, the allowed property is 2, otherwise 1
    if (copMask .eq. MSHREG_MASK_AND) then
      IallowedProp(1) = 2
    else
      IallowedProp(1) = 1
    end if

    ! Do we have any edges at all?
    if ((rmeshRegion%NMT .le. 0) .or. (rmeshRegion%h_IedgeIdx .eq. &
         ST_NOHANDLE)) then

      ! If the mesh region already had an vertex index array and we do not
      ! apply the OR-operator, kick the old vertex index array.
      if (copMask .ne. MSHREG_MASK_OR) then

        if ((rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) .and. &
            (iand(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .eq. 0)) then

            call storage_free(rmeshRegion%h_IvertexIdx)

        end if

        ! Remove anything that has to do with vertices from the mesh region
        rmeshRegion%NVT = 0
        rmeshRegion%h_IvertexIdx = ST_NOHANDLE
        rmeshRegion%cshareFlags = iand(rmeshRegion%cshareFlags, &
                                       not(MSHREG_IDX_VERTEX))

      end if

      ! Return here
      return

    end if

    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation

    ! Get the edge index array
    call storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)

    ! Get the vertices-at-edge array from the triangulation
    call storage_getbase_int2D(p_rtria%h_IverticesAtEdge, p_IvertsAtEdge)

    ! Allocate an vertex map
    allocate(Imap(p_rtria%NVT))
    do i=1, p_rtria%NVT
      Imap(i) = 0
    end do

    ! Go through all edges
    do i=1, rmeshRegion%NMT

      ! Get the edge index
      iedge = p_IedgeIdx(i)

      ! Go through all vertices adjacent to that edge
      Imap(p_IvertsAtEdge(1,iedge)) = 1
      Imap(p_IvertsAtEdge(2,iedge)) = 1

    end do

    ! Now let us see what to do with the old vertices
    if (rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) then

      ! Get the old vertices index array
      call storage_getbase_int(rmeshRegion%h_IvertexIdx, p_IvertIdx)

      select case(copMask)
      case (MSHREG_MASK_OR)
        ! Go through all vertices and set the map entry to 1
        do i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = 1
        end do

      case (MSHREG_MASK_AND)
        ! Go through all vertices and add 1 to the map entry
        do i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = Imap(p_IvertIdx(i)) + 1
        end do

        ! Now all vertices that have already been in the mesh region
        ! and which are adjacent to an edge in the mesh region have
        ! a value of 2 in the map.

      end select

      ! Now let us destroy the old vertex index array
      if ((rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) .and. &
          (iand(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .eq. 0)) then

          call storage_free(rmeshRegion%h_IvertexIdx)

      end if

      ! Remove anything that has to do with vertices from the mesh region
      rmeshRegion%NVT = 0
      rmeshRegion%h_IvertexIdx = ST_NOHANDLE
      rmeshRegion%cshareFlags = iand(rmeshRegion%cshareFlags, &
                                     not(MSHREG_IDX_VERTEX))

    end if

    ! Now we have the vertex index map, so create an index array from this.
    call mshreg_aux_calcIdxArray1(Imap, 0, p_rtria%NVT, &
        rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0, IallowedProp)

    ! And release the vertex map
    deallocate(Imap)

    ! That is it

  end subroutine mshreg_recalcVerticesFromEdges

!************************************************************************

!<subroutine>

  subroutine mshreg_recalcVerticesFromFaces(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the vertice index array from the face index array.
  ! The operator mask parameter decides which vertices will be added into
  ! the mesh region`s vertice index array:
  ! 1. coperatorMask = MSHREG_MASK_KICK:
  !    The old vertice index array is deleted (if it exists), and every
  !    vertice that belongs to a face in the mesh region will be added to
  !    the vertice index array.
  !
  ! 2. coperatorMask = MSHREG_MASK_OR:
  ! -> Every vertice that already exists in the mesh region or that belongs
  !    to a face in the mesh region will be added to vertice index array.
  !
  ! 3. coperatorMask = MSHREG_MASK_AND:
  ! -> Every vertice that already exists in the mesh region and that belongs
  !    to a face in the mesh region will be added to vertice index array.
  !    This implies that if the mesh region`s vertice index array is empty
  !    on entry, it will also be empty on exit.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old vertice index array.
  ! One of the MSHREG_MASK_XXXX constants defined above.
  ! If not given, MSHREG_MASK_KICK is used.
  integer(I32), intent(in), optional   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  type(t_meshRegion), intent(inout)    :: rmeshRegion
!</inputoutput>

!</subroutine>

  integer, dimension(:), allocatable :: Imap
  integer, dimension(:), pointer :: p_IfaceIdx, p_IvertIdx
  integer, dimension(:,:), pointer :: p_IvertsAtFace
  integer, dimension(1) :: IallowedProp
  integer(I32) :: copMask
  integer :: i, j, iface
  type(t_triangulation), pointer :: p_rtria

    ! Decide what operator to use
    copMask = MSHREG_MASK_KICK
    if (present(coperatorMask)) copMask = coperatorMask

    ! If we use the AND-operator, the allowed property is 2, otherwise 1
    if (copMask .eq. MSHREG_MASK_AND) then
      IallowedProp(1) = 2
    else
      IallowedProp(1) = 1
    end if

    ! Do we have any faces at all?
    if ((rmeshRegion%NAT .le. 0) .or. (rmeshRegion%h_IfaceIdx .eq. &
         ST_NOHANDLE)) then

      ! If the mesh region already had an vertex index array and we do not
      ! apply the OR-operator, kick the old vertex index array.
      if (copMask .ne. MSHREG_MASK_OR) then

        if ((rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) .and. &
            (iand(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .eq. 0)) then

            call storage_free(rmeshRegion%h_IvertexIdx)

        end if

        ! Remove anything that has to do with vertices from the mesh region
        rmeshRegion%NVT = 0
        rmeshRegion%h_IvertexIdx = ST_NOHANDLE
        rmeshRegion%cshareFlags = iand(rmeshRegion%cshareFlags, &
                                  not(MSHREG_IDX_VERTEX))

      end if

      ! Return here
      return

    end if

    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation

    ! Get the face index array
    call storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)

    ! Get the vertices-at-face array from the triangulation
    call storage_getbase_int2D(p_rtria%h_IverticesAtFace,&
                               p_IvertsAtFace)

    ! Allocate an vertex map
    allocate(Imap(p_rtria%NVT))
    do i=1, p_rtria%NVT
      Imap(i) = 0
    end do

    ! Go through all faces
    do i=1, rmeshRegion%NAT

      ! Get the face index
      iface = p_IfaceIdx(i)

      ! Go through all edges adjacent to that face
      do j=1, ubound(p_IvertsAtFace,1)

        if (p_IvertsAtFace(j,iface) .gt. 0) then
          Imap(p_IvertsAtFace(j,iface)) = 1
        end if

      end do

    end do

    ! Now let us see what to do with the old vertices
    if (rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) then

      ! Get the old vertices index array
      call storage_getbase_int(rmeshRegion%h_IvertexIdx, p_IvertIdx)

      select case(copMask)
      case (MSHREG_MASK_OR)
        ! Go through all vertices and set the map entry to 1
        do i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = 1
        end do

      case (MSHREG_MASK_AND)
        ! Go through all vertices and add 1 to the map entry
        do i=1, rmeshRegion%NVT
          Imap(p_IvertIdx(i)) = Imap(p_IvertIdx(i)) + 1
        end do

        ! Now all vertices that have already been in the mesh region
        ! and which are adjacent to a face in the mesh region have
        ! a value of 2 in the map.

      end select

      ! Now let us destroy the old vertex index array
      if ((rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) .and. &
          (iand(rmeshRegion%cshareFlags, MSHREG_IDX_VERTEX) .eq. 0)) then

          call storage_free(rmeshRegion%h_IvertexIdx)

      end if

      ! Remove anything that has to do with vertices from the mesh region
      rmeshRegion%NVT = 0
      rmeshRegion%h_IvertexIdx = ST_NOHANDLE
      rmeshRegion%cshareFlags = iand(rmeshRegion%cshareFlags, &
                                     not(MSHREG_IDX_VERTEX))

    end if

    ! Now we have the vertex index map, so create an index array from this.
    call mshreg_aux_calcIdxArray1(Imap, 0, p_rtria%NVT, &
        rmeshRegion%h_IvertexIdx, rmeshRegion%NVT, 0, IallowedProp)

    ! And release the vertex map
    deallocate(Imap)

    ! That is it

  end subroutine mshreg_recalcVerticesFromFaces

!************************************************************************

!<subroutine>

  subroutine mshreg_recalcEdgesFromFaces(rmeshRegion,coperatorMask)

!<description>
  ! Recalculates the edge index array from the face index array.
  ! The operator mask parameter decides which edges will be added into
  ! the mesh region`s edge index array:
  ! 1. coperatorMask = MSHREG_MASK_KICK:
  !    The old edge index array is deleted (if it exists), and every
  !    edge that belongs to a face in the mesh region will be added to
  !    the edge index array.
  !
  ! 2. coperatorMask = MSHREG_MASK_OR:
  ! -> Every edge that already exists in the mesh region or that belongs
  !    to a face in the mesh region will be added to edge index array.
  !
  ! 3. coperatorMask = MSHREG_MASK_AND:
  ! -> Every edge that already exists in the mesh region and that belongs
  !    to a face in the mesh region will be added to edge index array.
  !    This implies that if the mesh region`s edge index array is empty
  !    on entry, it will also be empty on exit.
!</description>

!<input>
  ! OPTIONAL: An operator mask for the old edge index array.
  ! One of the MSHREG_MASK_XXXX constants defined above.
  ! If not given, MSHREG_MASK_KICK is used.
  integer(I32), intent(in), optional   :: coperatorMask
!</input>

!<inputoutput>
  ! The mesh region.
  type(t_meshRegion), intent(inout)    :: rmeshRegion
!</inputoutput>

!</subroutine>

  integer, dimension(:), allocatable :: Imap
  integer, dimension(:), pointer :: p_IfaceIdx, p_IedgeIdx
  integer, dimension(:,:), pointer :: p_IedgesAtFace
  integer, dimension(1) :: IallowedProp
  integer(I32) :: copMask
  integer :: i, j, iface
  type(t_triangulation), pointer :: p_rtria

    ! Decide what operator to use
    copMask = MSHREG_MASK_KICK
    if (present(coperatorMask)) copMask = coperatorMask

    ! If we use the AND-operator, the allowed property is 2, otherwise 1
    if (copMask .eq. MSHREG_MASK_AND) then
      IallowedProp(1) = 2
    else
      IallowedProp(1) = 1
    end if

    ! Do we have any faces at all?
    if ((rmeshRegion%NAT .le. 0) .or. (rmeshRegion%h_IfaceIdx .eq. &
         ST_NOHANDLE)) then

      ! If the mesh region already had an edge index array and we do not
      ! apply the OR-operator, kick the old edge index array.
      if (copMask .ne. MSHREG_MASK_OR) then

        if ((rmeshRegion%h_IedgeIdx .ne. ST_NOHANDLE) .and. &
            (iand(rmeshRegion%cshareFlags, MSHREG_IDX_EDGE) .eq. 0)) then

            call storage_free(rmeshRegion%h_IedgeIdx)

        end if

        ! Remove anything that has to do with edges from the mesh region
        rmeshRegion%NMT = 0
        rmeshRegion%h_IedgeIdx = ST_NOHANDLE
        rmeshRegion%cshareFlags = iand(rmeshRegion%cshareFlags, &
                                       not(MSHREG_IDX_EDGE))

      end if

      ! Return here
      return

    end if

    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation

    ! Get the face index array
    call storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)

    ! Get the edges-at-face array from the triangulation
    call storage_getbase_int2D(p_rtria%h_IedgesAtFace, p_IedgesAtFace)

    ! Allocate an edge map
    allocate(Imap(p_rtria%NMT))
    do i = 1, p_rtria%NMT
      Imap(i) = 0
    end do

    ! Go through all faces
    do i=1, rmeshRegion%NAT

      ! Get the face index
      iface = p_IfaceIdx(i)

      ! Go through all edges adjacent to that face
      do j=1, ubound(p_IedgesAtFace,1)

        if (p_IedgesAtFace(j,iface) .gt. 0) then
          Imap(p_IedgesAtFace(j,iface)) = 1
        end if

      end do

    end do

    ! Now let us see what to do with the old edges
    if (rmeshRegion%h_IedgeIdx .ne. ST_NOHANDLE) then

      ! Get the old edge index array
      call storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)

      select case(copMask)
      case (MSHREG_MASK_OR)
        ! Go through all edges and set the map entry to 1
        do i=1, rmeshRegion%NMT
          Imap(p_IedgeIdx(i)) = 1
        end do

      case (MSHREG_MASK_AND)
        ! Go through all edges and add 1 to the map entry
        do i=1, rmeshRegion%NMT
          Imap(p_IedgeIdx(i)) = Imap(p_IedgeIdx(i)) + 1
        end do

        ! Now all edges that have already been in the mesh region
        ! and which are adjacent to a face in the mesh region have
        ! a value of 2 in the map.

      end select

      ! Now let us destroy the old edge index array
      if (iand(rmeshRegion%cshareFlags, MSHREG_IDX_EDGE) .eq. 0) then
          call storage_free(rmeshRegion%h_IedgeIdx)
      end if

      ! Remove anything that has to do with edges from the mesh region
      rmeshRegion%NMT = 0
      rmeshRegion%h_IedgeIdx = ST_NOHANDLE
      rmeshRegion%cshareFlags = iand(rmeshRegion%cshareFlags, &
                                     not(MSHREG_IDX_EDGE))

    end if

    ! Now we have the edge index map, so create an index array from this.
    call mshreg_aux_calcIdxArray1(Imap, 0, p_rtria%NMT, &
        rmeshRegion%h_IedgeIdx, rmeshRegion%NMT, 0, IallowedProp)

    ! And release the edge map
    deallocate(Imap)

    ! That is it

  end subroutine mshreg_recalcEdgesFromFaces

!************************************************************************

!<subroutine>

  subroutine mshreg_calcBoundaryNormals2D(rmeshRegion,Dnormals)

!<description>
  ! Calculates the outer normal vectors for the edges of a 2D mesh
  ! region, i.e. the normal vectors which point "out of the domain".
  ! It is silently assumed that all edges of the mesh region lie on
  ! the boundary of the mesh.
!</description>

!<input>
  ! The mesh region for whose edges the normals are to be calcuated.
  type(t_meshRegion), intent(in)        :: rmeshRegion
!</input>

!<output>
  ! An array that recieves the normal vectors of the edges.
  ! Dimension: Dnormals(NDIM2D,rmeshRegion%NMT)
  real(DP), dimension(:,:), intent(out) :: Dnormals
!</output>

!</subroutine>

  ! Some local variables
  type(t_triangulation), pointer :: p_rtria
  real(DP), dimension(:,:), pointer :: p_Dcoords
  integer, dimension(:), pointer :: p_IedgeIdx
  integer, dimension(:,:), pointer :: p_IvertAtEdge
  real(DP) :: dx, dy, dt
  integer :: imt, iedge


    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation

    ! Is this a 2D triangulation?
    if (p_rtria%ndim .ne. 2) then
      call output_line('Triangulation must be 2D!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'mshreg_calcBoundaryNormals2D')
      call sys_halt()
    end if

    ! Do not we have any edges here?
    if (rmeshRegion%NMT .le. 0) return

    ! In 2D the task of calculating the normal vectors is quite easy,
    ! if we exploit some knowledge about the algorithm which calculates
    ! the edges in a 2D mesh.
    ! As we know that all edges on the boundary are oriented in mathematically
    ! positive direction, the only thing that we have to do is rotate the
    ! edge vector by 90deg counter-clock-wise (and normalise, of course)
    ! to get the inner normal vector of a boundary edge...

    ! Get a pointer to the vertices-at-edge array
    call storage_getbase_int2D(p_rtria%h_IverticesAtEdge, p_IvertAtEdge)

    ! Get the vertice coordinate array
    call storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dcoords)

    ! And get the edge index array of the mesh region
    call storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)

    ! Now loop through all edges in the mesh region
    do imt = 1, rmeshRegion%NMT

      ! Get the index of the edge
      iedge = p_IedgeIdx(imt)

      ! Calculate edge
      dx = p_Dcoords(1,p_IvertAtEdge(2,iedge)) &
         - p_Dcoords(1,p_IvertAtEdge(1,iedge))
      dy = p_Dcoords(2,p_IvertAtEdge(2,iedge)) &
         - p_Dcoords(2,p_IvertAtEdge(1,iedge))

      ! Calculate length of edge
      dt = 1.0_DP / sqrt(dx**2 + dy**2)

      ! Calculate normal
      Dnormals(1,imt) = -dy * dt
      Dnormals(2,imt) =  dx * dt

    end do

    ! That is it

  end subroutine mshreg_calcBoundaryNormals2D


!************************************************************************

!<subroutine>

  subroutine mshreg_calcBoundaryNormals3D(rmeshRegion,Dnormals)

!<description>
  ! Calculates the outer normal vectors for the faces of a 3D mesh
  ! region, i.e. the normal vectors which point "out of the domain".
  ! It is silently assumed that all faces of the mesh region lie on
  ! the boundary of the mesh.
!</description>

!<input>
  ! The mesh region for whose edges the normals are to be calcuated.
  type(t_meshRegion), intent(in)        :: rmeshRegion
!</input>

!<output>
  ! An array that recieves the normal vectors of the faces.
  ! Dimension: Dnormals(NDIM3D,rmeshRegion%NAT)
  real(DP), dimension(:,:), intent(out) :: Dnormals
!</output>

!</subroutine>

  ! Some local variables
  type(t_triangulation), pointer :: p_rtria
  real(DP), dimension(:,:), pointer :: p_Dcoords
  integer, dimension(:), pointer :: p_IfaceIdx
  integer, dimension(:,:), pointer :: p_IvertAtElem,&
    p_IfaceAtElem, p_IelemAtFace
  real(DP), dimension(3) :: Du,Dv,Dn
  real(DP), dimension(3,8) :: Dcorners
  real(DP) :: dt
  integer :: iat,ivt,iel,ifae,iface


    ! Get a pointer to the triangulation
    p_rtria => rmeshRegion%p_rtriangulation

    ! Is this a 3D triangulation?
    if (p_rtria%ndim .ne. 3) then
      call output_line('Triangulation must be 3D!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'mshreg_calcBoundaryNormals3D')
      call sys_halt()
    end if

    ! Do not we have any faces here?
    if (rmeshRegion%NAT .le. 0) return

    ! Get a pointer to the faces-at-element array
    call storage_getbase_int2D(p_rtria%h_IfacesAtElement, p_IfaceAtElem)

    ! Get a pointer to the elements-at-face array
    call storage_getbase_int2D(p_rtria%h_IelementsAtFace, p_IelemAtface)

    ! Get a pointer to the vertices-at-element array
    call storage_getbase_int2D(p_rtria%h_IverticesAtElement, p_IvertAtElem)

    ! Get the vertice coordinate array
    call storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dcoords)

    ! And get the edge index array of the mesh region
    call storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)

    ! Now loop through all face in the mesh region
    do iat = 1, rmeshRegion%NAT

      ! Get the index of the face
      iface = p_IfaceIdx(iat)

      ! Get the index of the element that is adjacent to that face
      iel = p_IelemAtFace(1,iface)

      ! Now go through all faces of the element and search for this one
      do ifae = 1, 6
        if (p_IfaceAtElem(ifae,iel) .eq. iface) exit
      end do

      ! Get the eight corner vertices of the hexahedron
      do ivt = 1, 8
        Dcorners(1:3,ivt) = p_Dcoords(1:3,p_IvertAtElem(ivt,iel))
      end do

      ! Calculate the tangentials of the face
      select case(ifae)
      case (1)
        ! (1->3) x (4->2)
        Du(:) = Dcorners(:,3) - Dcorners(:,1)
        Dv(:) = Dcorners(:,2) - Dcorners(:,4)
      case (2)
        ! (1->6) x (2->5)
        Du(:) = Dcorners(:,6) - Dcorners(:,1)
        Dv(:) = Dcorners(:,5) - Dcorners(:,2)
      case (3)
        ! (2->7) x (6->3)
        Du(:) = Dcorners(:,7) - Dcorners(:,2)
        Dv(:) = Dcorners(:,3) - Dcorners(:,6)
      case (4)
        ! (3->8) x (4->7)
        Du(:) = Dcorners(:,8) - Dcorners(:,3)
        Dv(:) = Dcorners(:,7) - Dcorners(:,4)
      case (5)
        ! (4->5) x (1->8)
        Du(:) = Dcorners(:,5) - Dcorners(:,4)
        Dv(:) = Dcorners(:,8) - Dcorners(:,1)
      case (6)
        ! (5->7) x (6->8)
        Du(:) = Dcorners(:,7) - Dcorners(:,5)
        Dv(:) = Dcorners(:,8) - Dcorners(:,6)
      end select

      ! Calculate normal by 3D cross product
      Dn(1) = Du(2)*Dv(3) - Du(3)*Dv(2)
      Dn(2) = Du(3)*Dv(1) - Du(1)*Dv(3)
      Dn(3) = Du(1)*Dv(2) - Du(2)*Dv(1)
      dt = 1.0_DP / sqrt(Dn(1)**2 + Dn(2)**2 + Dn(3)**2)
      Dn = dt * Dn

      ! Store the normal
      if (gaux_isFlipped_hexa3D(Dcorners)) then
        Dnormals(:,iat) = Dn(:)
      else
        Dnormals(:,iat) = -Dn(:)
      end if

    end do

    ! That is it

  end subroutine mshreg_calcBoundaryNormals3D

!************************************************************************

!<subroutine>

  subroutine mshreg_auxMarkElementsAtVerts(rtriangulation,Ivertices,rset)

!<description>
  ! Marks all elements adjacent to the vertices in Ivertices.
!</description>

!<input>
  ! Underlying triangulation.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of vertices
  integer, dimension(:), intent(in) :: Ivertices
!</input>

!<inputoutput>
  ! Set where all elements adjacent to the vertices are marked.
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer :: i,ivt,iel1,iel2
    
    ! Get data from the triangulation
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)
    call storage_getbase_int (rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
    
    ! Loop through the vertices
    do i=1,size(Ivertices)
      ivt = Ivertices(i)
      iel1 = p_IelementsAtVertexIdx(ivt)
      iel2 = p_IelementsAtVertexIdx(ivt+1)-1
      
      ! Add the adjacent elements to the set
      call nsets_putElements (rset,p_IelementsAtVertex(iel1:iel2))
    end do

  end subroutine

!************************************************************************

!<subroutine>

  subroutine mshreg_auxMarkElementsAtEdges(rtriangulation,Iedges,rset)

!<description>
  ! Marks all elements adjacent to the edges in Iedges.
!</description>

!<input>
  ! Underlying triangulation.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of edges
  integer, dimension(:), intent(in) :: Iedges
!</input>

!<inputoutput>
  ! Set where all elements adjacent to the edges are marked.
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: idim
    integer, dimension(:), pointer :: p_IelementsAtEdgeIdx
    integer, dimension(:), pointer :: p_IelementsAtEdge
    integer, dimension(:,:), pointer :: p_IelementsAtEdge2D
    integer :: i,imt,iel1,iel2

    ! Get the dimension of the triangulation
    idim = rtriangulation%ndim
    
    select case (idim)
    case (NDIM2D)
      ! Get data from the triangulation
      call storage_getbase_int2d (rtriangulation%h_IelementsAtEdge,p_IelementsAtEdge2D)
      
      ! Add adjacent elements
      do i=1,size(Iedges)
        imt = Iedges(i)
      
        ! The first element is always there, the 2nd element not necessarily
        call nsets_addToDASet (rset,p_IelementsAtEdge2D(1,imt))
        
        if (p_IelementsAtEdge2D(2,i) .ne. 0) then
          call nsets_addToDASet (rset,p_IelementsAtEdge2D(2,imt))
        end if
      end do
      
    case (NDIM3D)
      ! Get data from the triangulation
      call storage_getbase_int (rtriangulation%h_IelementsAtEdgeIdx3d,p_IelementsAtEdgeIdx)
      call storage_getbase_int (rtriangulation%h_IelementsAtEdge3d,p_IelementsAtEdge)
      
      ! Loop through the vertices
      do i=1,size(Iedges)
        imt = Iedges(i)
        iel1 = p_IelementsAtEdgeIdx(imt)
        iel2 = p_IelementsAtEdgeIdx(imt+1)-1
        
        ! Add the adjacent elements to the set
        call nsets_putElements (rset,p_IelementsAtEdge(iel1:iel2))
      end do
    end select

  end subroutine

!************************************************************************

!<subroutine>

  subroutine mshreg_auxMarkElementsAtFaces(rtriangulation,Ifaces,rset)

!<description>
  ! Marks all elements adjacent to the faces in Ifaces.
!</description>

!<input>
  ! Underlying triangulation.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of faces
  integer, dimension(:), intent(in) :: Ifaces
!</input>

!<inputoutput>
  ! Set where all elements adjacent to the faces are marked.
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IelementsAtFace
    integer :: i,ifc

    ! Get data from the triangulation
    call storage_getbase_int2d (rtriangulation%h_IelementsAtFace,p_IelementsAtFace)
    
    ! Add adjacent elements
    do i=1,size(Ifaces)
      ifc = Ifaces(i)

      ! The first element is always there, the 2nd element not necessarily
      call nsets_addToDASet (rset,p_IelementsAtFace(1,ifc))
      
      if (p_IelementsAtFace(2,i) .ne. 0) then
        call nsets_addToDASet (rset,p_IelementsAtFace(2,ifc))
      end if
    end do

  end subroutine

!************************************************************************
! to be implemented... does not work until now.

!# 10.) mshreg_getAssocElements
!#      -> Creates a mesh region that contains all elements associated to
!#         the vertices/edges/faces in a "source" mesh region
!#

!!<subroutine>
!
!  subroutine mshreg_getAssocElements(rmeshRegionSrc, rmeshRegionDst, cidxCalc)
!
!!<description>
!  ! Creates a new mesh region rmesgRegionDst that contains all elements
!  ! associated to the vertices/edges/faces in the mesh region
!  ! rmeshRegionSrc. cidxCalc defines which parts of rmeshRegionSrc
!  ! are used for determining the elements.
!!</description>
!
!!<remarks>
!  ! A value MSHREG_IDX_ELEMENT in cidxCalc is currently ignored. It would have 
!  ! to be defined whether the adjacency of an element to its neighbours is via 
!  ! vertices, edges or faces.
!!</remarks>
!
!!<input>
!  ! A combination of MSHREG_IDX_XXXX constants defined above which
!  ! specifies which index arrays in rmeshRegionSrc are used to determine
!  ! the elements. A value MSHREG_IDX_ALL, e.g., computes the associated
!  ! to vertices, edges and faces. A value MSHREG_IDX_EDGE computes
!  ! only the elements associated to the edges in rmesgRegionSrc, etc.
!  !
!  ! If not specified, MSHREG_IDX_ALL is assumed.
!  integer(I32), intent(in) :: cidxCalc
!
!  ! "Source" mesh region that defines vertices/edges/faces to be analysed.
!  type(t_meshRegion), intent(in) :: rmeshRegionSrc
!!</input>
!
!!<output>
!  ! The mesh region to be created based on rmeshRegionSrc.
!  ! Receives a sorted list of all elements adjacent to vertices/edges/faces
!  ! in rmesgRegionSrc, depending on cidxCalc.
!  type(t_meshRegion), intent(out) :: rmeshRegionDst
!!</output>
!
!!</subroutine>
!
!  ! Some local variables
!  type(t_directAccessIntSet) :: rset
!  integer :: nelInSet
!  integer, dimension(:), pointer :: p_Iarray
!  integer :: idim
!
!    ! Get the dimension of the triangulation
!    idim = rmeshRegionSrc%p_rtriangulation%ndim
!  
!    ! To prevent elements from being collected twice, set up a flag array
!    ! that mark all elements we have.
!    call nsets_initDASet (rset,rmeshRegionSrc%p_rtriangulation%NEL)
!    
!    ! Mark elements adjacent to vertices
!    if ((iand(cidxCalc, MSHREG_IDX_VERTEX) .ne. 0) .and. &
!        (rmeshRegionSrc%h_IvertexIdx .ne. ST_NOHANDLE)) then
!
!      call storage_getbase_int (rmeshRegionSrc%h_IvertexIdx,p_Iarray)
!      call mshreg_auxMarkElementsAtVerts(rmeshRegionSrc%p_rtriangulation,p_Iarray,rset)
!
!    end if
!    
!    ! Mark elements adjacent to edges
!    if ((iand(cidxCalc, MSHREG_IDX_EDGE) .ne. 0) .and. &
!        (idim .ge. NDIM2D) .and. (rmeshRegionSrc%h_IedgeIdx .ne. ST_NOHANDLE)) then
!
!      call storage_getbase_int (rmeshRegionSrc%h_IedgeIdx,p_Iarray)
!      call mshreg_auxMarkElementsAtEdges(rmeshRegionSrc%p_rtriangulation,p_Iarray,rset)
!
!    end if
!    
!    ! Mark elements adjacent to faces
!    if ((iand(cidxCalc, MSHREG_IDX_FACE) .ne. 0) .and. &
!        (idim .ge. NDIM3D) .and. (rmeshRegionSrc%h_IfaceIdx .ne. ST_NOHANDLE)) then
!
!      call storage_getbase_int (rmeshRegionSrc%h_IfaceIdx,p_Iarray)
!      call mshreg_auxMarkElementsAtFaces(rmeshRegionSrc%p_rtriangulation,p_Iarray,rset)
!
!    end if
!    
!    ! Form a list of elements from the set
!    nelInSet = nsets_nelementsInSet (rset)
!
!    call storage_new("mshreg_getAssocElements", "p_IidxArray", nelInSet, &
!        ST_INT, rmeshRegionDst%h_IelementIdx, ST_NEWBLOCK_NOINIT)
!    call storage_getbase_int(rmeshRegionDst%h_IelementIdx, p_Iarray)
!
!    call nsets_getElements (rset,p_Iarray)
!        
!    ! Release the element flags
!    call nsets_doneDASet (rset)
!
!  end subroutine

end module meshregion

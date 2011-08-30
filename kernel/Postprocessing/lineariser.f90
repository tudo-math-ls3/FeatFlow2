!##############################################################################
!# ****************************************************************************
!# <name> lineariser </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains various routines to convert an arbitrary
!# finite element vector into its linearised representation which can
!# be used for visualuzation. Given a FE-solution and its underlying
!# triangulation a new triangulation is generated which consists of
!# simplex elements (lines in 1D, triangles in 2D, tetrahedras in
!# 3D). A piecewise linear representation of the FE-solution is
!# generated, whereby the new triangulation is refined if required.
!#
!# The module contains the following subroutines:
!#
!# 1.) lin_lineariseVectorGlobal
!#     -> Linearises an FE-vector by performing global refinement
!#        of the underlying triangulation structure
!#
!# 2.) lin_refineTriaGlobal
!#     -> Refines triangulation structure globally
!#
!# 3.) lin_evalAtDOFScalar
!#     -> Evaluate scalar FE-vector at degrees of freedom
!#
!#
!# The module contains the following auxiliary subroutines:
!#
!# 1.) lin_prepareUnshareVertices
!#      -> Compute auxiliary quantities as preparation for unsharing
!#        unsharing vertices which belong to more than one element
!#
!# 2.) lin_unshareVertices
!#     -> Duplicate vertices which shared by more than one element
!#
!# 3.) lin_convertQuad2Tri
!#     -> Convert all quadrilaterals into two triangles
!#
!# </purpose>
!##############################################################################

module lineariser

  
  use basicgeometry
  use boundary
  use derivatives
  use element
  use feevaluation
  use fsystem
  use genoutput
  use geometryaux
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private

  public :: lin_lineariseVectorGlobal
  public :: lin_refineTriaGlobal
  public :: lin_evalAtDOFScalar

  interface lin_lineariseVectorGlobal
    module procedure lin_lineariseVecScalarGlobal
    module procedure lin_lineariseVecBlockGlobal
  end interface
  

!<constants>
!<constantblock description="Linearisation Types">
! Simple evaluation of DOF`s
  integer(I32), parameter, public :: LIN_TYPE_DOFEVAL = 0

  ! L2-projection
  integer(I32), parameter, public :: LIN_TYPE_L2_PROJ = 1
!</constantblock>
!</constants>

contains

!<subroutine>
  subroutine lin_lineariseVecScalarGlobal(rvectorSrc, rspatialDiscrDest,&
      rtriangulationDest, rvectorDest, nrefsteps, clintype, bdiscontinuous)

!<description>
    ! This subroutine converts the source FE-vector rvectorSrc into a
    ! finite element vector revectorDest which uses simplex elements.
    ! The underlying triangulation is refined globally.
!</description>

!<input>
    ! Source FE-vector
    type(t_vectorScalar), intent(in) :: rvectorSrc

    ! Number of global refinement steps
    ! If nrefsteps=0 then no additional refinement is performed.
    integer, intent(in) :: nrefsteps

    ! Type of linearisation
    integer, intent(in) :: clintype

    ! OPTIONAL: Allow discontinuities at element boundaries
    ! If not present, then continuity is assumed at element boundaries
    logical, intent(in), optional :: bdiscontinuous
!</input>

!<output>
    ! Destination spatial discretisation structure
    type(t_spatialDiscretisation), intent(out) :: rspatialDiscrDest

    ! Destination triangulation structure
    type(t_triangulation), intent(out) :: rtriangulationDest

    ! Destination FE-vector
    type(t_vectorScalar), intent(out) :: rvectorDest
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrSrc
    type(t_triangulation), pointer :: p_rtriangulationSrc
    integer, dimension(:), pointer :: p_ImacroElements

    
    ! Get discretization of source vector
    if (associated(rvectorSrc%p_rspatialDiscr)) then
      p_rspatialDiscrSrc => rvectorSrc%p_rspatialDiscr
    else
      call output_line('Source vector must provide spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecScalarGlobal')
      call sys_halt()
    end if

    ! Get triangulation of source vector
    if (associated(p_rspatialDiscrSrc%p_rtriangulation)) then
      p_rtriangulationSrc => p_rspatialDiscrSrc%p_rtriangulation
    else
      call output_line('Source vector must provide triangulation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecScalarGlobal')
      call sys_halt()
    end if

    ! Generate simplex mesh
    call lin_refineTriaGlobal(p_rtriangulationSrc, rtriangulationDest,&
        nrefsteps, p_ImacroElements, bdiscontinuous)
    
    ! Generate simplex discretisation
    select case(rtriangulationDest%ndim)
    case (NDIM1D)
      call spdiscr_initDiscr_simple(rspatialDiscrDest, EL_P1_1D,&
          SPDISC_CUB_AUTOMATIC, rtriangulationDest)
    case (NDIM2D)
      call spdiscr_initDiscr_simple(rspatialDiscrDest, EL_P1_2D,&
          SPDISC_CUB_AUTOMATIC, rtriangulationDest)
    case (NDIM3D)
      call spdiscr_initDiscr_simple(rspatialDiscrDest, EL_P1_3D,&
          SPDISC_CUB_AUTOMATIC, rtriangulationDest)
    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecScalarGlobal')
      call sys_halt()
    end select
    
    ! Create scalar vector
    call lsyssc_createVector(rspatialDiscrDest, rvectorDest,&
        .false., rvectorSrc%cdataType)
    
    ! Linearise source vector
    select case(clintype)
    case (LIN_TYPE_DOFEVAL)
      call lin_evalAtDOFScalar(rvectorSrc, rspatialDiscrDest,&
          rtriangulationDest, rvectorDest, p_ImacroElements)

    case default
      call output_line('Unsupported type of linearisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecScalarGlobal')
      call sys_halt()
    end select
    
    ! Deallocate temporal memory
    deallocate(p_ImacroElements)

  end subroutine lin_lineariseVecScalarGlobal

  !*****************************************************************************

!<subroutine>

  subroutine lin_lineariseVecBlockGlobal(rvectorSrc, rblockDiscrDest,&
      rtriangulationDest, rvectorDest, nrefsteps, clintype, bdiscontinuous)

!<description>
    ! This subroutine converts the source FE-vector rvectorSrc into a
    ! finite element vector revectorDest which uses simplex elements.
    ! The underlying triangulation is refined globally.
!</description>

!<input>
    ! Source FE-vector
    type(t_vectorBlock), intent(in) :: rvectorSrc

    ! Number of global refinement steps
    ! If nrefsteps=0 then no additional refinement is performed.
    integer, intent(in) :: nrefsteps

    ! Type of linearisation
    integer, intent(in) :: clintype

    ! OPTIONAL: Allow discontinuities at element boundaries
    ! If not present, then continuity is assumed at element boundaries
    logical, intent(in), optional :: bdiscontinuous
!</input>

!<output>
    ! Destination block discretisation structure
    type(t_blockDiscretisation), intent(out) :: rblockDiscrDest

    ! Destination triangulation structure
    type(t_triangulation), intent(out) :: rtriangulationDest

    ! Destination FE-vector
    type(t_vectorBlock), intent(out) :: rvectorDest
!</output>
!</subroutine>

    ! local variables
    type(t_blockDiscretisation), pointer :: p_rblockDiscrSrc
    type(t_triangulation), pointer :: p_rtriangulationSrc
    integer, dimension(:), pointer :: p_ImacroElements
    integer :: iblock

    ! Get block discretization of source vector
    if (associated(rvectorSrc%p_rblockDiscr)) then
      p_rblockDiscrSrc => rvectorSrc%p_rblockDiscr
    else
      call output_line('Source vector must provide block discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecBlockGlobal')
      call sys_halt()
    end if

    ! Get triangulation of source vector
    if (associated(p_rblockDiscrSrc%p_rtriangulation)) then
      p_rtriangulationSrc => p_rblockDiscrSrc%p_rtriangulation
    else
      call output_line('Source vector must provide triangulation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecBlockGlobal')
      call sys_halt()
    end if

    ! Generate simplex mesh
    call lin_refineTriaGlobal(p_rtriangulationSrc, rtriangulationDest,&
        nrefsteps, p_ImacroElements, bdiscontinuous)

    ! Generate simplex block discretisation
    call spdiscr_initBlockDiscr(rblockDiscrDest,&
        p_rblockDiscrSrc%ncomponents, p_rtriangulationSrc)
    select case(rtriangulationDest%ndim)
    case (NDIM1D)
      call spdiscr_initDiscr_simple(rblockDiscrDest%RspatialDiscr(1),&
          EL_P1_1D, SPDISC_CUB_AUTOMATIC, rtriangulationDest)
    case (NDIM2D)
      call spdiscr_initDiscr_simple(rblockDiscrDest%RspatialDiscr(1),&
          EL_P1_2D, SPDISC_CUB_AUTOMATIC, rtriangulationDest)
    case (NDIM3D)
      call spdiscr_initDiscr_simple(rblockDiscrDest%RspatialDiscr(1),&
          EL_P1_3D, SPDISC_CUB_AUTOMATIC, rtriangulationDest)
    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecBlockGlobal')
      call sys_halt()
    end select

    do iblock = 2, p_rblockDiscrSrc%ncomponents
      call spdiscr_duplicateDiscrSc(rblockDiscrDest%RspatialDiscr(1),&
          rblockDiscrDest%RspatialDiscr(iblock), .true.)
    end do

    ! Create block vector
    call lsysbl_createVectorBlock(rblockDiscrDest, rvectorDest,&
        .false., rvectorSrc%cdataType)

    ! Linearise source vector
    select case(clintype)
    case (LIN_TYPE_DOFEVAL)
      do iblock = 1, p_rblockDiscrSrc%ncomponents
        call lin_evalAtDOFScalar(rvectorSrc%RvectorBlock(iblock),&
            rblockDiscrDest%RspatialDiscr(iblock), rtriangulationDest,&
            rvectorDest%RvectorBlock(iblock), p_ImacroElements)
      end do
      
    case default
      call output_line('Unsupported type of linearisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lin_lineariseVecBlockGlobal')
      call sys_halt()
    end select
    
    ! Deallocate temporal memory
    deallocate(p_ImacroElements)
    
  end subroutine lin_lineariseVecBlockGlobal

  !*****************************************************************************

!<subroutine>

  subroutine lin_refineTriaGlobal(rtriangulationSrc, rtriangulationDest,&
      nrefsteps, p_ImacroElements, bdiscontinuous)

!<description>
    ! This subroutine refines the source triangulation rtriaSrc
    ! nrefsteps time and stores the resulting triangulation in rtriaDest
!</description>

!<input>
    ! Source triangulation
    type(t_triangulation), intent(in) :: rtriangulationSrc

    ! Number of global refinement steps
    ! If nrefsteps=0 then no additional refinement is performed.
    integer, intent(in), optional :: nrefsteps

    ! OPTIONAL: Allow discontinuities at element boundaries
    ! If not present, then continuity is assumed at element boundaries
    logical, intent(in), optional :: bdiscontinuous
!</input>

!<output>
    ! Destination triangulation
    type(t_triangulation), intent(out) :: rtriangulationDest

    ! OPTIONAL: Mapping of element numbers to macro elements
    ! If present, this pointer will by allocated with correct size
    ! and the mapping will be returned by this subroutine
    integer, dimension(:), intent(out), pointer, optional :: p_ImacroElements
!</output>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoordsDest,p_DvertexCoordsSrc
    integer, dimension(:,:), pointer :: p_IverticesAtElementDest,p_IverticesAtElementSrc
    integer, dimension(:,:), pointer :: p_IverticesAtElementTmp
    integer, dimension(:,:), pointer :: p_IneighboursAtElementDest,p_IneighboursAtElementSrc
    integer, dimension(:), pointer :: p_InodalPropertyDest,p_InodalPropertySrc
    integer, dimension(:), pointer :: p_IverticesAtBoundaryDest,p_IverticesAtBoundarySrc
    integer, dimension(:), pointer :: p_ImacroElementsTmp
    integer, dimension(2) :: Isize
    integer :: h_IverticesAtElement
    integer :: iel,istep,ive,nel,nvt
    logical :: bisDiscontinuous

    bisDiscontinuous = .false.
    if (present(bdiscontinuous)) bisDiscontinuous = bdiscontinuous
    
    ! What spatial dimension are we?
    select case(rtriangulationSrc%ndim)
      
    case(NDIM1D)
      ! Set pointers
      call storage_getbase_double2d(&
          rtriangulationSrc%h_DvertexCoords, p_DvertexCoordsSrc)
      call storage_getbase_int2d(&
          rtriangulationSrc%h_IverticesAtElement, p_IverticesAtElementSrc)
      call storage_getbase_int(&
          rtriangulationSrc%h_InodalProperty, p_InodalPropertySrc)
      
      ! Most attributes are initialised by INTENT(OUT) with standard values
      rtriangulationDest%ndim = NDIM1D
      rtriangulationDest%NNVE = 2
      rtriangulationDest%NBCT = rtriangulationSrc%NBCT
      rtriangulationDest%NVBD = rtriangulationSrc%NVBD

      ! The number of elements/vertices will be updated below
      rtriangulationDest%NVT  = rtriangulationSrc%NVT
      rtriangulationDest%NEL  = rtriangulationSrc%NEL
      
      ! Duplicate shared vertices? 
      if (bisDiscontinuous)&
          call lin_prepareUnshareVertices(p_IverticesAtElementSrc,&
          p_InodalPropertySrc, rtriangulationDest%NVT,&
          rtriangulationDest%NVBD, rtriangulationDest%NMT,&
          rtriangulationDest%NEL)

      ! Compute number of vertices/elements after global refinement
      do istep = 1, nrefsteps
        ! New vertices are inserted at element midpoints and
        ! each element is subdivided into two new elements
        rtriangulationDest%NVT = rtriangulationDest%NVT + rtriangulationDest%NEL
        rtriangulationDest%NEL = 2*rtriangulationDest%NEL
      end do

      ! Destination triangulation has only line segments
      rtriangulationDest%InelOfType(TRIA_NVELINE1D) = rtriangulationDest%NEL

      ! Allocate memory for coordinate vector
      Isize = (/1, rtriangulationDest%NVT/)
      call storage_new ('lin_refineTriaGlobal', 'DCORVG', Isize,&
          ST_DOUBLE, rtriangulationDest%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2d(&
          rtriangulationDest%h_DvertexCoords, p_DvertexCoordsDest)
      
      ! Allocate memory for vertices at element list
      Isize = (/2, rtriangulationDest%NEL/)
      call storage_new ('lin_refineTriaGlobal', 'KVERT', Isize,&
          ST_INT, rtriangulationDest%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2d(&
          rtriangulationDest%h_IverticesAtElement, p_IverticesAtElementDest)
      
      ! Allocate memory for nodal property array
      call storage_new ('lin_refineTriaGlobal', 'KNPR', rtriangulationDest%NVT,&
          ST_INT, rtriangulationDest%h_InodalProperty, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(&
          rtriangulationDest%h_InodalProperty, p_InodalPropertyDest)
      
      ! Copy arrays which do not change
      call storage_copy(rtriangulationSrc%h_IboundaryCpIdx,&
          rtriangulationDest%h_IboundaryCpIdx)
      call storage_copy(rtriangulationSrc%h_IverticesAtBoundary,&
          rtriangulationDest%h_IverticesAtBoundary)

      ! Copy data from source to destination triangulation
      call lalg_copyVector(&
          p_InodalPropertySrc, p_InodalPropertyDest, rtriangulationSrc%NVT)
      call lalg_copyVector(&
          p_DvertexCoordsSrc, p_DvertexCoordsDest, 1, rtriangulationSrc%NVT)
      call lalg_copyVector(&
          p_IverticesAtElementSrc, p_IverticesAtElementDest, 2, rtriangulationSrc%NEL)
      
      ! Make a backup of some quantities
      nel = rtriangulationSrc%NEL
      nvt = rtriangulationSrc%NVT

      ! Duplicate shared vertices?
        if (bisDiscontinuous)&
            call lin_unshareVertices(p_DvertexCoordsDest,&
            p_IverticesAtElementDest, p_InodalPropertyDest, nvt, nel)
      
      ! Refine mesh globally
      if (present(p_ImacroElements)) then
        allocate(p_ImacroElements(rtriangulationDest%NEL))
        call doRefine1D(p_DvertexCoordsDest, p_IverticesAtElementDest,&
            p_InodalPropertyDest, p_ImacroElements, nvt, nel, .true., nrefsteps)
      else
        allocate(p_ImacroElementsTmp(rtriangulationDest%NEL))
        call doRefine1D(p_DvertexCoordsDest, p_IverticesAtElementDest,&
            p_InodalPropertyDest, p_ImacroElementsTmp, nvt, nel, .true., nrefsteps)
        deallocate(p_ImacroElementsTmp)
      end if
      

    case (NDIM2D)
      ! Set pointers
      call storage_getbase_double2d(&
          rtriangulationSrc%h_DvertexCoords, p_DvertexCoordsSrc)
      call storage_getbase_int2d(&
          rtriangulationSrc%h_IverticesAtElement, p_IverticesAtElementSrc)
      call storage_getbase_int2d(&
          rtriangulationSrc%h_IneighboursAtElement, p_IneighboursAtElementSrc)
      call storage_getbase_int(&
          rtriangulationSrc%h_InodalProperty, p_InodalPropertySrc)
      call storage_getbase_int(&
          rtriangulationSrc%h_IverticesAtBoundary, p_IverticesAtBoundarySrc)
     
      ! Most attributes are initialised by INTENT(OUT) with standard values.
      rtriangulationDest%ndim = NDIM2D
      rtriangulationDest%NNVE = 3
      rtriangulationDest%NBCT = rtriangulationSrc%NBCT
      
      ! The number of elements/vertices/edges will be updated below
      rtriangulationDest%NVBD = rtriangulationSrc%NVBD
      rtriangulationDest%NVT  = rtriangulationSrc%NVT
      rtriangulationDest%NEL  = rtriangulationSrc%NEL
      rtriangulationDest%NMT  = rtriangulationSrc%NMT
      
      ! Duplicate shared vertices?
      if (bisDiscontinuous)&
          call lin_prepareUnshareVertices(p_IverticesAtElementSrc,&
          p_InodalPropertySrc, rtriangulationDest%NVT,&
          rtriangulationDest%NVBD, rtriangulationDest%NMT,&
          rtriangulationDest%NEL)

      ! If the source triangulation consists only of triangles, then
      ! we can simply copy and refine it without further checks
      if ((rtriangulationSrc%NNVE .eq. 3) .or.&
          (rtriangulationSrc%NEL  .eq. &
          rtriangulationSrc%InelOfType(TRIA_NVETRI2D))) then
        
        ! Compute number of vertices/elements after global refinement
        do istep = 1, nrefsteps
          ! New vertices are inserted at edge midpoints
          rtriangulationDest%NVT = rtriangulationDest%NVT + rtriangulationDest%NMT
          
          ! Each edge is subdivided into two edges (also on the boundary)
          ! Moreover, each element gives rise to three new internal edges
          rtriangulationDest%NMT  = 2*rtriangulationDest%NMT+3*rtriangulationDest%NEL
          rtriangulationDest%NVBD = 2*rtriangulationDest%NVBD
          
          ! Each element is subdivided into four elements
          rtriangulationDest%NEL = 4*rtriangulationDest%NEL
        end do

        ! Destination triangulation has only triangles
        rtriangulationDest%InelOfType(TRIA_NVETRI2D) = rtriangulationDest%NEL

        ! Allocate memory for coordinate vector
        Isize = (/2, rtriangulationDest%NVT/)
        call storage_new ('lin_refineTriaGlobal', 'DCORVG', Isize,&
            ST_DOUBLE, rtriangulationDest%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2d(&
            rtriangulationDest%h_DvertexCoords, p_DvertexCoordsDest)

        ! Allocate memory for vertices at element list
        Isize = (/3, rtriangulationDest%NEL/)
        call storage_new ('lin_refineTriaGlobal', 'KVERT', Isize,&
            ST_INT, rtriangulationDest%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2d(&
            rtriangulationDest%h_IverticesAtElement, p_IverticesAtElementDest)

        ! Allocate memory for neighbours at element list
        Isize = (/3, rtriangulationDest%NEL/)
        call storage_new ('lin_refineTriaGlobal', 'KADJ', Isize,&
            ST_INT, rtriangulationDest%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2d(&
            rtriangulationDest%h_IneighboursAtElement, p_IneighboursAtElementDest)

        ! Allocate memory for nodal property array
        call storage_new ('lin_refineTriaGlobal', 'KNPR', rtriangulationDest%NVT,&
            ST_INT, rtriangulationDest%h_InodalProperty, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(&
            rtriangulationDest%h_InodalProperty, p_InodalPropertyDest)

        ! Allocate memory for vertices at boundary list
        call storage_new ('lin_refineTriaGlobal', 'KVBD', rtriangulationDest%NVBD,&
            ST_INT, rtriangulationDest%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(&
            rtriangulationDest%h_IverticesAtBoundary, p_IverticesAtBoundaryDest)
                
        ! Copy arrays which do not change its size
        call storage_copy(rtriangulationSrc%h_IboundaryCpIdx,&
            rtriangulationDest%h_IboundaryCpIdx)
        
        ! Copy data from source to destination triangulation
        call lalg_copyVector(&
            p_InodalPropertySrc, p_InodalPropertyDest, rtriangulationSrc%NVT)
        call lalg_copyVector(&
            p_DvertexCoordsSrc, p_DvertexCoordsDest, 2, rtriangulationSrc%NVT)

        ! If the source triangulation results from dividing some
        ! quadrilaterals into triangles then we have to copy the
        ! list of vertices at each element by hand
        if (rtriangulationSrc%NNVE .eq. 3) then
          call lalg_copyVector(&
              p_IverticesAtElementSrc, p_IverticesAtElementDest,&
              3, rtriangulationSrc%NEL)
          if (bisDiscontinuous) then
            call lalg_clearVector(p_IneighboursAtElementDest)
          else
            call lalg_copyVector(&
                p_IneighboursAtElementSrc, p_IneighboursAtElementDest,&
                3, rtriangulationSrc%NEL)
          end if
        else
          ! Copy arrays by hand
          do iel = 1, rtriangulationDest%NEL
            p_IverticesAtElementDest(1:3,iel) = p_IverticesAtElementSrc(1:3,iel)
          end do

          if (bisDiscontinuous) then
            call lalg_clearVector(p_IneighboursAtElementDest)
          else
            do iel = 1, rtriangulationDest%NEL
              p_IneighboursAtElementDest(1:3,iel) = p_IneighboursAtElementSrc(1:3,iel)
            end do
          end if
        end if
        
        ! Make a backup of some quantities
        nel = rtriangulationSrc%NEL
        nvt = rtriangulationSrc%NVT
        
        ! Duplicate shared vertices?
        if (bisDiscontinuous)&
            call lin_unshareVertices(p_DvertexCoordsDest,&
            p_IverticesAtElementDest, p_InodalPropertyDest, nvt, nel)
        
        ! Refine mesh globally
        if (present(p_ImacroElements)) then
          allocate(p_ImacroElements(rtriangulationDest%NEL))
          call doRefine2D(p_DvertexCoordsDest, p_IverticesAtElementDest,&
              p_IneighboursAtElementDest, p_InodalPropertyDest,&
              p_ImacroElements, nvt, nel, .true., nrefsteps)
        else
          allocate(p_ImacroElementsTmp(rtriangulationDest%NEL))
          call doRefine2D(p_DvertexCoordsDest, p_IverticesAtElementDest,&
              p_IneighboursAtElementDest, p_InodalPropertyDest,&
              p_ImacroElementsTmp, nvt, nel, .true., nrefsteps)
          deallocate(p_ImacroElementsTmp)
        end if
        
      else
        
        ! Triangles are kept as is and quadrilaterals are subdivided
        ! into two triangles without introducing new vertices.
        ! Hence, each quadrilateral gives rise to one new edge
        rtriangulationDest%NMT = rtriangulationDest%NMT +&
                                 rtriangulationSrc%InelOfType(TRIA_NVEQUAD2D)
        rtriangulationDest%NEL = rtriangulationSrc%InelOfType(TRIA_NVETRI2D) +&
                                 rtriangulationSrc%InelOfType(TRIA_NVEQUAD2D)*2
        
        ! Compute number of vertices/elements after global refinement
        do istep = 1, nrefsteps
          ! New vertices are inserted at edge midpoints
          rtriangulationDest%NVT = rtriangulationDest%NVT + rtriangulationDest%NMT
          
          ! Each edge is subdivided into two edges (also on the boundary)
          ! Moreover, each element gives rise to three new internal edges
          rtriangulationDest%NMT  = 2*rtriangulationDest%NMT+3*rtriangulationDest%NEL
          rtriangulationDest%NVBD = 2*rtriangulationDest%NVBD
          
          ! Each element is subdivided into four elements
          rtriangulationDest%NEL = 4*rtriangulationDest%NEL
        end do

        ! Destination triangulation has only triangles
        rtriangulationDest%InelOfType(TRIA_NVETRI2D) = rtriangulationDest%NEL
        
        ! Allocate memory for coordinate vector
        Isize = (/2, rtriangulationDest%NVT/)
        call storage_new ('lin_refineTriaGlobal', 'DCORVG', Isize,&
            ST_DOUBLE, rtriangulationDest%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2d(&
            rtriangulationDest%h_DvertexCoords, p_DvertexCoordsDest)
        
        ! Allocate auxiliary memory for vertices at element list
        Isize = (/rtriangulationSrc%NNVE, rtriangulationDest%NEL/)
        h_IverticesAtElement = ST_NOHANDLE
        call storage_new ('lin_refineTriaGlobal', 'KVERT', Isize,&
            ST_INT, h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2D(h_IverticesAtElement, p_IverticesAtElementTmp)
        
        ! Allocate memory for vertices at element list
        Isize = (/3, rtriangulationDest%NEL/)
        call storage_new ('lin_refineTriaGlobal', 'KVERT', Isize,&
            ST_INT, rtriangulationDest%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2D(rtriangulationDest%h_IverticesAtElement,&
            p_IverticesAtElementDest)

        ! Allocate memory for neighbours at element list
        Isize = (/rtriangulationSrc%NNVE, rtriangulationDest%NEL/)
        call storage_new ('lin_refineTriaGlobal', 'KADJ', Isize,&
            ST_INT, rtriangulationDest%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2d(&
            rtriangulationDest%h_IneighboursAtElement, p_IneighboursAtElementDest)
        
        ! Allocate memory for nodal property array
        call storage_new ('lin_refineTriaGlobal', 'KNPR', rtriangulationDest%NVT,&
            ST_INT, rtriangulationDest%h_InodalProperty, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(&
            rtriangulationDest%h_InodalProperty, p_InodalPropertyDest)
        
        ! Allocate memory for vertices at boundary list
        call storage_new ('lin_refineTriaGlobal', 'KVBD', rtriangulationDest%NVBD,&
            ST_INT, rtriangulationDest%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(&
            rtriangulationDest%h_IverticesAtBoundary, p_IverticesAtBoundaryDest)
        
        ! Copy arrays which do not change its size
        call storage_copy(rtriangulationSrc%h_IboundaryCpIdx,&
            rtriangulationDest%h_IboundaryCpIdx)
        
        ! Copy data from source to destination triangulation and to
        ! the auxiliary list of vertices at elements
        call lalg_copyVector(&
            p_InodalPropertySrc, p_InodalPropertyDest, rtriangulationSrc%NVT)
        call lalg_copyVector(&
            p_DvertexCoordsSrc, p_DvertexCoordsDest, 2, rtriangulationSrc%NVT)
        call lalg_copyVector(&
            p_IverticesAtElementSrc, p_IverticesAtElementTmp,&
            rtriangulationSrc%NNVE, rtriangulationSrc%NEL)
        
        if (bisDiscontinuous) then
          call lalg_clearVector(p_IneighboursAtElementDest)
        else
          call lalg_copyVector(&
              p_IneighboursAtElementSrc, p_IneighboursAtElementDest,&
              rtriangulationSrc%NNVE, rtriangulationSrc%NEL)
        end if
        
        ! Make a backup of some quantities
        nel = rtriangulationSrc%NEL
        nvt = rtriangulationSrc%NVT
        
        ! Duplicate shared vertices?
        if (bisDiscontinuous)&
            call lin_unshareVertices(p_DvertexCoordsDest,&
            p_IverticesAtElementTmp, p_InodalPropertyDest, nvt, nel)
        
        ! Convert quadrilaterals into triangles and perform global refinement
        if (present(p_ImacroElements)) then
          allocate(p_ImacroElements(rtriangulationDest%NEL))
          call lin_convertQuad2Tri(p_DvertexCoordsDest,&
              p_IverticesAtElementTmp, p_IverticesAtElementDest,&
              p_IneighboursAtElementDest, p_ImacroElements, nel, .true.)
          call doRefine2D(p_DvertexCoordsDest, p_IverticesAtElementDest,&
              p_IneighboursAtElementDest, p_InodalPropertyDest,&
              p_ImacroElements, nvt, nel, .false., nrefsteps)
        else
          allocate(p_ImacroElementsTmp(rtriangulationDest%NEL))
          call lin_convertQuad2Tri(p_DvertexCoordsDest,&
              p_IverticesAtElementTmp, p_IverticesAtElementDest,&
              p_IneighboursAtElementDest, p_ImacroElementsTmp, nel, .true.)
          call doRefine2D(p_DvertexCoordsDest, p_IverticesAtElementDest,&
              p_IneighboursAtElementDest, p_InodalPropertyDest,&
              p_ImacroElementsTmp, nvt, nel, .false., nrefsteps)
          deallocate(p_ImacroElementsTmp)
        end if
        
        ! Free temporal memory
        call storage_free(h_IverticesAtElement)
        call storage_free(rtriangulationDest%h_IneighboursAtElement)
                
      end if
      
    end select

    ! Regenerate boundary information

  contains

    ! Here, some working routines follow
    
    !**************************************************************
    ! Global refinement of triangulation in 1D
    
    subroutine doRefine1D(DvertexCoords, IverticesAtElement, InodalProperty,&
        ImacroElements, NVT, NEL, binit, nrefsteps)

      ! On input: vertex coordinates of source triangulation
      ! On output: vertex coordinates of destination triangulation
      ! DIMENSION(NDIM1D,NVT)
      real(DP), dimension(:,:), intent(inout) :: DvertexCoords
      
      ! Vertices at element list
      ! On input: vertices at element list of source triangulation
      ! On output: vertices at element list of destination triangulation
      ! DIMENSION(NVE,NEL)
      integer, dimension(:,:), intent(inout) :: IverticesAtElement

      ! On input: nodal property list of source triangulation
      ! On output: nodal property list of destination triangulation
      ! DIMENSION(NVT)
      integer, dimension(:), intent(inout) :: InodalProperty
      
      ! On output: mapping between macro element in source triangulation
      ! and element in destination triangulation
      integer, dimension(:), intent(inout) :: ImacroElements

      ! On input: number of vertices of source triangulation
      ! On output: number of vertices of destination triangulation
      integer, intent(inout) :: NVT
      
      ! On input: number of elements of source triangulation
      ! On output: number of elements of destination triangulation
      integer, intent(inout) :: NEL

      ! Flag to enable initialisation of mapping
      logical, intent(in) :: binit

      ! Number of global refinement steps
      integer, intent(in) :: nrefsteps


      ! local variables
      integer :: iel,istep,NEL0,NNEL,i1,i2

      ! Make a backup of some quantities
      NEL0 = NEL
      NNEL = NEL
      
      ! Initialse element mapping
      if (binit) then
        do iel = 1, NEL0
          ImacroElements(iel) = iel
        end do
      end if

      ! Perform global refinement steps
      do istep = 1, nrefsteps

        ! Loop over all elements of the source triangulation
        do iel = 1, NNEL
          
          ! Get vertex numbers
          i1 = IverticesAtElement(1,iel)
          i2 = IverticesAtElement(2,iel)
          
          ! Increment number of vertices
          NVT = NVT+1
          
          ! Insert new interior vertex at the midpoint
          DvertexCoords(1,NVT) = 0.5_DP*(DvertexCoords(1,i1)+DvertexCoords(1,i2))
          InodalProperty(NVT)  = 0

          ! Increment number of elements
          NEL = NEL+1

          ! Update vertices at element lists
          IverticesAtElement(1,iel) = i1
          IverticesAtElement(2,iel) = NVT
          IverticesAtElement(1,NEL) = NVT
          IverticesAtElement(2,NEL) = i2

          ! Set mapping to macro element
          ImacroElements(NEL) = iel
        end do

        ! Update total number of elements
        NNEL = NEL
      end do

      ! Convert element mapping to original source triangulation
      do iel = NEL0+1, NEL
        ImacroElements(iel) = ImacroElements(ImacroElements(iel))
      end do

    end subroutine doRefine1D
   
    !**************************************************************
    ! Global refinement of triangulation in 2D
    
    subroutine doRefine2D(DvertexCoords, IverticesAtElement, IneighboursAtElement,&
        InodalProperty, ImacroElements, NVT, NEL, binit, nrefsteps)

      ! On input: vertex coordinates of source triangulation
      ! On output: vertex coordinates of destination triangulation
      ! DIMENSION(NDIM2D,NVT)
      real(DP), dimension(:,:), intent(inout) :: DvertexCoords
      
      ! Vertices at element list
      ! On input: vertices at element list of source triangulation
      ! On output: vertices at element list of destination triangulation
      ! DIMENSION(NVE,NEL)
      integer, dimension(:,:), intent(inout) :: IverticesAtElement

      ! Neighbours at element list
      ! On input: neighbours at element list of source triangulation
      ! On output: neighbours at element list of destination triangulation
      ! DIMENSION(NVE,NEL)
      integer, dimension(:,:), intent(inout) :: IneighboursAtElement

      ! On input: nodal property list of source triangulation
      ! On output: nodal property list of destination triangulation
      ! DIMENSION(NVT)
      integer, dimension(:), intent(inout) :: InodalProperty
      
      ! On output: mapping between macro element in source triangulation
      ! and element in destination triangulation
      integer, dimension(:), intent(inout) :: ImacroElements

      ! On input: number of vertices of source triangulation
      ! On output: number of vertices of destination triangulation
      integer, intent(inout) :: NVT
      
      ! On input: number of elements of source triangulation
      ! On output: number of elements of destination triangulation
      integer, intent(inout) :: NEL

      ! Flag to enable initialisation of mapping
      logical, intent(in) :: binit

      ! Number of global refinement steps
      integer, intent(in) :: nrefsteps


      ! local variables
      real(DP), dimension(NDIM2D) :: Dcoords
      integer :: NEL0,NNEL,i1,i2,i3,i4,i5,i6,iel,istep,jel,kel

      
      ! Make a backup of some quantities
      NEL0 = NEL
      NNEL = NEL
      
      ! Initialse element mapping
      if (binit) then
        do iel = 1, NEL0
          ImacroElements(iel) = iel
        end do
      end if

      ! Perform global refinement steps
      do istep = 1, nrefsteps
        
        ! Loop over all elements of the source triangulation
        do iel = 1, NNEL
          
          ! Get vertex numbers
          i1 = IverticesAtElement(1,iel)
          i2 = IverticesAtElement(2,iel)
          i3 = IverticesAtElement(3,iel)
          
          ! New vertices are inserted at the edge midpoints if
          ! (a) the edge is located at the boundary, i.e. there is no
          !     adjacent element neighbour
          ! (b) the edge has not been subdivided before, i.e. this is
          !     the first visit of the edge

          ! Process first edge (I1,I2)
          jel = IneighboursAtElement(1,iel)
          if ((iel .lt. jel) .or. (jel .eq. 0)) then
            ! Insert new vertex           
            NVT = NVT+1; i4 = NVT
            DvertexCoords(:,i4) = 0.5_DP*(DvertexCoords(:,i1)+DvertexCoords(:,i2))
            if (InodalProperty(i1)*InodalProperty(i2) .ne. 0) then
              InodalProperty(i4) = InodalProperty(i1)
            else
              InodalProperty(i4) = 0
            end if
            
            ! Update neighbouring element list along the edge
            IneighboursAtElement(1,NEL+1) = jel
            IneighboursAtElement(3,NEL+2) = jel
          else
            ! Seach vertex in neighbouring element
            Dcoords = 0.5_DP*(DvertexCoords(:,i1)+DvertexCoords(:,i2))
            do ive = 1, 3
              i4 = IverticesAtElement(ive,jel)
              if (all(Dcoords .eq. DvertexCoords(:,i4))) exit
            end do

            ! We have found the midpoint vertex in the adjacent
            ! element. Update the neighbours at element list for both
            ! triangles located at the current edge (I1,I2)
            adj12: do ive = 1, 3
              kel = IneighboursAtElement(ive,jel)
              if (IneighboursAtElement(1,kel) .eq. iel) then
                IneighboursAtElement(1,kel) = NEL+2
                IneighboursAtElement(3,NEL+2) = kel
                cycle adj12
              end if
              if (IneighboursAtElement(3,kel) .eq. iel) then
                IneighboursAtElement(3,kel) = NEL+1
                IneighboursAtElement(1,NEL+1 ) = kel
                cycle adj12
              end if
            end do adj12
          end if

          ! Process second edge (I2,I3)
          jel = IneighboursAtElement(2,iel)
          if ((iel .lt. jel) .or. (jel .eq. 0)) then
            ! Insert new vertex
            NVT = NVT+1; i5 = NVT
            DvertexCoords(:,i5) = 0.5_DP*(DvertexCoords(:,i2)+DvertexCoords(:,i3))
            if (InodalProperty(i2)*InodalProperty(i3) .ne. 0) then
              InodalProperty(i5) = InodalProperty(i2)
            else
              InodalProperty(i5) = 0
            end if
            
            ! Update neighbouring element list along the edge
            IneighboursAtElement(1,NEL+2) = jel
            IneighboursAtElement(3,NEL+3) = jel
          else
            ! Search vertex in neighbouring element
            Dcoords = 0.5_DP*(DvertexCoords(:,i2)+DvertexCoords(:,i3))
            do ive = 1, 3
              i5 = IverticesAtElement(ive,jel)
              if (all(Dcoords .eq. DvertexCoords(:,i5))) exit
            end do
            
            ! We have found the midpoint vertex in the adjacent
            ! element. Update the neighbours at element list for both
            ! triangles located at the current edge (I2,I3)
            adj23: do ive = 1, 3
              kel = IneighboursAtElement(ive,jel)
              if (IneighboursAtElement(1,kel) .eq. iel) then
                IneighboursAtElement(1,kel) = NEL+3
                IneighboursAtElement(3,NEL+3) = kel
                cycle adj23
              end if
              if (IneighboursAtElement(3,kel) .eq. iel) then
                IneighboursAtElement(3,kel) = NEL+2
                IneighboursAtElement(1,NEL+2) = kel
                cycle adj23
              end if
            end do adj23
          end if

          ! Process third edge (I3,I1)
          jel = IneighboursAtElement(3,iel)
          if ((iel .lt. jel) .or. (jel .eq. 0)) then
            ! Insert new vertex
            NVT = NVT+1; i6 = NVT
            DvertexCoords(:,i6) = 0.5_DP*(DvertexCoords(:,i3)+DvertexCoords(:,i1))
            if (InodalProperty(i3)*InodalProperty(i1) .ne. 0) then
              InodalProperty(i6) = InodalProperty(i3)
            else
              InodalProperty(i6) = 0
            end if
            
            ! Update neighbouring element list along the edge
            IneighboursAtElement(1,NEL+3) = IneighboursAtElement(3,iel)
            IneighboursAtElement(3,NEL+1) = IneighboursAtElement(3,iel)
          else
            ! Search vertex in neighbouring element
            Dcoords = 0.5_DP*(DvertexCoords(:,i3)+DvertexCoords(:,i1))
            do ive = 1, 3
              i6 = IverticesAtElement(ive,jel)
              if (all(Dcoords .eq. DvertexCoords(:,i6))) exit
            end do

            ! We have found the midpoint vertex in the adjacent
            ! element. Update the neighbours at element list for both
            ! triangles located at the current edge (I3,I1)
            adj31: do ive = 1, 3
              kel = IneighboursAtElement(ive,jel)
              if (IneighboursAtElement(1,kel) .eq. iel) then
                IneighboursAtElement(1,kel) = NEL+1
                IneighboursAtElement(3,NEL+1) = kel
                cycle adj31
              end if
              if (IneighboursAtElement(3,kel) .eq. iel) then
                IneighboursAtElement(3,kel) = NEL+3
                IneighboursAtElement(1,NEL+3) = kel
                cycle adj31
              end if
            end do adj31
          end if
          
          ! Update vertices at element lists of the refined triangle
          IverticesAtElement(1,iel) = i4
          IverticesAtElement(2,iel) = i5
          IverticesAtElement(3,iel) = i6
          
          IverticesAtElement(1,NEL+1) = i1
          IverticesAtElement(2,NEL+1) = i4
          IverticesAtElement(3,NEL+1) = i6

          IverticesAtElement(1,NEL+2) = i2
          IverticesAtElement(2,NEL+2) = i5
          IverticesAtElement(3,NEL+2) = i4

          IverticesAtElement(1,NEL+3) = i3
          IverticesAtElement(2,NEL+3) = i6
          IverticesAtElement(3,NEL+3) = i5

          ! Set mappings to macro element
          ImacroElements(NEL+1:NEL+3) = iel
          
          ! Update neighbours at element list of the newly created
          ! elements and update list of the interior element
          IneighboursAtElement(2,NEL+1) = iel
          IneighboursAtElement(2,NEL+2) = iel
          IneighboursAtElement(2,NEL+3) = iel
          
          IneighboursAtElement(1,iel) = NEL+2
          IneighboursAtElement(2,iel) = NEL+3
          IneighboursAtElement(3,iel) = NEL+1

          ! Increment number of elements
          NEL = NEL+3
        end do
        
        ! Update total number of elements
        NNEL = NEL
      end do
      
      ! Convert element mapping to original source triangulation
      do iel = NEL0+1, NEL
        ImacroElements(iel) = ImacroElements(ImacroElements(iel))
      end do
    end subroutine doRefine2D

  end subroutine lin_refineTriaGlobal

  !*****************************************************************************

!<subroutine>

  subroutine lin_evalAtDOFScalar(rvectorSrc, rspatialDiscrDest,&
      rtriangulationDest, rvectorDest, ImacroElements)


!<description>
    ! This subroutine linearises the given FE-solution rvectorSrc 
    ! by evaluating it at the degrees of freedom of the destination
    ! discretisation structure rspatialDiscreDest and stores the
    ! values in the destination FE-vector rvectorDest
!</description>

!<input>
    ! Source FE-vector
    type(t_vectorScalar), intent(in) :: rvectorSrc

    ! Destination spatial discretisation structure
    type(t_spatialDiscretisation), intent(in) :: rspatialDiscrDest

    ! Destination triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationDest

    ! Mapping of element numbers to the macro element
    integer, dimension(:), intent(in) :: ImacroElements
!</input>

!<inputoutput>
    ! Destination FE-vector
    type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:,:), allocatable :: Dpoints
    real(DP), dimension(:), allocatable :: Dvalues
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: Icount
    integer :: iel,ive,ivt
    
    ! Allocate temporal memory
    allocate(Dpoints(rtriangulationDest%ndim,rtriangulationDest%NNVE))
    allocate(Dvalues(rtriangulationDest%NNVE*rvectorDest%NVAR))
    allocate(Icount(rtriangulationDest%NVT)); Icount=0

    ! Set pointers
    call lsyssc_getbase_double(rvectorDest, p_Ddata)
    call storage_getbase_double2d(&
        rtriangulationDest%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2d(&
        rtriangulationDest%h_IverticesAtElement, p_IverticesAtElement)
    
    ! Loop over all elements
    do iel = 1, rtriangulationDest%NEL

      ! Collection physical coordinates of vertices
      do ive = 1, rtriangulationDest%NNVE
        Dpoints(:,ive) = p_DvertexCoords(:,p_IverticesAtElement(ive,iel))
      end do
      
      ! Evaluet FE-solution
      call fevl_evaluate_mult(DER_FUNC, Dvalues, rvectorSrc,&
          ImacroElements(iel), Dpoints=Dpoints)
      
      ! Store points in destination vector
      do ive = 1, rtriangulationDest%NNVE
        ivt = p_IverticesAtElement(ive,iel)
        p_Ddata(ivt) = p_Ddata(ivt) + Dvalues(ive)
        Icount(ivt)  = Icount(ivt) + 1
      end do
    end do

    ! Average nodal values
    do ivt = 1, rtriangulationDest%NVT
      p_Ddata(ivt) = p_Ddata(ivt)/real(Icount(ivt),DP)
    end do
      
    ! Deallocate temporal memory
    deallocate(Dpoints, Dvalues, Icount)
      
  end subroutine lin_evalAtDOFScalar

  !*****************************************************************************

!<subroutine>

  subroutine lin_prepareUnshareVertices(IverticesAtElement, InodalProperty,&
      NVT, NVBD, NMT, NEL)

!<description>
    ! This subroutine is used to prepare a triangulation for the duplication
    ! of all vertices which are shared by more than one element.
!</description>

!<input>
    ! Vertices at element list of source triangulation
    integer, dimension(:,:), intent(in) :: IverticesAtElement

    ! Nodal property list of source triangulation
    integer, dimension(:), intent(in) :: InodalProperty

    ! Number of elements of source triangulation
    integer, intent(in) :: NEL
!</input>

!<inputoutput>
    ! Total number of vertices
    ! On input: number of vertices of source triangulation
    ! On output: number of vertices of destinationtriangulation
    integer, intent(inout) :: NVT

    ! Number of vertices at the boundary
    ! On input: number of boundary vertices of source triangulation
    ! On output: number of boundary vertices of destination triangulation
    integer, intent(inout) :: NVBD
    
    ! Total number of edges
    ! On input: number of edges of source triangulation
    ! On output: number of edges of destinationtriangulation
    integer, intent(inout) :: NMT
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: iel,ive,ivt,NNVE
    
    ! All midpoints (=edges) located in the interior are
    ! duplicated. Thus, we can multiply them by two and substract
    ! those located at the boundary which are not duplicated
    if (NMT .gt. 0) NMT = 2*NMT-NVBD
    
    ! Get number of vertices per element
    NNVE  = size(IverticesAtElement,1)
    
    ! The total number of vertices is computed by summing the number
    ! of vertices per element (=NVE) over all elements
    NVT = 0; NVBD = 0
    
    ! Loop over all elements
    element: do iel = 1, NEL
      
      ! Loop over all vertices at the element
      do ive = 1, NNVE
        
        ! Get vertex number
        ivt = IverticesAtElement(ive,iel)
        if (ivt .eq. 0) cycle element
        
        ! Increment number of boundary vertices?
        if (InodalProperty(ivt) .gt. 0) NVBD = NVBD+1

        ! Increment number of total vertices
        NVT = NVT+1
      end do
    end do element
    
  end subroutine lin_prepareUnshareVertices

  !*****************************************************************************

!<subroutine>

  subroutine lin_unshareVertices(DvertexCoords, IverticesAtElement,&
      InodalProperty, NVT, NEL)

!<description>
    ! This subroutine duplicates all vertices of the triangulation
    ! which are shared by more than one element
!</description>

!<input>
    ! Number of elements
    integer, intent(in) :: NEL
!</input>

!<inputoutput>
    ! Vertex coordinates
    ! On input: vertex coordinates of source triangulation
    ! On output: vertex coordinates of destination triangulation
    ! DIMENSION(1:NDIM,1:NVT)
    real(DP), dimension(:,:), intent(inout) :: DvertexCoords
    
    ! Vertices at element list
    ! On input: vertices at element list of source triangulation
    ! On output: vertices at element list of destination triangulation
    ! DIMENSION(1:NVE,1:NEL)
    integer, dimension(:,:), intent(inout) :: IverticesAtElement
    
    ! Nodal property list
    ! On input: nodal property list of source triangulation
    ! On output: nodal property list of destination triangulation
    ! DIMENSION(1:NVT)
    integer, dimension(:), intent(inout) :: InodalProperty
    
    ! Number of vertices
    ! On input: number of vertices of source triangulation
    ! On output: number of vertices of destination triangulation
    integer, intent(inout) :: NVT
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), allocatable :: Icount
    integer :: iel,ive,ivt,NNVE
    
    ! Get the total number of vertices per element
    NNVE = size(IverticesAtElement,1)
    
    ! Allocate temporal memory for vertex counter
    allocate(Icount(NVT)); Icount = 0

    ! Loop over all elements of the source triangulation, duplicate
    ! all interior vertices and update the vertices-at-element list
    element: do iel = 1, NEL
      
      ! Loop over all vertices at the element
      do ive = 1, NNVE
        
        ! Get vertex number
        ivt = IverticesAtElement(ive,iel)
        if (ivt .eq. 0) cycle element
        
        ! Check if this vertex which has been visited before
        if (Icount(ivt) .gt. 0) then
          
          ! Increment number of vertices
          NVT = NVT+1
          
          ! Duplicate vertex
          InodalProperty(NVT)  = InodalProperty(ivt)
          DvertexCoords(:,NVT) = DvertexCoords(:,ivt)
          
          ! Update vertex in element list
          IverticesAtElement(ive,iel) = NVT
        end if
        
        ! Update vertex counter
        Icount(ivt) = Icount(ivt)+1
      end do
    end do element
    
    ! Deallocate temporal memory
    deallocate(Icount)
    
  end subroutine lin_unshareVertices

  !*****************************************************************************

!<subroutine>

  subroutine lin_convertQuad2Tri(DvertexCoords, IverticesAtElementSrc,&
      IverticesAtElementDest, IneighboursAtElement, ImacroElements,&
      NEL, binit)

!<description>
    ! This subroutine converts all quadrilaterals in the source
    ! triangulation into two triangles.
!</description>

!<input>
    ! Vertex coordinates
    real(DP), dimension(:,:), intent(in) :: DvertexCoords

    ! Vertices at element list of source triangulation
    integer, dimension(:,:), intent(in) :: IverticesAtElementSrc

    ! Flag to enable initialisation of mapping
    logical, intent(in) :: binit
!</input>

!<inputoutput>
    ! Vertices at element list
    ! On input: vertices of element list of source triangulation
    ! On output: vertices of element list of converted source triangulation
    integer, dimension(:,:), intent(inout) :: IneighboursAtElement

    ! Mapping between macro element in source triangulation and
    ! element in destination triangulation
    integer, dimension(:), intent(inout) :: ImacroElements
    
    ! Total number of elements
    ! On input: number of elements of source triangulation
    ! On output: number of elements of destination triangulation
    integer, intent(inout) :: NEL
!</inputoutput>

!<output>
    ! Vertices at element list of destination triangulation
    integer, dimension(:,:), intent(out) :: IverticesAtElementDest
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(2,3) :: Dpoints
    real(DP) :: darea1,darea2,dratio1,dratio2
    integer :: iel,jel,ive,jve,nve,NEL0

    ! Initialisation
    NEL0 = NEL
    NEL  = 0

    ! Initialse element mapping
    if (binit) then
      do iel = 1, NEL0
        ImacroElements(iel) = iel
      end do
    end if

    ! Loop over all elements of the source triangulation
    do iel = 1, NEL0
      
      ! Get number of vertices per element
      nve = merge(3, 4, IverticesAtElementSrc(4,iel) .eq. 0)
      
      ! Are we a quadrilateral element?
      if (nve .eq. 4) then
        
        ! Subdivide quadrilateral into two triangles such that
        ! the ratio of their areas is as close to unity as possible
        
        ! First triangle 1-2-3
        Dpoints(1:2,1) = DvertexCoords(1:2, IverticesAtElementSrc(1,iel))
        Dpoints(1:2,2) = DvertexCoords(1:2, IverticesAtElementSrc(2,iel))
        Dpoints(1:2,3) = DvertexCoords(1:2, IverticesAtElementSrc(3,iel))
        
        darea1 = gaux_getArea_tria2D(Dpoints)
        
        ! Second triangle 1-3-4
        Dpoints(1:2,1) = DvertexCoords(1:2, IverticesAtElementSrc(1,iel))
        Dpoints(1:2,2) = DvertexCoords(1:2, IverticesAtElementSrc(3,iel))
        Dpoints(1:2,3) = DvertexCoords(1:2, IverticesAtElementSrc(4,iel))
        
        darea2 = gaux_getArea_tria2D(Dpoints)
        
        if (darea1 .gt. darea2) then
          dratio1 = darea2/darea1
        else
          dratio1 = darea1/darea2
        end if

        ! First triangle 1-2-4
        Dpoints(1:2,1) = DvertexCoords(1:2, IverticesAtElementSrc(1,iel))
        Dpoints(1:2,2) = DvertexCoords(1:2, IverticesAtElementSrc(2,iel))
        Dpoints(1:2,3) = DvertexCoords(1:2, IverticesAtElementSrc(4,iel))
        
        darea1 = gaux_getArea_tria2D(Dpoints)
        
        ! Second triangle 2-3-4
        Dpoints(1:2,1) = DvertexCoords(1:2, IverticesAtElementSrc(2,iel))
        Dpoints(1:2,2) = DvertexCoords(1:2, IverticesAtElementSrc(3,iel))
        Dpoints(1:2,3) = DvertexCoords(1:2, IverticesAtElementSrc(4,iel))
        
        darea2 = gaux_getArea_tria2D(Dpoints)
        
        if (darea1 .gt. darea2) then
          dratio2 = darea2/darea1
        else
          dratio2 = darea1/darea2
        end if

        ! Increment counter
        NEL = NEL+1

        ! Set mapping to macro element
        ImacroElements(NEL0+NEL) = iel

        if (dratio1 .ge. dratio2) then
          ! Subdivide quadrilateral into triangles 1-2-3 and 1-3-4
          IverticesAtElementDest(1:3,iel) = IverticesAtElementSrc(1:3,iel)
          
          IverticesAtElementDest(1,NEL0+NEL) = IverticesAtElementSrc(1,iel)
          IverticesAtElementDest(2,NEL0+NEL) = IverticesAtElementSrc(3,iel)
          IverticesAtElementDest(3,NEL0+NEL) = IverticesAtElementSrc(4,iel)

          ! Update element numbe NEL0+NEL in the list of neighbouring
          ! elements in both possible element neighbours along edges
          ! (I3,I4) and (I4,I1)
          do ive = 3,4
            jel = IneighboursAtElement(ive,iel)
            if (jel .gt. 0) then
              do jve = 1,4
                if (IneighboursAtElement(jve,jel) .eq. iel) then
                  IneighboursAtElement(jve,jel) = NEL0+NEL
                  exit
                end if
              end do
            end if
          end do
          
          ! Update neighbours at elements
          IneighboursAtElement(1,NEL0+NEL) = iel
          IneighboursAtElement(2,NEL0+NEL) = IneighboursAtElement(3,iel)
          IneighboursAtElement(3,NEL0+NEL) = IneighboursAtElement(4,iel)
          IneighboursAtElement(4,NEL0+NEL) = 0

          IneighboursAtElement(3,iel) = NEL0+NEL
          IneighboursAtElement(4,iel) = 0

        else
          ! Subdivide quadrilateral into triangles 1-2-4 and 2-3-4
          IverticesAtElementDest(1,iel) = IverticesAtElementSrc(1,iel)
          IverticesAtElementDest(2,iel) = IverticesAtElementSrc(2,iel)
          IverticesAtElementDest(3,iel) = IverticesAtElementSrc(4,iel)
          
          IverticesAtElementDest(1:3,NEL0+NEL) = IverticesAtElementSrc(2:4,iel)

          ! Update element numbe NEL0+NEL in the list of neighbouring
          ! elements in both possible element neighbours along edges
          ! (I2,I3) and (I3,I4)
          do ive = 2,3
            jel = IneighboursAtElement(ive,iel)
            if (jel .gt. 0) then
              do jve = 1,4
                if (IneighboursAtElement(jve,jel) .eq. iel) then
                  IneighboursAtElement(jve,jel) = NEL0+NEL
                  exit
                end if
              end do
            end if
          end do

          ! Update neighbours at elements
          IneighboursAtElement(1,NEL0+NEL) = IneighboursAtElement(2,iel)
          IneighboursAtElement(2,NEL0+NEL) = IneighboursAtElement(3,iel)
          IneighboursAtElement(3,NEL0+NEL) = iel
          IneighboursAtElement(4,NEL0+NEL) = 0

          IneighboursAtElement(2,iel) = NEL0+NEL
          IneighboursAtElement(3,iel) = IneighboursAtElement(4,iel)
          IneighboursAtElement(4,iel) = 0
        end if
        
      else
        ! Copy the three corner vertices
        IverticesAtElementDest(1:3,iel) = IverticesAtElementSrc(1:3,iel)
      end if
      
    end do

    ! Update total number of elements
    NEL = NEL0 + NEL
    
  end subroutine lin_convertQuad2Tri

end module lineariser

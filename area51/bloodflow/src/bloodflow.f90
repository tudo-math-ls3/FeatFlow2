!##############################################################################
!# ****************************************************************************
!# <name> bloodflow </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main bloodflow module
!# </purpose>
!##############################################################################

module bloodflow

  use bilinearformevaluation
  use boundary
  use genoutput
  use hadaptaux
  use hadaptivity
  use linearsystemscalar
  use linearsystemblock
  use list
  use paramlist
  use spatialdiscretisation
  use storage
  use triangulation
  use ucd

  implicit none

  ! Every subroutine is declared private by default
  private

  ! Subroutine which are accessable from outside 
  ! this module are explicitly declared public
  public :: t_bloodflow
  public :: bloodflow_init
  public :: bloodflow_done
  public :: bloodflow_outputStructure
  public :: bloodflow_evalObject
  public :: bloodflow_adaptObject 

  !*****************************************************************************

!<constants>

!<constantblock description="Global constants for mesh spacing tolerances">

  ! Tolerance for considering two points as equivalent
  real(DP), parameter :: POINT_COLLAPSE_TOLERANCE = 1.0_DP/3.0_DP
  
!</constantblock>

!<constantblock description="Bitfield patterns for element states">

  ! Bitfield to identify the corners
  integer, parameter :: BITFIELD_POINT1 = ibset(0,1)
  integer, parameter :: BITFIELD_POINT2 = ibset(0,2)
  integer, parameter :: BITFIELD_POINT3 = ibset(0,3)

  ! Bitfield to identify the edges
  integer, parameter :: BITFIELD_EDGE1 = ibset(0,4)
  integer, parameter :: BITFIELD_EDGE2 = ibset(0,5)
  integer, parameter :: BITFIELD_EDGE3 = ibset(0,6)

  ! Bitfield used to check point intersection
  integer, parameter :: BITFIELD_POINT_INTERSECTION = BITFIELD_POINT1 +&
                                                      BITFIELD_POINT2 +&
                                                      BITFIELD_POINT3

  ! Bitfield used to cehck edge intersection
  integer, parameter :: BITFIELD_EDGE_INTERSECTION = BITFIELD_EDGE1 +&
                                                     BITFIELD_EDGE2 +&
                                                     BITFIELD_EDGE3

!</constantblock>

!</constants>

  !*****************************************************************************

!<types>

!<typeblock>
  
  ! This data structure contains all sub-structures for the bloodflow application

  type t_bloodflow

    ! Parameterlist structure
    type(t_parlist) :: rparlist

    ! Triangulation structure
    type(t_triangulation) :: rtriangulation

    ! Boundary parametrization
    type(t_boundary) :: rboundary

    ! Adaptation structure
    type(t_hadapt) :: rhadapt

    ! Handle to array storing the points of the object
    integer :: h_DobjectCoords = ST_NOHANDLE

    ! Indicator vector
    type(t_vectorScalar) :: rindicator    

    ! Number of object points which have no corresponding mesh point
    integer :: nunresolvedObjectPoints = 0
  end type t_bloodflow

!</typeblock>

!</types>

  !*****************************************************************************
   
contains

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_init(rbloodflow, sfilename, sdirectoryname)

!<description>

    ! This subroutine performs 
    !
    ! - all low-level initialization tasks of Featflow subsystems
    ! - reads in the parameter file(s) describing the bloodflow structure
    ! - reads the triangulation in the Featflow2 TRI format
    ! - reads the boundary parametrization in the Featflow2 PRM format

!</description>

!<input>

    ! File name
    character(len=*), intent(IN) :: sfilename

    ! Directory name
    character(len=*), intent(IN) :: sdirectoryname

!</input>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(INOUT) :: rbloodflow

!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(NDIM2D) :: Dorigin, Dlength
    character(len=SYS_STRLEN) :: strifilename, sprmfilename
    character(LEN=10) :: stime
    character(LEN=8)  :: sdate
    integer :: iglobRefLevel

    ! Initialize Feat2 subsystem
    call system_init()
  
    ! System will stop if critical error occurs
    sys_haltmode = SYS_HALT_THROWFPE

    ! Initialize the output system
    call date_and_time(sdate, stime)
    call output_init('./log/flagship_'//sdate//'_'//stime(1:4)//'.log')

    ! Initialize storage subsystem. The initial number of handles is
    ! set to 500. If this is not enough, then it will be increased
    ! repeatedly by 100 new handles.
    call storage_init(500, 100)

    ! Read the parameter file(s)
    call parlst_init(rbloodflow%rparlist)
    call parlst_readfromfile(rbloodflow%rparlist, sfilename, sdirectoryname)
    
    ! Read the boundary parametrization
    call parlst_getvalue_string(rbloodflow%rparlist, 'Input', 'prmfilename', sprmfilename)
    call boundary_read_prm(rbloodflow%rboundary, trim(adjustl(sprmfilename)))

    ! Read the triangulation
    call parlst_getvalue_string(rbloodflow%rparlist, 'Input', 'trifilename', strifilename)
    call tria_readTriFile2D(rbloodflow%rtriangulation,&
        trim(adjustl(strifilename)), rbloodflow%rboundary)

    ! Perform global refinement
    call parlst_getvalue_int(rbloodflow%rparlist, 'Input', 'iglobRefLevel', iglobRefLevel, 1)
    if (iglobRefLevel .gt. 1) then
      call tria_quickRefine2LevelOrdering(iglobRefLevel-1,&
          rbloodflow%rtriangulation, rbloodflow%rboundary)
    end if

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rbloodflow%rtriangulation,&
        rbloodflow%rboundary)

    ! Initialize adaptation structure
    call hadapt_initFromParameterlist(rbloodflow%rhadapt, rbloodflow%rparlist, 'Adaptation')
    
  end subroutine bloodflow_init

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_done(rbloodflow)

!<description>

    ! This subroutine performs all clean-up tasks

!</description>

!<inputoutput>

    ! OPTIONAL: Bloodflow structure
    type(t_bloodflow), intent(INOUT), optional :: rbloodflow
!</inputoutput>

!</subroutine>
    
    ! Release the bloodflow structure
    if (present(rbloodflow)) then

      ! Release the parameter list
      call parlst_done(rbloodflow%rparlist)

      ! Release the boundary parametrization
      call boundary_release(rbloodflow%rboundary)
      
      ! Release the triangulation structure
      call tria_done(rbloodflow%rtriangulation)

      ! Release the object coordinates
      call storage_free(rbloodflow%h_DobjectCoords)

      ! Release indicator vector
      call lsyssc_releaseVector(rbloodflow%rindicator)

      ! Release adaptation structure
      call hadapt_releaseAdaptation(rbloodflow%rhadapt)

    end if

    ! Release storage
    call storage_info(.true.)
    call storage_done()
    call output_lbrk()
    
    ! Close logfile
    call output_done()
    
  end subroutine bloodflow_done

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_outputStructure(rbloodflow)

!<description>

    ! This subroutine writes the content of the structure to a GMV file

!</description>

!<input>

    ! Bloodflow structure
    type(t_bloodflow), intent(IN) :: rbloodflow

!</input>
!</subroutine>

    ! UCD export structure
    type(t_ucdExport) :: rexport
    real(DP), dimension(:,:), pointer :: p_Ddata2d
    real(DP), dimension(:), pointer :: p_Ddata, Dtracer
    character(len=SYS_STRLEN) :: sucdfilename

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! Get filename and start GMV output
    call parlst_getvalue_string(rbloodflow%rparlist, 'Output', 'ucdfilename', sucdfilename)
    call ucd_startGMV(rexport, UCD_FLAG_STANDARD, rbloodflow%rtriangulation,&
        trim(sucdfilename)//'.'//trim(sys_si0(ifilenumber,5)))

    ! Attach the indicator vector
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Ddata)
    call ucd_addVariableElementBased(rexport, 'Indicator', UCD_VAR_STANDARD, p_Ddata)

    ! Attach the thin object as tracer
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_Ddata2d)
    call ucd_setTracers (rexport, p_Ddata2d)
    
    allocate(Dtracer(size(p_Ddata2d, 2)))
    Dtracer = 1
    call ucd_addTracerVariable (rexport,'objectpoints',Dtracer)
    deallocate(Dtracer)
    
    ! Write UCD exporter to file
    call ucd_write(rexport)
    call ucd_release(rexport)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

  end subroutine bloodflow_outputStructure

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_evalObject(rbloodflow, dtime)

!<description>

    ! This subroutine evaluates the location of the thin object

!</description>

!<input>

    ! Simulation time
    real(DP), intent(IN) :: dtime

!</input>

!<inputoutput>
    ! Bloodflow structure
    type(t_bloodflow), intent(INOUT) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DobjectCoords
    real(DP) :: c,L
    integer, dimension(2) :: Isize
    integer :: ipoint, npoints

    ! Release thin object from previous evaluation
    if (rbloodflow%h_DobjectCoords .ne. ST_NOHANDLE) then
      call storage_free(rbloodflow%h_DobjectCoords)
    end if

    ! Get values from parameter list
    call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'L',       L)
    call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'c',       c)
    call parlst_getvalue_int(rbloodflow%rparlist,    'Object', 'npoints', npoints)

    ! Generate list of vertices
    Isize = (/2, npoints/)
    call storage_new('bloodflow_evalObject', 'DobjectCoords', Isize,&
                     ST_DOUBLE, rbloodflow%h_DobjectCoords, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)
    
    do ipoint = 0, npoints-1
      
      ! Compute x-coordinate
      p_DobjectCoords(1,ipoint+1) = L*ipoint/real(npoints-1, DP)

      ! Compute y-coordinate
      p_DobjectCoords(2,ipoint+1) = c*sin(dtime)*(3*L*p_DobjectCoords(1,ipoint+1)**2 -&
                                                      p_DobjectCoords(1,ipoint+1)**3)
    end do
    
  end subroutine bloodflow_evalObject

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_adaptObject(rbloodflow)

!<description>
    
    ! This subroutine adapts the computational grid to the object.

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(INOUT) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iadapt

    ! Check if adaptation structure has been prepared
    if (rbloodflow%rhadapt%iSpec .eq. HADAPT_HAS_PARAMETERS) then
      ! Initialize adaptation structure from triangulation
      call hadapt_initFromTriangulation(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    else
      ! Refresh adaptation structure
      call hadapt_refreshAdaptation(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    end if

!!$    do iadapt = 1, 100
!!$
!!$      ! Evaluate the indicator
    call bloodflow_evalIndicator(rbloodflow)
!!$      
!!$      ! Check if there are still unresolved object points
!!$      if (rbloodflow%nunresolvedObjectPoints .eq. 0) return
!!$      
!!$      ! If so, then adapt the computational mesh
!!$      call hadapt_performAdaptation(rbloodflow%rhadapt, rbloodflow%rindicator)
!!$      
!!$      ! Generate raw mesh from adaptation structure
!!$      call hadapt_generateRawMesh(rbloodflow%rhadapt, rbloodflow%rtriangulation)
!!$    end do

!!$    ! If we end up here, then the maximum admissible refinement level
!!$    ! does not suffice to resolve all object points
!!$    call output_line('Some objects points could not be resolved!',&
!!$        OU_CLASS_ERROR,OU_MODE_STD,'bloodflow_adaptObject')
!!$    call output_line('Increase the maximum admissible refinement level!',&
!!$        OU_CLASS_ERROR,OU_MODE_STD,'bloodflow_adaptObject')
!!$    call sys_halt()

  end subroutine bloodflow_adaptObject

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_evalIndicator(rbloodflow)

!<description>

    ! This subroutine evaluates the indicator function

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(INOUT) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    type(t_list) :: relementList
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dindicator
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex, p_Iindicator
    integer, dimension(:), pointer :: IelementPatch
    integer, dimension(2) :: Isize
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(TRIA_NVETRI2D) :: DbarycCoords
    real(DP), dimension(NDIM2D) :: Dpoint0Coords, Dpoint1Coords,DIntersectionCoords
    integer :: ive,jve,iel,jel,i,i1,i2,i3,istatus,idx,ipoint,ipatch,npatch,ipos
  
    
    ! Release the indicator vector and create new one
    call lsyssc_releaseVector(rbloodflow%rindicator)
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Dindicator)
    
    ! Create new indicator vector as integer array
    call lsyssc_createVector(rindicator, rbloodflow%rtriangulation%NEL, .true., ST_INT)
    call lsyssc_getbase_int(rindicator, p_Iindicator)

    ! Allocate temporal memory for local element patch
    allocate(IelementPatch(rbloodflow%rtriangulation%NNelAtVertex))

    ! Create linked list for storing the elements adjacent to the object
    call list_createList(relementList, ceiling(0.1*rbloodflow%rtriangulation%NEL), ST_INT, 0, 0, 0)

    ! Set pointers
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call storage_getbase_int(&
        rbloodflow%rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
    call storage_getbase_int(&
        rbloodflow%rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
    call storage_getbase_double2d(&
        rbloodflow%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)

    !---------------------------------------------------------------------------
    ! (1) Find element surrounding/meeting at first point of the object:
    !
    !     The algorithm is really simply. An extensive search over all
    !     element of the triangulation is performed and the first
    !     element which either surrounds the first point of the thin
    !     object or is connected to this point via a corner certex or
    !     an edge is selected.
    !---------------------------------------------------------------------------
    
    ! Initialize point number
    ipoint = 1
    
    ! Get global coordinates of first point
    Dpoint0Coords =  p_DobjectCoords(:,ipoint)

    ! Loop over all elements in triangulation
    findElementOfFirstVertex: do iel = 1, rbloodflow%rtriangulation%NEL

      ! Get vertices at element
      i1 = p_IverticesAtElement(1, iel)
      i2 = p_IverticesAtElement(2, iel)
      i3 = p_IverticesAtElement(3, iel)

      ! Get global coordinates of corner vertices
      DtriaCoords(:,1) = p_DvertexCoords(:,i1)
      DtriaCoords(:,2) = p_DvertexCoords(:,i2)
      DtriaCoords(:,3) = p_DvertexCoords(:,i3)
      
      ! Check if the point is 'inside' or on the boundary of the element
      call PointInTriangleTest(DtriaCoords, Dpoint0Coords, 0.0_DP, istatus)
      if (istatus .ge. 0) exit findElementOfFirstVertex
      
    end do findElementOfFirstVertex

    ! Update point number
    ipoint = ipoint+1

    !---------------------------------------------------------------------------
    ! (2) Find elements surrounding/meeting at all others points of the object:
    !
    !     This algorithm is slightly more complicated. The list of
    !     points on the thin object is visited segment-by-segment. If
    !     the starting point (aka previous point) of the segment is
    !     located inside an element, then nothing needs to be done. If
    !     the previous point was located on an edge, then the opposite
    !     element is also marked. Finally, if the previous point
    !     coincides with some corner vertex, then all elements meeting
    !     at that point are marked.
    !
    !     Next, we proceed to the endpoint of the segment. If it lies
    !     inside the same element, then the above procedure applies to
    !     the endpoint. Otherwise, the segment must either intersect
    !     one edge of the current element or run through a corner
    !     vertex. In any case, we proceed to the adjacent element(s)
    !     and check if the end point is located inside that element or
    !     coincides with one of its edges/corners. This process
    !     continues until the endpoint of the last segment is reached.
    !
    !---------------------------------------------------------------------------

    ! Loop over all remaining points
    findElementsOfOtherVertices: do while (ipoint .le. size(p_DobjectCoords, 2))
      
      ! Create single element patch
      npatch = 1
      IelementPatch(1) = iel

      ! Append element to the list of elements adjacent to the object
      if (p_Iindicator(iel) .eq. 0)&
          call list_appendToList(relementList, iel, ipos)
      
      ! Mark element for potential refinement
      p_Iindicator(iel) = ibset(p_Iindicator(iel), istatus)
      
      ! Check status of previous point
      select case(istatus)       
      case (1:3)
        ! Previous point was one of the corner vertices, thus, mark
        ! all elements meeting at that corner and create local patch
        i = p_IverticesAtElement(istatus, iel); ipatch = 0

        do idx = p_IelementsAtVertexIdx(i),&
                 p_IelementsAtVertexIdx(i+1)-1
         
          jel = p_IelementsAtVertex(idx)
          if (p_Iindicator(jel) .eq. 0)&
              call list_appendToList(relementList, jel, ipos)
          
          ipatch = ipatch+1
          IelementPatch(ipatch) = jel

          do ive = 1, 3
            if (p_IverticesAtElement(ive, jel) .eq. i) then
              p_Iindicator(jel) = ibset(p_Iindicator(jel), ive)
              exit
            end if
          end do

        end do
        npatch = ipatch
        
      case (4:6)               
        ! Previous point was located on the edge of the element, thus,
        ! mark adjacent element and create two element patch
        jel = p_IneighboursAtElement(istatus-3, iel)
        
        if (jel .ne. 0) then

          if (p_Iindicator(jel) .eq. 0)&
              call list_appendToList(relementList, jel, ipos)

          npatch = 2
          IelementPatch(2) = jel
          
          do ive = 1, 3
            if (p_IneighboursAtElement(ive, jel) .eq. iel) then
              p_Iindicator(jel) = ibset(p_Iindicator(jel), ive+3)
              exit
            end if
          end do

        end if
        
      end select
      
      ! Get global coordinates of the endpoint
      Dpoint1Coords =  p_DobjectCoords(:,ipoint)
      
      ! Loop over all elements in patch
      findInPatchDirect: do ipatch = 1, npatch
        
        ! Get element number
        iel = IelementPatch(ipatch)
        
        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        ! Check if the endpoint is 'inside' or on the boundary of the current element
        call PointInTriangleTest(DtriaCoords, Dpoint1Coords, 0.0_DP, istatus)
        if (istatus .ge. 0) then
          
          ! If so, use the current point as new starting point of the segment
          Dpoint0Coords = Dpoint1Coords
          
          ! Update point number
          ipoint = ipoint+1
          
          ! That's it, we can forget about the current patch of elements
          cycle findElementsOfOtherVertices
          
        end if
      end do findInPatchDirect
      
      ! Ok, the endpoint was not in one of the elements belonging to
      ! the current patch. Hence, the segment connecting the start and
      ! endpoint must cross some boundary edge of the patch

      ! Loop over all elements in patch
      findInPatchIndirect: do ipatch = 1, npatch
        
        ! Get element number
        iel = IelementPatch(ipatch)

        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        ! Loop over all edges of the element
        findIntersectionEdge: do ive = 1, 3

          ! Check if the edge intersects with the line between start
          ! and endpoint and compute coordinates of intersection point
          call LinesIntersectionTest(Dpoint0Coords, Dpoint1Coords,&
                                     DtriaCoords(:,ive), DtriaCoords(:,mod(ive,3)+1),&
                                     istatus, DintersectionCoords)
          if (istatus .eq. 0) then
            
            ! Get global element number of adjacent element
            iel = p_IneighboursAtElement(ive, iel)
            
            ! Get vertices at element
            i1 = p_IverticesAtElement(1, iel)
            i2 = p_IverticesAtElement(2, iel)
            i3 = p_IverticesAtElement(3, iel)
            
            ! Get global coordinates of corner vertices
            DtriaCoords(:,1) = p_DvertexCoords(:,i1)
            DtriaCoords(:,2) = p_DvertexCoords(:,i2)
            DtriaCoords(:,3) = p_DvertexCoords(:,i3)
            
            ! Check if the endpoint is 'inside' or on the boundary of the adjacent element
            call PointInTriangleTest(DtriaCoords, Dpoint1Coords, 0.0_DP, istatus)
            
            if (istatus .eq. -1) then
              
              ! If the point is 'outside' the element, then the
              ! segment just passes through the element, whereby no
              ! point falls into it. In this case, the intersection
              ! point serves as new starting point of the segment
              Dpoint0Coords = DintersectionCoords

              ! Perform point-in-triangle test for intersection point.
              call PointInTriangleTest(DtriaCoords, Dpoint0Coords, 0.0_DP, istatus)
              
            else
              
              ! If the endpoint is 'inside' the adjacent element, then
              ! we can use it as new starting point of the next segment
              Dpoint0Coords = Dpoint1Coords

              ! Update point number
              ipoint = ipoint+1

            end if
            
            ! In any case, we have found a new element
            cycle findElementsOfOtherVertices
          end if
          
        end do findIntersectionEdge

      end do findInPatchIndirect
      
    end do findElementsOfOtherVertices
    
    ! Append element to the list of elements adjacent to the object
    if (p_Iindicator(iel) .eq. 0)&
        call list_appendToList(relementList, iel, ipos)
    
    ! Mark last element for potential refinement
    p_Iindicator(iel) = ibset(p_Iindicator(iel), istatus)
    
    ! Check status of previous point
    select case(istatus)
    case (1:3)
      ! Previous point was one of the corner vertices, thus, mark
      ! all elements meeting at that corner and create local patch
      i = p_IverticesAtElement(istatus, iel)
      
      do idx = p_IelementsAtVertexIdx(i),&
               p_IelementsAtVertexIdx(i+1)-1
        
        jel = p_IelementsAtVertex(idx)
        if (p_Iindicator(jel) .eq. 0)&
            call list_appendToList(relementList, jel, ipos)

        do ive = 1, 3
          if (p_IverticesAtElement(ive, jel) .eq. i) then
            p_Iindicator(jel) = ibset(p_Iindicator(jel), ive)
            exit
          end if
        end do

      end do
      
    case (4:6)
      ! Previous point was located on the edge of the element, thus,
      ! mark adjacent element and create two element patch
      jel = p_IneighboursAtElement(istatus-3, iel)
      
      if (p_Iindicator(jel) .eq. 0)&
          call list_appendToList(relementList, jel, ipos)
      
      do ive = 1, 3
        if (p_IneighboursAtElement(ive, jel) .eq. iel) then
          p_Iindicator(jel) = ibset(p_Iindicator(jel), ive+3)
          exit
        end if
      end do

    end select

    ! Deallocate temporal memory
    deallocate(IelementPatch)
    
    goto 999

    !---------------------------------------------------------------------------
    ! (3) Loop over all segments/elements and decide whether to refine or
    !     to reposition mesh points so that elements get aligned to the object
    !---------------------------------------------------------------------------

    ! Loop over all elements adjacent to the object
    ipos = list_getNextInList(relementList, .true.)
    do while(ipos .ne. LNULL)
      
      ! Get element number and proceed to next position
      call list_getByPosition(relementList, ipos, iel)
      ipos = list_getNextInList(relementList, .false.)

      ! Check status of elements
      select case(iand(p_Iindicator(iel), BITFIELD_EDGE_INTERSECTION))
        
        case(BITFIELD_EDGE1)
          
          ! Get vertices at element
          i1 = p_IverticesAtElement(1, iel)
          i2 = p_IverticesAtElement(2, iel)
          i3 = p_IverticesAtElement(3, iel)

          ! Get global coordinates of corner vertices
          DtriaCoords(:,1) = p_DvertexCoords(:,i1)
          DtriaCoords(:,2) = p_DvertexCoords(:,i2)
          DtriaCoords(:,3) = p_DvertexCoords(:,i3)

          ! Loop over all segments and check intersection point
          edge1: do ipoint = 1, size(p_DobjectCoords,2)-1

            ! Perform line intersection test
            call LinesIntersectionTest(p_DobjectCoords(:,ipoint),&
                                       p_DobjectCoords(:,ipoint+1),&
                                       DtriaCoords(:,1), DtriaCoords(:,2),&
                                       istatus, DintersectionCoords)
            if (istatus .eq. 1) then
              ! Remove the refinement indicator
              p_Iindicator(iel) = 0

              ! That's it
              exit edge1
            elseif (istatus .eq. 0) then

              ! Perform point in triangle test
              call PointInTriangleTest(DtriaCoords, DintersectionCoords,&
                                       POINT_COLLAPSE_TOLERANCE, istatus)
              select case(istatus)
              case(1)
                p_DvertexCoords(:, i1) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case(2)
                p_DvertexCoords(:, i2) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case(3)
                p_DvertexCoords(:, i3) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case default
                p_Iindicator(iel)      = 1
              end select

              ! That's it
              exit edge1
            end if
          end do edge1

        case(BITFIELD_EDGE2)
          
          ! Get vertices at element
          i1 = p_IverticesAtElement(1, iel)
          i2 = p_IverticesAtElement(2, iel)
          i3 = p_IverticesAtElement(3, iel)

          ! Get global coordinates of corner vertices
          DtriaCoords(:,1) = p_DvertexCoords(:,i1)
          DtriaCoords(:,2) = p_DvertexCoords(:,i2)
          DtriaCoords(:,3) = p_DvertexCoords(:,i3)

          ! Loop over all segments and check intersection point
          edge2: do ipoint = 1, size(p_DobjectCoords,2)-1

            ! Perform line intersection test
            call LinesIntersectionTest(p_DobjectCoords(:,ipoint),&
                                       p_DobjectCoords(:,ipoint+1),&
                                       DtriaCoords(:,2), DtriaCoords(:,3),&
                                       istatus, DintersectionCoords)
            if (istatus .eq. 1) then
              ! Remove the refinement indicator
              p_Iindicator(iel) = 0

              ! That's it
              exit edge2
            elseif (istatus .eq. 0) then

              ! Perform point in triangle test
              call PointInTriangleTest(DtriaCoords, DintersectionCoords,&
                                       POINT_COLLAPSE_TOLERANCE, istatus)
              select case(istatus)
              case(1)
                p_DvertexCoords(:, i1) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case(2)
                p_DvertexCoords(:, i2) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case(3)
                p_DvertexCoords(:, i3) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case default
                p_Iindicator(iel)      = 1
              end select
             
              ! That's it
              exit edge2
            end if
          end do edge2

        case(BITFIELD_EDGE3)

          ! Get vertices at element
          i1 = p_IverticesAtElement(1, iel)
          i2 = p_IverticesAtElement(2, iel)
          i3 = p_IverticesAtElement(3, iel)

          ! Get global coordinates of corner vertices
          DtriaCoords(:,1) = p_DvertexCoords(:,i1)
          DtriaCoords(:,2) = p_DvertexCoords(:,i2)
          DtriaCoords(:,3) = p_DvertexCoords(:,i3)

          ! Loop over all segments and check intersection point
          edge3: do ipoint = 1, size(p_DobjectCoords,2)-1

            ! Perform line intersection test
            call LinesIntersectionTest(p_DobjectCoords(:,ipoint),&
                                       p_DobjectCoords(:,ipoint+1),&
                                        DtriaCoords(:,1), DtriaCoords(:,3),&
                                       istatus, DintersectionCoords)
            if (istatus .eq. 1) then
              ! Remove the refinement indicator
              p_Iindicator(iel) = 0

              ! That's it
              exit edge3
            elseif (istatus .eq. 0) then

              ! Perform point in triangle test
              call PointInTriangleTest(DtriaCoords, DintersectionCoords,&
                                       POINT_COLLAPSE_TOLERANCE, istatus)
              select case(istatus)
              case(1)
                p_DvertexCoords(:, i1) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case(2)
                p_DvertexCoords(:, i2) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case(3)
                p_DvertexCoords(:, i3) = DintersectionCoords
                p_Iindicator(iel)      = 0
              case default
                p_Iindicator(iel)      = 1
              end select

              ! That's it
              exit edge3
            end if
          end do edge3
                    
        case(BITFIELD_EDGE1 + BITFIELD_EDGE2,&
             BITFIELD_EDGE2 + BITFIELD_EDGE3,&
             BITFIELD_EDGE1 + BITFIELD_EDGE3)
          print *, "Two edges"

        case (BITFIELD_EDGE1 + BITFIELD_EDGE2 + BITFIELD_EDGE3)
          ! Set the refinement indicator to unity
          p_Iindicator(iel) = 1

        case default
          ! Remove the refinement indicator
          p_Iindicator(iel) = 0
        end select
      
      end do

    ! Release list of elements
999   call list_releaseList(relementList)

    ! Change the data type of the indicator
    do i = 1, size(p_Dindicator)
      p_Dindicator(i) = min(p_Iindicator(i), 1)
    end do
    call lsyssc_releaseVector(rindicator)

  contains
    
    !***************************************************************************
    
    pure subroutine PointInTriangleTest(TriaCoords, P, dtolerance, istatus)

      ! This subroutine calculates the relation of the given point P
      ! compare to the triangle which is defined by its three corners
      ! The meaning of the resulting istatus is as follows:
      !
      ! istatus:
      !  = -1 : point P is outside of the triangle
      !  =  0 : point P is located inside the triangle
      !  =  1 : point P is equivalent to point A of the triangle
      !  =  2 : point P is equivalent to point B of the triangle
      !  =  3 : point P is equivalent to point C of the triangle
      !  =  4 : point P is located in the edge between points A and B
      !  =  5 : point P is located in the edge between points B and C
      !  =  6 : point P is located in the edge between points C and A
      
      real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(IN) :: TriaCoords
      real(DP), dimension(NDIM2D), intent(IN) :: P
      real(DP), intent(IN) :: dtolerance

      integer, intent(OUT) :: istatus


      ! local variables
      real(DP) :: area,area1,area2,area3,u,v,w

      ! Compute area of global and sub-triangles
      area  = signedArea(TriaCoords(:,1), TriaCoords(:,2), TriaCoords(:,3))
      area1 = signedArea(P,               TriaCoords(:,2), TriaCoords(:,3))
      area2 = signedArea(P,               TriaCoords(:,3), TriaCoords(:,1))
      area3 = signedArea(P,               TriaCoords(:,1), TriaCoords(:,2))

      ! Compute barycentric coordinates
      u = area1/area                             
      v = area2/area
      w = area3/area

      ! Determine status
      if ((u .lt. 0) .or. (v .lt. 0) .or. (w .lt. 0) .or. (u+v+w .gt. 1)) then
        ! If the barycentric coordinates are negative of larger than
        ! one then the point is located outside of the triangle
        istatus = -1

      else
        ! Otherwise, the point is located inside the triangle. We have
        ! to considere several cases, e.g., the point is close to an
        ! edge or it is close to a vertex and collapses
        if (u .ge. 1-dtolerance) then
          ! Point P collapses with point A
          istatus = 1
        elseif (v .ge. 1-dtolerance) then
          ! Point P collapses with point B
          istatus = 2
        elseif (w .ge. 1-dtolerance) then
          ! Point P collapses with point C
          istatus = 3
        elseif (w .le. 0.5*dtolerance) then
          ! Point P collapses with edge AB
          istatus = 4
        elseif (u .le. 0.5*dtolerance) then
          ! Point P collapses with edge BC
          istatus = 5
        elseif (v .le. 0.5*dtolerance) then
          ! Point P collapses with edge AC
          istatus = 6
        else
          ! Point P is in the interior
          istatus = 0
        end if
      end if
      
    end subroutine PointInTriangleTest

    !***************************************************************************

    subroutine LinesIntersectionTest(P1, P2, P3, P4, istatus, S)

      ! This subroutine checks if the two line segments P1-P2 and
      ! P3-P4 intersect each other and returns the intersection point.
      ! Note that the start point is excluded from the
      ! segment. Otherwise, the starting point would by considered as
      ! intersection point if it was located on the edge.
      !
      ! The algorithm is taken from:
      ! Paul Bourke, Intersection Point of Two Lines (2 Dimensions)
      ! http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
      !
      ! The meaning of the resulting istatus is as follows:
      !
      ! istatus:
      !  = -1 : the two line segments do not intersect
      !  =  0 : the two line segments intersect each other
      !  =  1 : the two line segments partially coincide
            
      real(DP), dimension(NDIM2D), intent(IN) :: P1,P2,P3,P4

      real(DP), dimension(NDIM2D), intent(OUT) :: S
      integer, intent(OUT) :: istatus

      ! local variables
      real(DP) :: denom,nom1,nom2,u,v


      ! Compute (de-)nominators
      denom = (P4(2)-P3(2)) * (P2(1)-P1(1)) - (P4(1)-P3(1)) * (P2(2)-P1(2))
      nom1  = (P4(1)-P3(1)) * (P1(2)-P3(2)) - (P4(2)-P3(2)) * (P1(1)-P3(1))
      nom2  = (P2(1)-P1(1)) * (P1(2)-P3(2)) - (P2(2)-P1(2)) * (P1(1)-P3(1))

      ! Check if lines are parallel
      if (abs(denom) .le. SYS_EPSREAL) then

        ! The two lines are parallel
        
        if ( (abs(nom1) .le. SYS_EPSREAL) .or.&
             (abs(nom2) .le. SYS_EPSREAL) ) then
          istatus =  1
        else
          istatus = -1
        end if
        S = SYS_INFINITY
        
      else
        
        ! The two lines are not parallel, hence they must intersect each other

        ! Compute parameter values
        u = nom1 / denom
        v = nom2 / denom
        
        ! Check if the intersection point is located within the line segments
        if ((u .le. sqrt(SYS_EPSREAL)) .or. (u .gt. 1)&
            .or. (v .lt. 0) .or. (v .gt. 1)) then
          
          ! Intersection point does not belong to both line segments
          istatus = -1          
          S = SYS_INFINITY
          
        else
          
          ! Intersection point belongs to both line segments
          istatus = 0
          
          ! Compute coordinates of intersection point
          S(1) = P1(1) + u*(P2(1)-P1(1))
          S(2) = P1(2) + u*(P2(2)-P1(2))
        end if

      end if

    end subroutine LinesIntersectionTest

    !***************************************************************************
    
    pure subroutine getBarycentricCoords(TriaCoords, P, BarycentricCoords)

      ! This subroutine calculates the barycentric coordinates of the
      ! given point P using the given reference triangle.

      real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(IN) :: TriaCoords
      real(DP), dimension(NDIM2D), intent(IN) :: P

      real(DP), dimension(TRIA_NVETRI2D), intent(OUT) :: BarycentricCoords

      
      ! local variables
      real(DP) :: area,area1,area2,area3
      
      
      ! Compute area of global and sub-triangles
      area  = signedArea(TriaCoords(:,1), TriaCoords(:,2), TriaCoords(:,3))
      area1 = signedArea(P,               TriaCoords(:,2), TriaCoords(:,3))
      area2 = signedArea(P,               TriaCoords(:,3), TriaCoords(:,1))
      area3 = signedArea(P,               TriaCoords(:,1), TriaCoords(:,2))

      ! Compute barycentric coordinates
      BarycentricCoords(1) = area1/area                             
      BarycentricCoords(2) = area2/area
      BarycentricCoords(3) = area3/area

    end subroutine getBarycentricCoords

    !***************************************************************************
    
    pure function signedArea(P1,P2,P3) result(area)

      ! This subroutine computes the signed area of the triangle
      
      real(DP), dimension(NDIM2D), intent(IN) :: P1,P2,P3
      real(DP) :: area

      area = 0.5 * ( (P2(1)-P1(1))*(P3(2)-P1(2)) -&
                     (P3(1)-P1(1))*(P2(2)-P1(2)) )

    end function signedArea

  end subroutine bloodflow_evalIndicator

end module bloodflow

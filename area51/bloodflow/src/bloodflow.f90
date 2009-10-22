!##############################################################################
!# ****************************************************************************
!# <name> bloodflow </name>
!# ****************************************************************************
!#
!# <purpose>
!# This is the main bloodflow module
!# </purpose>
!#
!# The following routines are available:
!#
!# 1.) bloodflow_init
!#     -> Initializes the bloodflow structure
!#
!# 2.) bloodflow_done
!#     -> Releases the bloodflow structure
!#
!# 3.) bloodflow_outputStructure
!#     -> Exports the content of the bloodflow structure to UCD file
!#
!# 4.) bloodflow_evalObject
!#     -> Evaluates the object at a given time step
!#
!# 5.) bloodflow_adaptObject
!#     -> Adaps the computational mesh to the object
!#
!# 6.) bloodflow_evalIndicator
!#     -> Evaluates the indicator vector based on the object
!#
!# The following auxiliary routines are available:
!#
!# 1.) PointInTriangleTest
!#     -> Tests if a point is located inside a triangle
!#
!# 2.) LinesIntersectionTest
!#     -> Tests if two line segments intersect each other
!#
!# 3.) getBarycentricCoords
!#     -> Calculates the barycentric coordinates of a point w.r.t. triangle
!#
!# 4.) signedArea
!#     -> Calculates the signed area of a triangle
!#
!##############################################################################

module bloodflow

  use bilinearformevaluation
  use boundary
  use genoutput
  use hadaptaux
  use hadaptaux2d
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
  public :: bloodflow_evalIndicator

  !*****************************************************************************

!<constants>

!<constantblock description="Global constants for mesh spacing tolerances">

  ! Tolerance for considering two points as equivalent
  real(DP), parameter :: POINT_EQUAL_TOLERANCE = 1e-6_DP

  ! Tolerance for considering two points as equivalent
  real(DP), parameter :: POINT_COLLAPSE_TOLERANCE = 1.0/3.0
  
!</constantblock>

!<constantblock description="Bitfield patterns for element states">

  ! Bitfield to identify an interior vertex
  integer, parameter :: BITFIELD_INNER = ibset(0,0)

  ! Bitfield to identify the corners
  integer, parameter :: BITFIELD_POINT1 = ibset(0,1)
  integer, parameter :: BITFIELD_POINT2 = ibset(0,2)
  integer, parameter :: BITFIELD_POINT3 = ibset(0,3)

  ! Bitfield to identify the edges
  integer, parameter :: BITFIELD_EDGE1 = ibset(0,4)
  integer, parameter :: BITFIELD_EDGE2 = ibset(0,5)
  integer, parameter :: BITFIELD_EDGE3 = ibset(0,6)

  ! Bitfield to identify elements in list
  integer, parameter :: BITFIELD_INLIST = ibset(0,7)

  ! Bitfield to identify multi-intersected edges
  integer, parameter :: BITFIELD_MULTI_INTERSECTION = ibset(0,8)

  ! Bitfield used to check point intersection
  integer, parameter :: BITFIELD_POINT_INTERSECTION = BITFIELD_POINT1 +&
                                                      BITFIELD_POINT2 +&
                                                      BITFIELD_POINT3

  ! Bitfield used to check edge intersection
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

    ! Triangulation structure of base mesh
    type(t_triangulation) :: rtriangulationBase

    ! Boundary parametrization
    type(t_boundary) :: rboundary

    ! Adaptation structure
    type(t_hadapt) :: rhadapt

    ! Handle to array storing the points of the object
    integer :: h_DobjectCoords = ST_NOHANDLE

    ! Indicator vector
    type(t_vectorScalar) :: rindicator    

    ! List of elements intersected by object
    type(t_list) :: relementList

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
    character(len=*), intent(in) :: sfilename

    ! Directory name
    character(len=*), intent(in) :: sdirectoryname

!</input>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>
    
    ! local variables
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
    call output_init('./log/bloodflow_'//sdate//'_'//stime(1:4)//'.log')

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
    call tria_readTriFile2D(rbloodflow%rtriangulationBase,&
        trim(adjustl(strifilename)), rbloodflow%rboundary)

    ! Perform global refinement
    call parlst_getvalue_int(rbloodflow%rparlist, 'Input', 'iglobRefLevel', iglobRefLevel, 1)
    if (iglobRefLevel .gt. 1) then
      call tria_quickRefine2LevelOrdering(iglobRefLevel-1,&
          rbloodflow%rtriangulationBase, rbloodflow%rboundary)
    end if

    ! Generate standard mesh from raw mesh
    call tria_initStandardMeshFromRaw(rbloodflow%rtriangulationBase,&
        rbloodflow%rboundary)
    
  end subroutine bloodflow_init

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_done(rbloodflow)

!<description>

    ! This subroutine performs all clean-up tasks

!</description>

!<inputoutput>

    ! OPTIONAL: Bloodflow structure
    type(t_bloodflow), intent(inout), optional :: rbloodflow
!</inputoutput>

!</subroutine>
    
    ! Release the bloodflow structure
    if (present(rbloodflow)) then

      ! Release the parameter list
      call parlst_done(rbloodflow%rparlist)

      ! Release the boundary parametrization
      call boundary_release(rbloodflow%rboundary)
      
      ! Release the triangulation structures
      call tria_done(rbloodflow%rtriangulationBase)
      call tria_done(rbloodflow%rtriangulation)

      ! Release the object coordinates
      call storage_free(rbloodflow%h_DobjectCoords)

      ! Release indicator vector
      call lsyssc_releaseVector(rbloodflow%rindicator)

      ! Release list of elements
      call list_releaseList(rbloodflow%relementList)

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
    type(t_bloodflow), intent(in) :: rbloodflow

!</input>
!</subroutine>

    ! UCD export structure
    type(t_ucdExport) :: rexport
    real(DP), dimension(:,:), pointer :: p_Ddata2d
    real(DP), dimension(:), pointer :: Dtracer
    character(len=SYS_STRLEN) :: sucdfilename

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! Get filename and start GMV output
    call parlst_getvalue_string(rbloodflow%rparlist, 'Output', 'ucdfilename', sucdfilename)
    call ucd_startGMV(rexport, UCD_FLAG_STANDARD, rbloodflow%rtriangulation,&
        trim(sucdfilename)//'.'//trim(sys_si0(ifilenumber,5)))

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
    real(DP), intent(in) :: dtime

!</input>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DobjectCoords
    real(DP) :: a,b,c,L,t
    integer, dimension(2) :: Isize
    integer :: icase,ipoint, npoints

    ! Release thin object from previous evaluation
    if (rbloodflow%h_DobjectCoords .ne. ST_NOHANDLE) then
      call storage_free(rbloodflow%h_DobjectCoords)
    end if

    ! Get values from parameter list
    call parlst_getvalue_int(rbloodflow%rparlist, 'Object', 'npoints', npoints)
    call parlst_getvalue_int(rbloodflow%rparlist, 'Object', 'icase', icase)

    ! Generate list of vertices
    Isize = (/2, npoints/)
    call storage_new('bloodflow_evalObject', 'DobjectCoords', Isize,&
                     ST_DOUBLE, rbloodflow%h_DobjectCoords, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)
    
    select case(icase)
      
    case (1)
      ! Thin object
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'L', L)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'c', c)
      
      do ipoint = 0, npoints-1
        ! Compute x-coordinate
        p_DobjectCoords(1,ipoint+1) = L*ipoint/real(npoints-1, DP)
        
        ! Compute y-coordinate
        p_DobjectCoords(2,ipoint+1) = c*sin(dtime)*(3*L*p_DobjectCoords(1,ipoint+1)**2 -&
                                                        p_DobjectCoords(1,ipoint+1)**3)
      end do

    case (2)
      ! Rotating box
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'a', a)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'b', b)

      do ipoint = 0, npoints-1
        
        ! Compute parameter value
        t = 2*SYS_PI*ipoint/real(npoints-1, DP)

        ! Compute x-coordinate
        p_DobjectCoords(1,ipoint+1) = cos(dtime)*a*cos(t)-sin(dtime)*b*sin(t)

        ! Compute y-coordinate
        p_DobjectCoords(2,ipoint+1) = sin(dtime)*a*cos(t)+cos(dtime)*b*sin(t)
      end do
      
    case default
      print *, 'Invalid test case!'
      stop

    end select
    
  end subroutine bloodflow_evalObject

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_adaptObject(rbloodflow)

!<description>
    
    ! This subroutine adapts the computational grid to the object.

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iadapt,nadapt
    
    

    ! Initialize adaptation structure from triangulation
    call hadapt_initFromParameterlist(rbloodflow%rhadapt, rbloodflow%rparlist, 'Adaptation')
    call hadapt_initFromTriangulation(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    
    ! Perform prescribed number of standard h-adaptation steps
    call parlst_getvalue_int(rbloodflow%rparlist, 'Adaptation', 'nadapt', nadapt)    
    do iadapt = 1, nadapt
      
      ! Evaluate the indicator
      call bloodflow_evalIndicator(rbloodflow)
      
      ! Convert it to standard h-adaptation indicator
      call bloodflow_convertRefIndicator(rbloodflow)
      
      ! Adapt the computational mesh
      call hadapt_performAdaptation(rbloodflow%rhadapt, rbloodflow%rindicator)
      
      ! Generate raw mesh from adaptation structure
      call hadapt_generateRawMesh(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    end do


    do iadapt = 1, 1
      ! Evaluate the indicator
      call bloodflow_evalIndicator(rbloodflow)
      
      ! Perform rh-adaptation
      call bloodflow_performAdaptation(rbloodflow)
      
      ! Generate raw mesh from adaptation structure
      call hadapt_generateRawMesh(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    end do

    
    ! Release adaptation structure
    call hadapt_releaseAdaptation(rbloodflow%rhadapt)

  end subroutine bloodflow_adaptObject

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_evalIndicator(rbloodflow)

!<description>

    ! This subroutine evaluates the indicator function

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex, p_Iindicator
    integer, dimension(:), pointer :: IelementPatch
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Dpoint0Coords, Dpoint1Coords,Daux
    real(DP) :: dscale, dscaleMax
    integer :: ive,jve,iel,jel,i,i1,i2,i3,istatus,idx,ipoint,ipatch,npatch,ipos
  
    
    ! Release the indicator vector (if any)
    call lsyssc_releaseVector(rbloodflow%rindicator)

    ! Release the list of elements (if any)
    call list_releaseList(rbloodflow%relementList)
    
    ! Create new indicator vector as integer array
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true., ST_INT)
    call lsyssc_getbase_int(rbloodflow%rindicator, p_Iindicator)
    
    ! Allocate temporal memory for local element patch
    allocate(IelementPatch(rbloodflow%rtriangulation%NNelAtVertex))

    ! Create linked list for storing the elements adjacent to the object
    call list_createList(rbloodflow%relementList,&
        ceiling(0.1*rbloodflow%rtriangulation%NEL), ST_INT, 0, 0, 0)

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
      call PointInTriangleTest(DtriaCoords, Dpoint0Coords, POINT_EQUAL_TOLERANCE, istatus)
      if (istatus .ge. 0) exit findElementOfFirstVertex
      
    end do findElementOfFirstVertex

    ! Mark first element for potential refinement
    p_Iindicator(iel) = ibset(p_Iindicator(iel), istatus)

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
    !     vertex. In any case, we compute the coordinates of the
    !     intersection point which serves as new starting point of the
    !     segment. In practice, the original starting point is still
    !     used but the decision how to proceed is based on the
    !     coordinates of the intersection point. This process
    !     continues until the endpoint of the last segment is reached.
    !
    !---------------------------------------------------------------------------

    ! Loop over all remaining points
    findElementsOfOtherVertices: do while (ipoint .le. size(p_DobjectCoords, 2))
      
      ! Get global coordinates of the endpoint
      Dpoint1Coords = p_DobjectCoords(:,ipoint)

      ! Create single element patch
      npatch = 1
      IelementPatch(npatch) = iel

      ! Append element to the list of elements adjacent to the object
      if (iand(p_Iindicator(iel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
        p_Iindicator(iel) = ior(p_Iindicator(iel), BITFIELD_INLIST)
        call list_appendToList(rbloodflow%relementList, iel, ipos)
      end if
      
      ! Check status of starting point, aka, previous point
      select case(istatus)
      case (1:3)
        ! Previous point was one of the corner vertices
        i = p_IverticesAtElement(istatus, iel); ipatch = 0

        ! Create local patch surrounding this point
        do idx = p_IelementsAtVertexIdx(i),&
                 p_IelementsAtVertexIdx(i+1)-1
          
          ! Get global element number
          jel = p_IelementsAtVertex(idx)
          
          ! Append element to the list of elements adjacent to the object
          if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
            p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_INLIST)
            call list_appendToList(rbloodflow%relementList, jel, ipos)
          end if

          ! Mark element for potential refinement
          do jve = 1, TRIA_NVETRI2D
            if (p_IverticesAtElement(jve, jel) .eq. i) then
              p_Iindicator(jel) = ibset(p_Iindicator(jel), jve)
              exit
            end if
          end do
          
          ipatch = ipatch+1
          IelementPatch(ipatch) = jel

        end do
        npatch = ipatch
        
      case (4:6)               
        ! Previous point was located on the edge of the element
        jel = p_IneighboursAtElement(istatus-3, iel)

        ! Create two element patch       
        if (jel .ne. 0) then
          
          ! Append element to the list of elements adjacent to the object
          if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
            p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_INLIST)
            call list_appendToList(rbloodflow%relementList, jel, ipos)
          end if
          
          ! Mark element for potential refinement
          do jve = 1, TRIA_NVETRI2D
            if (p_IneighboursAtElement(jve, jel) .eq. iel) then

              ! Check if edge has been intersected previously
              if (btest(p_Iindicator(jel), jve+3))&
                  p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_MULTI_INTERSECTION)

              ! Mark edge as intersected
              p_Iindicator(jel) = ibset(p_Iindicator(jel), jve+3)
              exit
            end if
          end do

          npatch = 2
          IelementPatch(npatch) = jel

        end if
        
      end select
      
      ! Loop over all elements in patch and try to find 
      ! element which contains the endpoint
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
        call PointInTriangleTest(DtriaCoords, Dpoint1Coords, POINT_EQUAL_TOLERANCE, istatus)
        if (istatus .ge. 0) then
          
          ! Check if edge has been intersected previously
          if ((istatus .ge. 4) .and. btest(p_Iindicator(iel), istatus))&
              p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_MULTI_INTERSECTION)

          ! Mark element for potential refinement
          p_Iindicator(iel) = ibset(p_Iindicator(iel), istatus)

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

      ! Initialize scaling parameter
      dscaleMax = 0.0_DP

      ! Loop over all elements in patch and find the intersection
      ! point of the current line segment with the patch
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
        findIntersectionEdge: do ive = 1, TRIA_NVETRI2D
          
          ! Check if the edge intersects with the line between start
          ! and endpoint and compute coordinates of intersection point
          call LinesIntersectionTest(Dpoint0Coords, Dpoint1Coords, DtriaCoords(:,ive),&
                                     DtriaCoords(:,mod(ive,3)+1), istatus, dscale)
          if (istatus .eq. 0) then
            
            ! Check if the distance between the starting point and the
            ! intersection point exceeds that of an intersection point
            ! which may have been found previously for a different edge
            if (dscale .gt. dscaleMax) then
              ! Store parameter for new intersection point
              dscaleMax = dscale
              
              ! Store element number and edge position
              jel = iel; ipos = ive
            end if

          end if
          
        end do findIntersectionEdge
      end do findInPatchIndirect

      ! If we end up here, the scaling parameter for the "best"
      ! intersection point has been determined and we can proceed
      iel = jel

      ! Check if edge has been intersected previously
      if (btest(p_Iindicator(iel), ipos+3))&
          p_Iindicator(iel) = ior(p_Iindicator(iel), BITFIELD_MULTI_INTERSECTION)

      ! Mark element for potential refinement
      p_Iindicator(iel) = ibset(p_Iindicator(iel), ipos+3)
      
      ! Get vertices at element
      i1 = p_IverticesAtElement(1, iel)
      i2 = p_IverticesAtElement(2, iel)
      i3 = p_IverticesAtElement(3, iel)
      
      ! Get global coordinates of corner vertices
      DtriaCoords(:,1) = p_DvertexCoords(:,i1)
      DtriaCoords(:,2) = p_DvertexCoords(:,i2)
      DtriaCoords(:,3) = p_DvertexCoords(:,i3)
      
      ! Compute intersection point
      Daux = Dpoint0Coords + dscaleMax * (Dpoint1Coords-Dpoint0Coords)
        
      ! Perform point-in-triangle test for intersection point
      call PointInTriangleTest(DtriaCoords, Daux, POINT_EQUAL_TOLERANCE, istatus)
      
    end do findElementsOfOtherVertices
    

    ! Append element to the list of elements adjacent to the object
    if (iand(p_Iindicator(iel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
      p_Iindicator(iel) = ior(p_Iindicator(iel), BITFIELD_INLIST)
      call list_appendToList(rbloodflow%relementList, iel, ipos)
    end if
    
    ! Check status of previous point
    select case(istatus)
    case (1:3)
      ! Previous point was one of the corner vertices, thus, mark
      ! all elements meeting at that corner and create local patch
      i = p_IverticesAtElement(istatus, iel)
      
      do idx = p_IelementsAtVertexIdx(i),&
               p_IelementsAtVertexIdx(i+1)-1
        
        jel = p_IelementsAtVertex(idx)

        ! Append element to the list of elements adjacent to the object
        if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
          p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_INLIST)
          call list_appendToList(rbloodflow%relementList, jel, ipos)
        end if

        ! Mark element for potential refinement
        do jve = 1, TRIA_NVETRI2D
          if (p_IverticesAtElement(jve, jel) .eq. i) then
            p_Iindicator(jel) = ibset(p_Iindicator(jel), jve)
            exit
          end if
        end do

      end do
      
    case (4:6)
      ! Previous point was located on the edge of the element, thus,
      ! mark adjacent element and create two element patch
      jel = p_IneighboursAtElement(istatus-3, iel)
      
      ! Append element to the list of elements adjacent to the object
      if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
        p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_INLIST)
        call list_appendToList(rbloodflow%relementList, jel, ipos)
      end if

      ! Mark element for potential refinement
      do jve = 1, TRIA_NVETRI2D
        if (p_IneighboursAtElement(jve, jel) .eq. iel) then

          ! Check if edge has been intersected previously
          if (btest(p_Iindicator(jel), jve+3))&
              p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_MULTI_INTERSECTION)
          
          ! Mark edge as intersected
          p_Iindicator(jel) = ibset(p_Iindicator(jel), jve+3)
          exit
        end if
      end do
      
    end select

    ! Deallocate temporal memory
    deallocate(IelementPatch)

  end subroutine bloodflow_evalIndicator

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_convertRefIndicator(rbloodflow)

!<description>

    ! This subroutine converts the indicator function into a pure
    ! refinement indicator used in the h-adaptation procedure

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dindicator
    integer :: ipos, iel

    ! Release the indicator vector (if any)
    call lsyssc_releaseVector(rbloodflow%rindicator)

    ! Create new indicator vector as double
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Dindicator)
    
    !---------------------------------------------------------------------------
    ! Loop over all segments/elements and decide if the element needs
    ! refinement. This decision is simply based on the fact whether
    ! the element is intersected by the object or not.
    !---------------------------------------------------------------------------

    ! Loop over all elements adjacent to the object
    ipos = list_getNextInList(rbloodflow%relementList, .true.)
    list: do while(ipos .ne. LNULL)
      
      ! Get element number and proceed to next position
      call list_getByPosition(rbloodflow%relementList, ipos, iel)
      ipos = list_getNextInList(rbloodflow%relementList, .false.)
    
      ! Set indicator for element
      p_Dindicator(iel) = 1.0
    end do list

  end subroutine bloodflow_convertRefIndicator

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_performAdaptation(rbloodflow)

!<description>

    ! This subroutine performs rh-adaptation based on the indicator

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    integer, dimension(:), pointer :: p_Iindicator, p_Imarker
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Daux
    real(DP) :: dscale
    integer, dimension(2) :: Isize
    integer :: ive,jve,iel,i1,i2,i3,istatus,ipoint,ipos,iresult,nvt,nel


    ! Initialize marker structure
    if (rbloodflow%rhadapt%h_Imarker .ne. ST_NOHANDLE)&
        call storage_free(rbloodflow%rhadapt%h_Imarker)
    call storage_new('bloodflow_convertMoveRefIndicator', 'Imarker',&
        rbloodflow%rindicator%NEQ, ST_INT, rbloodflow%rhadapt%h_Imarker, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rbloodflow%rhadapt%h_Imarker, p_Imarker)

    ! Set pointers
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call storage_getbase_double2d(&
        rbloodflow%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)
    call lsyssc_getbase_int(rbloodflow%rindicator, p_Iindicator)


    ! Initialize initial dimensions
    call storage_getsize(rbloodflow%rhadapt%h_IverticesAtElement, Isize)
    rbloodflow%rhadapt%NELMAX = Isize(2)
    call storage_getsize(rbloodflow%rhadapt%h_IneighboursAtElement, Isize)
    rbloodflow%rhadapt%NELMAX = min(rbloodflow%rhadapt%NELMAX, Isize(2))

    rbloodflow%rhadapt%InelOfType0 = rbloodflow%rhadapt%InelOfType
    rbloodflow%rhadapt%NVT0        = rbloodflow%rhadapt%NVT
    rbloodflow%rhadapt%NEL0        = rbloodflow%rhadapt%NEL
    rbloodflow%rhadapt%NVBD0       = rbloodflow%rhadapt%NVBD
    rbloodflow%rhadapt%increaseNVT = 0
    

    !---------------------------------------------------------------------------
    ! (1) Loop over all segments/elements and decide whether to refine or to
    !     reposition mesh points so that elements get aligned to the object
    !---------------------------------------------------------------------------

    ! Loop over all elements adjacent to the object
    ipos = list_getNextInList(rbloodflow%relementList, .true.)
    list: do while(ipos .ne. LNULL)
      
      ! Get element number and proceed to next position
      call list_getByPosition(rbloodflow%relementList, ipos, iel)
      ipos = list_getNextInList(rbloodflow%relementList, .false.)
    
      ! Check if multi-intersections are present
      if (iand(p_Iindicator(iel), BITFIELD_MULTI_INTERSECTION) .eq.&
                                  BITFIELD_MULTI_INTERSECTION) then

        ! Set the refinement indicator to unity
        p_Imarker(iel) = 14
        print *, "Multiple intersection"

        ! Compute number of new vertices
        do ive = 1, TRIA_NVETRI2D
          if (p_IneighboursAtElement(ive, iel) .eq. 0) then
            
            ! Edge is adjacent to boundary
            rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1

          elseif(p_Imarker(p_IneighboursAtElement(ive, iel)) .eq. 0) then

            ! Edge has not been marked in previous steps
            rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1

          end if
        end do
        
        ! That's it
        cycle list
      end if
  
      ! Remove "in list" flag from indicator
      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_INLIST))
      
      ! Remove "corner vertices" flags from indicator.
      ! This may be used in future versions of this code.
      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_POINT_INTERSECTION))
      
      ! Remove "interior vertex" flag from indicator.
      ! This may be used in future versions of this code.
      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_INNER))


      ! Check status of intersected element edges
      select case(iand(p_Iindicator(iel), BITFIELD_EDGE_INTERSECTION))
      case (0)
        ! No edge is intersected, hence, do nothing
        
      case (BITFIELD_EDGE1 + BITFIELD_EDGE2 + BITFIELD_EDGE3)
        ! All three edges are intersected!
        ! Set the refinement indicator to unity
        p_Imarker(iel) = 14
        print *, "All three edges intersected"

        ! Compute number of new vertices
        do ive = 1, TRIA_NVETRI2D
          if (p_IneighboursAtElement(ive, iel) .eq. 0) then
            
            ! Edge is adjacent to boundary
            rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1

          elseif(p_Imarker(p_IneighboursAtElement(ive, iel)) .eq. 0) then

            ! Edge has not been marked in previous steps
            rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1

          end if
        end do
        
      case default
        
        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        ! Process each edge individually
        edge: do ive = 1, TRIA_NVETRI2D

          ! Check if ive-th edge is intersected  
          if (btest(p_Iindicator(iel), ive+3)) then

            ! Loop over all segments and check intersection point
            point: do ipoint = 1, size(p_DobjectCoords,2)-1
              
              ! Perform line intersection test
              call LinesIntersectionTest(p_DobjectCoords(:,ipoint), p_DobjectCoords(:,ipoint+1),&
                                         DtriaCoords(:,ive), DtriaCoords(:,mod(ive,3)+1), istatus, dscale)
              if (istatus .eq. 1) then
                
                ! Remove the refinement indicator for current edge
                p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+3)
                
              elseif (istatus .eq. 0) then
                
                ! Compute coordinates of intersection point
                Daux = p_DobjectCoords(:,ipoint) + dscale *&
                      (p_DobjectCoords(:,ipoint+1)-p_DobjectCoords(:,ipoint))
                
                ! Perform point in triangle test
                call PointInTriangleTest(DtriaCoords, Daux, POINT_COLLAPSE_TOLERANCE, istatus)
                select case(istatus)
                case(1)
                  ! Update coordinate of first vertex in quadtree
                  iresult = qtree_moveInQuadtree(rbloodflow%rhadapt%rVertexCoordinates2D,&
                                                 p_DvertexCoords(:, i1), Daux)
                  
                  ! Adjust first vertex
                  p_DvertexCoords(:, i1) = Daux
                  
                  ! Remove the refinement indicator for current edge
                  p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+3)

                  
                case(2)
                  ! Update coordinate of second vertex in quadtree
                  iresult = qtree_moveInQuadtree(rbloodflow%rhadapt%rVertexCoordinates2D,&
                                                 p_DvertexCoords(:, i2), Daux)

                  ! Adjust second vertex
                  p_DvertexCoords(:, i2) = Daux
                  
                  ! Remove the refinement indicator for current edge
                  p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+3)


                case(3)
                  ! Update coordinate of third vertex in quadtree
                  iresult = qtree_moveInQuadtree(rbloodflow%rhadapt%rVertexCoordinates2D,&
                                                 p_DvertexCoords(:, i3), Daux)

                  ! Adjust third vertex
                  p_DvertexCoords(:, i3) = Daux
                  
                  ! Remove the refinement indicator for current edge
                  p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+3)

                case default
                  p_Imarker(iel) = ibset(p_Imarker(iel), ive)

                  ! Compute number of new vertices
                  do jve = 1, TRIA_NVETRI2D
                    if (p_IneighboursAtElement(jve, iel) .eq. 0) then
                      
                      ! Edge is adjacent to boundary
                      rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1
                      
                    elseif(p_Imarker(p_IneighboursAtElement(jve, iel)) .eq. 0) then
                      
                      ! Edge has not been marked in previous steps
                      rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1
                      
                    end if
                  end do

                  ! Keep this edge
                  cycle edge
                end select
             
              end if
              
            end do point

            ! Marker for this edge can be removed
            p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+3)

          end if
        end do edge

      end select
    end do list
    
    !---------------------------------------------------------------------------
    ! (2) Perform h-adaptation based on the generated marker array
    !---------------------------------------------------------------------------

    ! Set specifier to "marked for refinement"
    rbloodflow%rhadapt%iSpec = ior(rbloodflow%rhadapt%iSpec, HADAPT_MARKEDREFINE)

    ! Mark additional elements to restore conformity
    call hadapt_markRedgreenRefinement2D(rbloodflow%rhadapt)

    ! Since no coarsening is performed also
    ! set specifier to "marked for coarsening"
    rbloodflow%rhadapt%iSpec = ior(rbloodflow%rhadapt%iSpec, HADAPT_MARKEDCOARSEN)

    
    ! Compute new dimensions
    nvt = rbloodflow%rhadapt%NVT+rbloodflow%rhadapt%increaseNVT
    nel = hadapt_CalcNumberOfElements2D(rbloodflow%rhadapt)

    ! Adjust array IvertexAge
    call storage_realloc('hadapt_performAdaptation', nvt,&
                         rbloodflow%rhadapt%h_IvertexAge, ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int(rbloodflow%rhadapt%h_IvertexAge, rbloodflow%rhadapt%p_IvertexAge)

    ! Adjust array InodalProperty
    if (iand(rbloodflow%rhadapt%iSpec, HADAPT_HAS_NODALPROP) .eq.&
                                       HADAPT_HAS_NODALPROP) then
      call storage_realloc('hadapt_performAdaptation', nvt,&
                           rbloodflow%rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int(rbloodflow%rhadapt%h_InodalProperty, rbloodflow%rhadapt%p_InodalProperty)
    end if
    
    ! Adjust array IverticesAtElement
    if (iand(rbloodflow%rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .eq.&
                                       HADAPT_HAS_VERTATELEM) then
      call storage_realloc('hadapt_performAdaptation', nel,&
                           rbloodflow%rhadapt%h_IverticesAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int2D(rbloodflow%rhadapt%h_IverticesAtElement,&
                                 rbloodflow%rhadapt%p_IverticesAtElement)
    end if

    ! Adjust array IneighboursAtElement
    if (iand(rbloodflow%rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM) .eq.&
                                       HADAPT_HAS_NEIGHATELEM) then
      call storage_realloc('hadapt_performAdaptation', nel,&
                           rbloodflow%rhadapt%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int2D(rbloodflow%rhadapt%h_IneighboursAtElement,&
                                 rbloodflow%rhadapt%p_IneighboursAtElement)
    end if

    ! Adjust array ImidneighboursAtElement
    if (iand(rbloodflow%rhadapt%iSpec, HADAPT_HAS_MIDNEIGH) .eq.&
                                       HADAPT_HAS_MIDNEIGH) then
      call storage_realloc('hadapt_performAdaptation', nel,&
                           rbloodflow%rhadapt%h_ImidneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int2D(rbloodflow%rhadapt%h_ImidneighboursAtElement,&
                                 rbloodflow%rhadapt%p_ImidneighboursAtElement)
    end if

    ! Perform element refinement in 2D
    call hadapt_refine2D(rbloodflow%rhadapt)
    
    ! Adjust nodal property array
    call storage_realloc('hadapt_performAdaptation', rbloodflow%rhadapt%NVT,&
                         rbloodflow%rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)

  end subroutine bloodflow_performAdaptation
  
  !*****************************************************************************

!<subroutine>
  
  pure subroutine PointInTriangleTest(TriaCoords, P, dtolerance, istatus)

!<description>
    
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

!</description>
    
!<input>
    
    ! Coordinates of the triangle
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(in) :: TriaCoords

    ! Coordinates of the point
    real(DP), dimension(NDIM2D), intent(in) :: P

    ! Tolerance for point collapse
    real(DP), intent(in) :: dtolerance

!</input>

!<output>

    ! Status of the test
    integer, intent(out) :: istatus

!</output>

!</subroutine>
        
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
    if ((u .lt. -dtolerance) .or. (v .lt. -dtolerance) .or.&
        (w .lt. -dtolerance) .or. (u+v+w .gt. 1+dtolerance)) then
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
  
  !*****************************************************************************

!<subroutine>
  
  subroutine LinesIntersectionTest(P1, P2, P3, P4, istatus, u)

!<description>
    
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
    !
    ! u: is the parameter so that 
    !    S = P1 + u*(P2-P1) is the intersection point

!</description>
    
!<input>

    ! Coordinates of the four points
    real(DP), dimension(NDIM2D), intent(in) :: P1,P2,P3,P4

!</input>

!<output>
    
    ! Parameter of the intersection point w.r.t. the first line
    real(DP), intent(out) :: u

    ! Status of the test
    integer, intent(out) :: istatus

!</output>

!</subroutine>
    
    ! local variables
    real(DP) :: denom,nom1,nom2,v
    
    
    ! Compute (de-)nominators
    denom = (P4(2)-P3(2)) * (P2(1)-P1(1)) - (P4(1)-P3(1)) * (P2(2)-P1(2))
    nom1  = (P4(1)-P3(1)) * (P1(2)-P3(2)) - (P4(2)-P3(2)) * (P1(1)-P3(1))
    nom2  = (P2(1)-P1(1)) * (P1(2)-P3(2)) - (P2(2)-P1(2)) * (P1(1)-P3(1))
    
    ! Check if lines are parallel
    if (abs(denom) .le. SYS_EPSREAL) then
      
      ! The two lines are parallel, hence there is no intersection point
      u = SYS_INFINITY
      
      ! Check of both lines coincide are not
      if ( (abs(nom1) .le. SYS_EPSREAL) .or.&
           (abs(nom2) .le. SYS_EPSREAL) ) then
        istatus =  1
      else
        istatus = -1
      end if
      
    else
      
      ! The two lines are not parallel, hence they must intersect each other
      
      ! Compute parameter values
      u = nom1 / denom
      v = nom2 / denom
      
      ! Check if the intersection point is located within the line segments
      if ((u .lt. 0) .or. (u .gt. 1) .or.&
          (v .lt. 0) .or. (v .gt. 1)) then
        
        ! Intersection point does not belong to both line segments
        istatus = -1          
        
      else
        
        ! Intersection point belongs to both line segments
        istatus = 0
        
      end if
      
    end if
    
  end subroutine LinesIntersectionTest
  
  !*****************************************************************************

!<subroutine>
  
  pure subroutine getBarycentricCoords(TriaCoords, P, BarycentricCoords)

!<description>
    
    ! This subroutine calculates the barycentric coordinates of the
    ! given point P using the given reference triangle.

!</description>
   
!<input>

    ! Coordinates of the triangle
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(in) :: TriaCoords

    ! Coordinates of the point
    real(DP), dimension(NDIM2D), intent(in) :: P

!</input>

!<output>
    
    ! Barycentric coordinates of the point
    real(DP), dimension(TRIA_NVETRI2D), intent(out) :: BarycentricCoords

!</output>

!</subroutine>
        
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
  
  !*****************************************************************************

!<function>
  
  pure function signedArea(P1,P2,P3) result(area)
    
!<description>

    ! This function computes the signed area of the triangle

!</description>
   
!<input>

    ! Coordinates of the triangle
    real(DP), dimension(NDIM2D), intent(in) :: P1,P2,P3

!</input>

!<result>
    
    ! Area of the triangle
    real(DP) :: area

!</result>

!</function>
    
    area = 0.5 * ( (P2(1)-P1(1))*(P3(2)-P1(2)) -&
        (P3(1)-P1(1))*(P2(2)-P1(2)) )
    
  end function signedArea

end module bloodflow

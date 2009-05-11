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

  use bloodflow_msd
  use boundary
  use genoutput
  use geometry
  use linearsystemscalar
  use paramlist
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
  public :: bloodflow_evalIndicator

  !*****************************************************************************

!<constants>

!<constantblock description="Global constants for mesh spacing tolerances">

  
  real(DP), parameter :: POINT_COLLAPSE_TOLERANCE = 1e-4_DP
  
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

    ! Bounding box which restricts the movement of vertices
    type(t_geometryObject) :: rboundingBox

    ! Handle to array storing the points of the object
    integer :: h_DobjectCoords = ST_NOHANDLE

    ! Indicator vector
    type(t_vectorScalar) :: rindicator
    
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
    ! - creates the bounding box which restricts the movement of nodes

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

    ! Generate bounding box
    call parlst_getvalue_double(rbloodflow%rparlist, 'BoundingBox', 'xorigin', Dorigin(1))
    call parlst_getvalue_double(rbloodflow%rparlist, 'BoundingBox', 'yorigin', Dorigin(2))
    call parlst_getvalue_double(rbloodflow%rparlist, 'BoundingBox', 'xlength', Dlength(1))
    call parlst_getvalue_double(rbloodflow%rparlist, 'BoundingBox', 'ylength', Dlength(2))
    call geom_init_rectangle(rbloodflow%rboundingBox, Dlength, Dorigin)   
    
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

      ! Release the bounding box
      call geom_done(rbloodflow%rboundingBox)
      
      ! Release the object coordinates
      call storage_free(rbloodflow%h_DobjectCoords)

      ! Release indicator vector
      call lsyssc_releaseVector(rbloodflow%rindicator)

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
    call ucd_startGMV(rexport, UCD_FLAG_STANDARD, rbloodflow%rtriangulation, trim(sucdfilename)//'.'//trim(sys_si0(ifilenumber,5)))

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
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dindicator
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
    integer, dimension(:), pointer :: IelementPatch
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Dpoint0Coords, Dpoint1Coords,DIntersectionCoords
    real(DP) :: daux
    integer :: ive,jve,iel,jel,i,i1,i2,i3,istatus,idx,ipoint,ipatch,npatch
    
    
    ! Release the indicator vector
    call lsyssc_releaseVector(rbloodflow%rindicator)

    ! Create new indicator vector
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Dindicator)

    ! Allocate temporal memory for local element patch
    allocate(IelementPatch(rbloodflow%rtriangulation%NNelAtVertex))

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
      call PointInTriangleTest(DtriaCoords, Dpoint0Coords, istatus)
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
      
      ! Mark element
      p_Dindicator(iel) = 1

      ! Create single element patch
      npatch = 1
      IelementPatch(1) = iel

      ! Check status of previous point
      select case(istatus)
        
      case (1:3)
        ! Previous point was one of the corner vertices, thus, mark
        ! all elements meeting at that corner and create local patch
        i = p_IverticesAtElement(istatus, iel); ipatch = 0
        
        do idx = p_IelementsAtVertexIdx(i), p_IelementsAtVertexIdx(i+1)-1
          
          jel = p_IelementsAtVertex(idx)
          p_Dindicator(jel) = 1
          
          ipatch = ipatch+1
          IelementPatch(ipatch) = jel
        end do
        npatch = ipatch
        
      case (4:6)
        ! Previous point was located on the edge of the element, thus,
        ! mark adjacent element and create two element patch
        jel = p_IneighboursAtElement(istatus-3, iel)

        if (jel .ne. 0) then
          p_Dindicator(jel) = 1
          
          npatch = 2
          IelementPatch(2) = jel
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
        call PointInTriangleTest(DtriaCoords, Dpoint1Coords, istatus)
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
            call PointInTriangleTest(DtriaCoords, Dpoint1Coords, istatus)
            
            if (istatus .eq. -1) then
              
              ! If the point is 'outside' the element, then the
              ! segment just passes through the element, whereby no
              ! point falls into it. In this case, the intersection
              ! point serves as new starting point of the segment
              Dpoint0Coords = DintersectionCoords

              ! Perform point-in-triangle test for intersection point.
              call PointInTriangleTest(DtriaCoords, Dpoint0Coords, istatus)
              
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
    
    ! Mark last element
    p_Dindicator(iel) = 1
    
    ! Check status of previous point
    select case(istatus)
      
    case (1:3)
      ! Previous point was one of the corner vertices, thus, mark
      ! all elements meeting at that corner and create local patch
      i = p_IverticesAtElement(istatus, iel)
      
      do idx = p_IelementsAtVertexIdx(i), p_IelementsAtVertexIdx(i+1)-1
        jel = p_IelementsAtVertex(idx)
        p_Dindicator(jel) = 1
      end do
      
    case (4:6)
      ! Previous point was located on the edge of the element, thus,
      ! mark adjacent element and create two element patch
      jel = p_IneighboursAtElement(istatus-3, iel)
      
      if (jel .ne. 0) p_Dindicator(jel) = 1
    end select

    ! Deallocate temporal memory
    deallocate(IelementPatch)


    !---------------------------------------------------------------------------
    ! (3) Find nearest neighbor of all points of the object:
    !
    !---------------------------------------------------------------------------
    
    

  contains
    
    !***************************************************************************
    
    pure subroutine PointInTriangleTest(TriaCoords, P, istatus)

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

      integer, intent(OUT) :: istatus

      ! local parameters
      real(DP), parameter :: PERCENTAGE_DISTANCE_TOLERACNE = 1e-3_DP

      ! local variables
      real(DP), dimension(NDIM2D) :: v0,v1,v2
      real(DP) :: dot00, dot01, dot02, dot11, dot12, invDenom, u, v

      ! Compute vectors
      v0 = TriaCoords(:,2) - TriaCoords(:,1)
      v1 = TriaCoords(:,3) - TriaCoords(:,1)
      v2 = P               - TriaCoords(:,1)
      
      ! Compute dot products
      dot00 = sum(v0*v0)
      dot01 = sum(v0*v1)
      dot02 = sum(v0*v2)
      dot11 = sum(v1*v1)
      dot12 = sum(v1*v2)
      
      ! Compute barycentric coordinates
      invDenom = 1.0_DP / (dot00 * dot11 - dot01 * dot01)
      u = (dot11 * dot02 - dot01 * dot12) * invDenom
      v = (dot00 * dot12 - dot01 * dot02) * invDenom
      
      ! Determine status
      if ((u .lt. 0) .or. (v .lt. 0) .or. (u+v .gt. 1)) then
        ! If the barycentric coordinates are negative of larger than
        ! one then the point is located outside of the triangle
        istatus = -1

        ! Otherwise, the point is located inside the triangle. We have
        ! to considere several cases, e.g., the point is close to an
        ! edge or it is close to a vertex and collapses
      elseif (u .le. POINT_COLLAPSE_TOLERANCE) then
        if (v .le. POINT_COLLAPSE_TOLERANCE) then
          ! Point P collapses with point A
          istatus = 1
        elseif (u+v .ge. 1-POINT_COLLAPSE_TOLERANCE) then
          ! Point P collapses with point C
          istatus = 3
        else
          ! Point P collapses with edge AC
          istatus = 6
        end if
      elseif (v .le. POINT_COLLAPSE_TOLERANCE) then
        if (u+v .ge. 1-POINT_COLLAPSE_TOLERANCE) then
          ! Point P collapses with point B
          istatus = 2
        else
          ! Point P collapses with edge AB
          istatus = 4
        end if
      elseif (u+v .ge. 1-POINT_COLLAPSE_TOLERANCE) then
        ! Point P collapses with edge BC
        istatus =5
      else
        ! Point P is in the interior
        istatus = 0
      end if
      
    end subroutine PointInTriangleTest

    !***************************************************************************

    subroutine LinesIntersectionTest(P1, P2, P3, P4, istatus, S)

      ! This subroutine checks if the two line segments P1-P2 and
      ! P3-P4 intersect each other and returns the intersection point.
      ! Note that the start point is excluded from the
      ! segment. Otherwise, the starting point would by considered as
      ! intersection point if it was located on the edge.
      ! The meaning of the resulting istatus is as follows:
      !
      ! istatus:
      !  = -1 : the two line segments do not intersect
      !  =  0 : the two line segments intersect each other
      
      real(DP), dimension(NDIM2D), intent(IN) :: P1,P2,P3,P4

      real(DP), dimension(NDIM2D), intent(OUT) :: S
      integer, intent(OUT) :: istatus

      ! local variables
      real(DP) :: aux,u,v


      ! Compute denominator
      aux = (P4(2)-P3(2)) * (P2(1)-P1(1)) - (P4(1)-P3(1)) * (P2(2)-P1(2))

      ! Check if lines are parallel
      if (aux .eq. 0) then
        istatus = -1;   S = SYS_INFINITY
        return
      end if
      
      ! Compute paremters
      u = ((P4(1)-P3(1)) * (P1(2)-P3(2)) - (P4(2)-P3(2)) * (P1(1)-P3(1))) / aux
      v = ((P2(1)-P1(1)) * (P1(2)-P3(2)) - (P2(2)-P1(2)) * (P1(1)-P3(1))) / aux
      
      ! Determine status
      if ((u .le. sqrt(SYS_EPSREAL)) .or. (u .gt. 1) .or. (v .lt. 0) .or. (v .gt. 1)) then

        istatus = -1

        S = SYS_INFINITY

      else
        
        istatus = 0
        
        ! Compute intersection point
        S(1) = P1(1) + u*(P2(1)-P1(1))
        S(2) = P1(2) + u*(P2(2)-P1(2))
      end if

    end subroutine LinesIntersectionTest

  end subroutine bloodflow_evalIndicator

end module bloodflow

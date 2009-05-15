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

!!$  use bloodflow_msd
  use bilinearformevaluation
  use boundary
  use genoutput
  use geometry
  use linearalgebra
  use linearsolver
  use linearsystemscalar
  use linearsystemblock
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
  public :: bloodflow_evalIndicator
  public :: bloodflow_redistMeshPoints

  !*****************************************************************************

!<constants>

!<constantblock description="Global constants for mesh spacing tolerances">

  ! Tolerance for considering two points as equivalent
  real(DP), parameter :: POINT_COLLAPSE_TOLERANCE = 1e-4_DP

  ! Stiffness parameter for mass spring system
  real(DP), parameter :: SPRING_STIFFNESS = 1._DP
  
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
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement, p_Isprings
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
    integer, dimension(:), pointer :: IelementPatch
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Dpoint0Coords, Dpoint1Coords,DIntersectionCoords
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

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_redistMeshPoints(rbloodflow)

!<description>

    ! This subroutine redistributes the mesh points based on the
    ! equilibration of a two-dimensional mass spring system

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(INOUT) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    type(t_blockDiscretisation) :: rdiscretisation
    type(t_matrixBlock), dimension(1) :: Rmatrix
    type(t_matrixScalar) :: rmatrixScalar
    type(t_vectorBlock) :: rparticles, rforces, rincrement, rtemp
    type(t_linsolNode), pointer :: p_rsolverNode
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dmatrix, p_Dparticles, p_Dforces
    real(DP), dimension(1) :: dnorm1, dnorm2
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, ierror, ite


    ! Create matrix structure for standard P1 discretization
    call spdiscr_initBlockDiscr(rdiscretisation, 1,&
        rbloodflow%rtriangulation, rbloodflow%rboundary)
    call spdiscr_initDiscr_simple(rdiscretisation%RspatialDiscr(1), &
        EL_E001, SPDISC_CUB_AUTOMATIC, rbloodflow%rtriangulation, rbloodflow%rboundary)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9, rmatrixScalar)
    
    ! Create scalar matrix with local 2x2 blocks
    rmatrixScalar%NVAR = 2
    rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9INTL
    rmatrixScalar%cinterleavematrixFormat = LSYSSC_MATRIX1
    call lsyssc_allocEmptyMatrix(rmatrixScalar, LSYSSC_SETM_UNDEFINED)

    ! Create 1-block system matrix and set points
    call lsysbl_createMatFromScalar(rmatrixScalar, Rmatrix(1), rdiscretisation)    
    call lsyssc_getbase_double(rmatrixScalar, p_Dmatrix)
    call lsyssc_getbase_Kld(rmatrixScalar, p_Kld)
    call lsyssc_getbase_Kcol(rmatrixScalar, p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmatrixScalar, p_Kdiagonal)
    
    ! Create diagonal separator
    h_Ksep = ST_NOHANDLE
    call storage_copy(rmatrixScalar%h_Kld, h_Ksep)
    call storage_getbase_int(h_Ksep, p_Ksep, rmatrixScalar%NEQ+1)
    
    ! Create vectors
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rparticles, .true.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rforces, .true.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rincrement, .false.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rtemp, .false.)
    call lsysbl_getbase_double(rparticles, p_Dparticles)
    call lsysbl_getbase_double(rforces, p_Dforces)

    ! Set pointers to coordinates vectors and triangulation data
    call storage_getbase_double2d(&
        rbloodflow%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)

    ! Create a linear BiCGSTAB solver
    call linsol_initBiCGStab(p_rsolverNode)
    p_rsolverNode%ioutputLevel       = 1
    p_rsolverNode%istoppingCriterion = LINSOL_STOP_ONEOF
    p_rsolverNode%depsRel            = 1e-3
    p_rsolverNode%depsAbs            = 1e-8
    
    ! Attach system matrix and initialize structures and data
    call linsol_setMatrices(p_RsolverNode, Rmatrix)
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop

    ! Initialize particle positions by point coordinates
    call setParticles(p_DvertexCoords, rbloodflow%rtriangulation%NVT, p_Dparticles)

!!$    ! Move point by hand
!!$    p_Dparticles(2*1013+1) = p_Dparticles(2*1013+1) + 0.05
!!$    p_Dparticles(2*1013+2) = p_Dparticles(2*1013+2) + 0.05
    
    ! Perform nonlinear iterations
    newton: do ite = 1, 100

      ! Clear matrix and force vector
      call lalg_clearVector(p_Dmatrix)
      call lalg_clearVector(p_Dforces)

      ! Assemble matrix and force vector
      call assembleForces(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_IverticesAtElement,&
          rmatrixScalar%NA, rmatrixScalar%NEQ, rbloodflow%rtriangulation%NEL, &
          p_Dmatrix, p_Dparticles, p_Dforces, p_DvertexCoords, p_DobjectCoords)
      
      ! Reset diagoanl separator
      call lalg_copyVector(p_Kld, p_Ksep)

      ! Set fixed points
      call setFixedPoints(p_Kld, p_Kcol, p_Kdiagonal, rmatrixScalar%NA,&
          rmatrixScalar%NEQ, p_Dmatrix, p_Dparticles, p_Dforces)
      
      ! Solve linear problem
      call lsysbl_clearVector(rincrement)
      call linsol_solveAdaptively (p_rsolverNode, rincrement, rforces, rtemp)
      
      ! Update particle positions
      call lsysbl_vectorLinearComb(rincrement, rparticles, -1.0_DP, 1.0_DP)

      ! Compute relative changes
      call lsysbl_vectorNormBlock(rincrement, (/LINALG_NORMEUCLID/), Dnorm1)
      call lsysbl_vectorNormBlock(rparticles, (/LINALG_NORMEUCLID/), Dnorm2)

      print *, Dnorm1/Dnorm2, p_rsolverNode%dfinalDefect, p_rsolverNode%iiterations, p_rsolverNode%dinitialDefect

      ! Check relative changes and exit Newton algorithm
!      if (Dnorm1(1)/Dnorm2(1) .le. 1e-6) exit newton
    end do newton

    ! Update point coordinates
    call setVertexCoords(rbloodflow%rtriangulation%NVT, p_Dparticles, p_DvertexCoords)

    ! Release temporal memory
    call storage_free(h_Ksep)
    call lsysbl_releaseVector(rparticles)
    call lsysbl_releaseVector(rforces)
    call lsysbl_releaseVector(rincrement)
    call lsysbl_releaseVector(rtemp)
    call lsysbl_releaseMatrix(Rmatrix(1))
    call lsyssc_releaseMatrix(rmatrixScalar)
    call spdiscr_releaseBlockDiscr(rdiscretisation, .true.)
    call linsol_doneData(p_rsolverNode)
    call linsol_doneStructure(p_rsolverNode)
    call linsol_releaseSolver(p_rsolverNode)

  contains

    !***************************************************************************

    subroutine setParticles(DvertexCoords, n, Dparticles)

      ! This subroutine initializes the particle positions by the
      ! physical coordinates of the grid points

      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      integer, intent(IN) :: n
      real(DP), dimension(NDIM2D, n), intent(OUT) :: Dparticles

      
      ! Copy data
      call lalg_copyVector(DvertexCoords, Dparticles)
      
    end subroutine setParticles

    !***************************************************************************

    subroutine setVertexCoords(n, Dparticles, DvertexCoords)

      ! This subroutine updates the physical coordinates of the grid
      ! points by the particle positions

      integer, intent(IN) :: n
      real(DP), dimension(NDIM2D, n), intent(IN) :: Dparticles
      real(DP), dimension(:,:), intent(INOUT) :: DvertexCoords
      
      ! Copy data
      call lalg_copyVector(Dparticles, DvertexCoords)
      
    end subroutine setVertexCoords
    
    !***************************************************************************
    
    subroutine assembleForces(Kld, Kcol, Kdiagonal, Ksep, IverticesAtElement,&
        na, neq, nel, Dmatrix, Dparticles, Dforces, DvertexCoords, DobjectCoords)

      ! This subroutine assembles the Jacobian matrix of the
      ! spring-mass system and the force vector which contains the
      ! forces for the linear edge springs satisfying Hooke's law.
      ! The force vector contains both the contributions from the edge
      ! springs and from the non-linear projection springs which are
      ! required to prevent collapsing of elements.
      
      real(DP), dimension(NDIM2D,neq), intent(IN) :: Dparticles
      real(DP), dimension(:,:), intent(IN) :: DvertexCoords, DobjectCoords
      integer, dimension(:,:), intent(IN) :: IverticesAtElement
      integer, dimension(:), intent(IN) :: Kld, Kcol, Kdiagonal
      integer, intent(IN) :: na, neq, nel
      
      real(DP), dimension(NDIM2D,NDIM2D,na), intent(INOUT) :: Dmatrix
      real(DP), dimension(NDIM2D,neq), intent(INOUT) :: Dforces
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D,NDIM2D) :: J_ij
      real(DP), dimension(NDIM2D) :: f_ij, d_ij, Dproj
      real(DP) :: dlength, dlength0, dprj, d2_ij, d3_ij
      integer :: i,j,k,ii,ij,ji,jj,iel,ive

      
      ! Loop over all rows
      rows: do i = 1, neq
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Loop over all off-diagonal matrix entries such that i<j
        cols: do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number j, the corresponding matrix positions ji
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute spring length and rest length
          dlength  = sqrt((Dparticles(1,i)-Dparticles(1,j))**2 +&
                          (Dparticles(2,i)-Dparticles(2,j))**2)
          dlength0 = sqrt((DvertexCoords(1,i)-DvertexCoords(1,j))**2 +&
                          (DvertexCoords(2,i)-DvertexCoords(2,j))**2)
          
          ! Compute coefficients and direction
          d_ij  = (Dparticles(:,j)-Dparticles(:,i))
          d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)
          d3_ij = d2_ij*sqrt(d2_ij)

          ! Compute local Jacobian matrix
          J_ij(1,1) = d2_ij - d_ij(1)*d_ij(1)
          J_ij(2,1) =       - d_ij(2)*d_ij(1)
          J_ij(1,2) =       - d_ij(1)*d_ij(2)
          J_ij(2,2) = d2_ij - d_ij(2)*d_ij(2)

          ! Scale local Jacobian matrix
          J_ij = J_ij / d3_ij

          ! Update local Jacobian matrix
          J_ij(1,1) = SPRING_STIFFNESS - SPRING_STIFFNESS * dlength0 * J_ij(1,1)
          J_ij(1,2) =                  - SPRING_STIFFNESS * dlength0 * J_ij(1,2)
          J_ij(2,1) =                  - SPRING_STIFFNESS * dlength0 * J_ij(2,1)
          J_ij(2,2) = SPRING_STIFFNESS - SPRING_STIFFNESS * dlength0 * J_ij(2,2)

          ! Apply off-diagonal matrix entries
          Dmatrix(:,:,ij) = J_ij
          Dmatrix(:,:,ji) = J_ij

          ! Update diagonal matrix entris
          Dmatrix(:,:,ii) = Dmatrix(:,:,ii) - J_ij
          Dmatrix(:,:,jj) = Dmatrix(:,:,jj) - J_ij
          
          ! Compute force vector by Hooke's law for edge (i,j)
          f_ij =  SPRING_STIFFNESS * (dlength-dlength0) * d_ij / dlength
    
          ! Apply force vector to nodes i and j
          Dforces(:,i) = Dforces(:,i) + f_ij
          Dforces(:,j) = Dforces(:,j) - f_ij
          
        end do cols
      end do rows
      
      i = 895

      d_ij = DobjectCoords(:,100) - Dparticles(:,i)
      d2_ij = sqrt(d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2))
      d3_ij = d2_ij**3

      ! Compute forces
      Dforces(:,i) = Dforces(:,i) + d_ij / d3_ij

      ! Compute local Jacobian matrix
      J_ij(1,1) = d_ij(1)*d_ij(1)
      J_ij(2,1) = d_ij(2)*d_ij(1)
      J_ij(1,2) = d_ij(1)*d_ij(2)
      J_ij(2,2) = d_ij(2)*d_ij(2)

      ! Scale local Jacobian matrix
      J_ij = 3 * J_ij / d2_ij**5

      ! Update local Jacobian matrix
      J_ij(1,1) = d2_ij**2 - J_ij(1,1)
      J_ij(1,2) =          - J_ij(1,2)
      J_ij(2,1) =          - J_ij(2,1)
      J_ij(2,2) = d2_ij**2 - J_ij(2,2)
      
      ii = Kdiagonal(i)
      Dmatrix(:,:,ii) = Dmatrix(:,:,ii) + J_ij



      d_ij = DobjectCoords(:,99) - Dparticles(:,i)
      d2_ij = sqrt(d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2))
      d3_ij = d2_ij**3

      ! Compute forces
      Dforces(:,i) = Dforces(:,i) + d_ij / d3_ij

      ! Compute local Jacobian matrix
      J_ij(1,1) = d_ij(1)*d_ij(1)
      J_ij(2,1) = d_ij(2)*d_ij(1)
      J_ij(1,2) = d_ij(1)*d_ij(2)
      J_ij(2,2) = d_ij(2)*d_ij(2)

      ! Scale local Jacobian matrix
      J_ij = 3 * J_ij / d2_ij**5

      ! Update local Jacobian matrix
      J_ij(1,1) = d2_ij**2 - J_ij(1,1)
      J_ij(1,2) =          - J_ij(1,2)
      J_ij(2,1) =          - J_ij(2,1)
      J_ij(2,2) = d2_ij**2 - J_ij(2,2)

      ii = Kdiagonal(i)
      Dmatrix(:,:,ii) = Dmatrix(:,:,ii) + J_ij


      return

      ! Loop over all elements to assemble the repulsion forces
      elems: do iel = 1, nel
        
        ! Loop over all corners
        do ive = 1, 3
          
          ! Get global vertex numbers
          i = p_IverticesAtElement(mod(ive-1, 3)+1, iel)
          j = p_IverticesAtElement(mod(ive,   3)+1, iel)
          k = p_IverticesAtElement(mod(ive+1, 3)+1, iel)

          ! Compute position of initial projection point
          d_ij  = DvertexCoords(:,k) - DvertexCoords(:,j)
          d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)
          dprj  = d_ij(1)*(DvertexCoords(1,i)-DvertexCoords(1,j)) +&
                  d_ij(2)*(DvertexCoords(2,i)-DvertexCoords(2,j))
          Dproj = DvertexCoords(:,j) + dprj/d2_ij * d_ij

          ! Compute square of rest length of the projection spring
          dlength0 = (DvertexCoords(1,i)-Dproj(1))**2 + (DvertexCoords(2,i)-Dproj(2))**2

          ! Compute position of effective projection point
          d_ij  = Dparticles(:,k) - Dparticles(:,j)
          d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)
          dprj  = d_ij(1)*(Dparticles(1,i)-Dparticles(1,j)) +&
                  d_ij(2)*(Dparticles(2,i)-Dparticles(2,j))
          Dproj = Dparticles(:,j) + dprj/d2_ij * d_ij

          ! Compute direction and effective length of the projection spring
          d_ij = Dproj-Dparticles(:,i)
          dlength = sqrt(d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2))

          ! Compute force vector for repulsion spring
          f_ij = SPRING_STIFFNESS * (dlength-dlength0/dlength) * d_ij / dlength
          
          ! Apply force vector to node i
          Dforces(:,i) = Dforces(:,i) + f_ij


          ! Update i-th row of matrix
          ii = Kdiagonal(i)
          Dmatrix(1,1,ii) = Dmatrix(1,1,ii) + SPRING_STIFFNESS * (dlength0/dlength**2 - 1.0_DP) * Dparticles(1,i)
          Dmatrix(2,2,ii) = Dmatrix(2,2,ii) + SPRING_STIFFNESS * (dlength0/dlength**2 - 1.0_DP) * Dparticles(2,i)

!!$          do ij = Kld(i), Kld(i+1)-1
!!$            
!!$            ! Update j-th column
!!$            if (Kcol(ij) .eq. j) then
!!$              Dmatrix(1,1,ij) = Dmatrix(1,1,ij) + SPRING_STIFFNESS * (dlength0/dlength**2 - 1.0_DP) * Dparticles(1,j) * dprj/d2_ij
!!$              Dmatrix(2,2,ij) = Dmatrix(2,2,ij) + SPRING_STIFFNESS * (dlength0/dlength**2 - 1.0_DP) * Dparticles(2,j) * dprj/d2_ij
!!$            end if
!!$
!!$            ! Update k-th column
!!$            if (Kcol(ij) .eq. k) then
!!$              Dmatrix(1,1,ij) = Dmatrix(1,1,ij) - SPRING_STIFFNESS * (dlength0/dlength**2 - 1.0_DP) * Dparticles(1,k) * dprj/d2_ij
!!$              Dmatrix(2,2,ij) = Dmatrix(2,2,ij) - SPRING_STIFFNESS * (dlength0/dlength**2 - 1.0_DP) * Dparticles(2,k) * dprj/d2_ij
!!$            end if
!!$          end do

        end do
      end do elems
      
    end subroutine assembleForces

    !***************************************************************************
    
    subroutine setFixedPoints(Kld, Kcol, Kdiagonal, na, neq, Dmatrix, Dparticles, Dforces)

      ! This subroutine nullifies the force vector for fixed points
      ! and replaces the corresponding row of the matrix by that of
      ! the identity matrix.

      real(DP), dimension(NDIM2D,neq), intent(IN) :: Dparticles
      integer, dimension(:), intent(IN) :: Kld, Kcol, Kdiagonal
      integer, intent(IN) :: na, neq
      
      real(DP), dimension(NDIM2D,NDIM2D,na), intent(INOUT) :: Dmatrix
      real(DP), dimension(NDIM2D,neq), intent(INOUT) :: Dforces

      ! local variables
      integer :: i,ii,ij


      ! Loop over all rows
      rows: do i = 1, neq

        if (Dparticles(2,i) .ge. 1.5_DP .or. Dparticles(2,i) .le. -1.5_DP .or.&
            Dparticles(1,i) .ge. 1.9_DP .or. Dparticles(1,i) .le.  0.1_DP) then

          ! Nullify forces
          Dforces(:,i) = 0.0_DP

          ! Nullify row in matrix
          do ij = Kld(i), Kld(i+1)-1
            Dmatrix(:,:,ij) = 0.0_DP
          end do

          ! Set identity matrix at diagonal block
          ii = Kdiagonal(i)
          Dmatrix(1,1,ii) = 1.0_DP
          Dmatrix(2,2,ii) = 1.0_DP

        end if
      end do rows
      
    end subroutine setFixedPoints

  end subroutine bloodflow_redistMeshPoints

end module bloodflow

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
  real(DP), parameter :: POINT_COLLAPSE_TOLERANCE = 1e-3_DP
  
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

    do iadapt = 1, 100

      ! Evaluate the indicator
      call bloodflow_evalIndicator(rbloodflow)
      
      ! Check if there are still unresolved object points
      if (rbloodflow%nunresolvedObjectPoints .eq. 0) return
      
      ! If so, then adapt the computational mesh
      call hadapt_performAdaptation(rbloodflow%rhadapt, rbloodflow%rindicator)
      
      ! Generate raw mesh from adaptation structure
      call hadapt_generateRawMesh(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    end do

    ! If we end up here, then the maximum admissible refinement level
    ! does not suffice to resolve all object points
    call output_line('Some objects points could not be resolved!',&
        OU_CLASS_ERROR,OU_MODE_STD,'bloodflow_adaptObject')
    call output_line('Increase the maximum admissible refinement level!',&
        OU_CLASS_ERROR,OU_MODE_STD,'bloodflow_adaptObject')
    call sys_halt()

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
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dindicator
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
    integer, dimension(:), pointer :: IelementPatch, IobjectNodes
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(TRIA_NVETRI2D) :: DbarycCoords
    real(DP), dimension(NDIM2D) :: Dpoint0Coords, Dpoint1Coords,DIntersectionCoords
    real(DP) :: distance
    integer :: ive,jve,iel,jel,i,i1,i2,i3,istatus,idx,ipoint,ipatch,npatch,ipos
    
    
    ! Release the indicator vector (if exists)
    call lsyssc_releaseVector(rbloodflow%rindicator)

    ! Create new indicator vector
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Dindicator)

    ! Allocate temporal memory for local element patch
    allocate(IelementPatch(rbloodflow%rtriangulation%NNelAtVertex))

    ! Initialize list of elements
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
      if (p_Dindicator(iel) .eq. 0) call list_appendToList(relementList, iel, ipos)
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
          if (p_Dindicator(jel) .eq. 0) call list_appendToList(relementList, jel, ipos)
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
          if (p_Dindicator(jel) .eq. 0) call list_appendToList(relementList, jel, ipos)
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
        if (p_Dindicator(jel) .eq. 0) call list_appendToList(relementList, jel, ipos)
        p_Dindicator(jel) = 1
      end do
      
    case (4:6)
      ! Previous point was located on the edge of the element, thus,
      ! mark adjacent element and create two element patch
      jel = p_IneighboursAtElement(istatus-3, iel)
      
      if (jel .ne. 0) then
        if (p_Dindicator(jel) .eq. 0) call list_appendToList(relementList, jel, ipos)
        p_Dindicator(jel) = 1
      end if
    end select

    ! Deallocate temporal memory
    deallocate(IelementPatch)

    !---------------------------------------------------------------------------
    ! (3) Find mesh points which can be projected onto the object and
    !     generate the list of elements which need to be refined.
    !---------------------------------------------------------------------------

    ! Initialize the number of unresolved object points
    rbloodflow%nunresolvedObjectPoints = 0

    ! Get position/content of first list entry
    ipos = list_getNextInList(relementList, .true.)

    ! Loop over all points of the object
    point: do ipoint = 1, size(p_DobjectCoords,2)
      
      ! Iterate over all elements in the vicinity of the object
      do while(ipos .ne. LNULL) 
        
        ! Get element number
        call list_getByPosition(relementList, ipos, iel)
        
        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        ! Compute barycentric coordinates of the object point
        ! using the current triangle as reference element 
        call getBarycentricCoords(DtriaCoords, p_DobjectCoords(:,ipoint), DbarycCoords)
        
        ! Check if point is inside the element
        if (any(DbarycCoords < -SYS_EPSREAL) .or. any(DbarycCoords > 1+SYS_EPSREAL)) then

          ! If not, then proceed to the next element
          ipos = list_getNextInList(relementList, .false.)

        else
          
          ! Compute the radius of the circumcircle for the point, that
          ! is, half the distance to its neighboring points on the object
          if (ipoint .eq. 1) then
            distance = 0.5 * sqrt( (p_DobjectCoords(1,ipoint)-p_DobjectCoords(1,ipoint+1))**2 +&
                                   (p_DobjectCoords(2,ipoint)-p_DobjectCoords(2,ipoint+1))**2 )
          elseif (ipoint .eq. size(p_DobjectCoords,2)) then
            distance = 0.5 * sqrt( (p_DobjectCoords(1,ipoint)-p_DobjectCoords(1,ipoint-1))**2 +&
                                   (p_DobjectCoords(2,ipoint)-p_DobjectCoords(2,ipoint-1))**2 )
          else
            distance = 0.5 * min( sqrt( (p_DobjectCoords(1,ipoint)-p_DobjectCoords(1,ipoint+1))**2 +&
                                        (p_DobjectCoords(2,ipoint)-p_DobjectCoords(2,ipoint+1))**2 ),&
                                  sqrt( (p_DobjectCoords(1,ipoint)-p_DobjectCoords(1,ipoint-1))**2 +&
                                        (p_DobjectCoords(2,ipoint)-p_DobjectCoords(2,ipoint-1))**2 ))
          end if
          
          ! Check if the element which contains the object point is
          ! completely contained in the circumcircle of the point
          do ive = 1, 3
            i = p_IverticesAtElement(ive, iel)
            if (sqrt( (p_DobjectCoords(1,ipoint)-p_DvertexCoords(1,i))**2 +&
                      (p_DobjectCoords(2,ipoint)-p_DvertexCoords(2,i))**2 ) > distance) then

              ! Increase number of unresolved object points
              rbloodflow%nunresolvedObjectPoints = rbloodflow%nunresolvedObjectPoints+1

              ! Mark element for h-refinement
              p_Dindicator(iel) = 2
              
              ! Proceed to the next object point
              cycle point
            end if
          end do

          ! Great, we have found a single element contained in the
          ! circumcircle of the object point. Now, we can project one
          ! of its corner vertices to the point and we are done.
          p_Dindicator(iel) = 0

          ! Determine the corner which largest barycentric coordinate
          ive = maxloc(DbarycCoords, 1); i = p_IverticesAtElement(ive, iel)

          ! Project the mesh point onto the object point
          p_DvertexCoords(:,i) = p_DobjectCoords(:,ipoint)
          
          ! Proceed to the next object point
          cycle point
        end if
      end do
    end do point
    
    ! Release list of elements
    call list_releaseList(relementList)
    
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
      area  = 0.5 * (-TriaCoords(1,2)*TriaCoords(2,1) + TriaCoords(1,3)*TriaCoords(2,1) +&
                      TriaCoords(1,1)*TriaCoords(2,2) - TriaCoords(1,3)*TriaCoords(2,2) -&
                      TriaCoords(1,1)*TriaCoords(2,3) + TriaCoords(1,2)*TriaCoords(2,3) )

      area1 = 0.5 * (-TriaCoords(1,2)*P(2) + TriaCoords(1,3)*P(2) +&
                      P(1)*TriaCoords(2,2) - TriaCoords(1,3)*TriaCoords(2,2) -&
                      P(1)*TriaCoords(2,3) + TriaCoords(1,2)*TriaCoords(2,3) )
      area2 = 0.5 * (-TriaCoords(1,3)*P(2) + TriaCoords(1,1)*P(2) +&
                      P(1)*TriaCoords(2,3) - TriaCoords(1,1)*TriaCoords(2,3) -&
                      P(1)*TriaCoords(2,1) + TriaCoords(1,3)*TriaCoords(2,1) )
      area3 = 0.5 * (-TriaCoords(1,1)*P(2) + TriaCoords(1,2)*P(2) +&
                      P(1)*TriaCoords(2,1) - TriaCoords(1,2)*TriaCoords(2,1) -&
                      P(1)*TriaCoords(2,2) + TriaCoords(1,1)*TriaCoords(2,2) )
      
      ! Compute barycentric coordinates
      BarycentricCoords(1) = area1/area                             
      BarycentricCoords(2) = area2/area
      BarycentricCoords(3) = area3/area

    end subroutine getBarycentricCoords

  end subroutine bloodflow_evalIndicator

end module bloodflow

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

    ! Handle to index array storing the first point of each object
    integer :: h_IobjectCoordsIdx = ST_NOHANDLE

    ! Handle to array storing the points of the object
    integer :: h_DobjectCoords = ST_NOHANDLE

    ! Indicator vector
    type(t_vectorScalar) :: rindicator    

    ! List of elements intersected by object
    type(t_list) :: relementList

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
      call storage_free(rbloodflow%h_IobjectCoordsIdx)
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

  subroutine bloodflow_outputStructure(rbloodflow, dtime)

!<description>

    ! This subroutine writes the content of the structure to a GMV file

!</description>

!<input>

    ! Bloodflow structure
    type(t_bloodflow), intent(in) :: rbloodflow

    ! Simulation time
    real(DP), intent(in), optional :: dtime

!</input>
!</subroutine>

    ! UCD export structure
    type(t_ucdExport) :: rexport
    real(DP), dimension(:,:), pointer :: p_DobjectCoords
    real(DP), dimension(:), pointer :: DtracerData
    integer, dimension(:), pointer :: p_IobjectCoordsIdx
    character(len=SYS_STRLEN) :: sucdfilename
    integer :: i,iobj

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! Get filename and start GMV output
    call parlst_getvalue_string(rbloodflow%rparlist,&
        'Output', 'ucdfilename', sucdfilename)
    call ucd_startGMV(rexport, UCD_FLAG_STANDARD,&
        rbloodflow%rtriangulation,&
        trim(sucdfilename)//'.'//trim(sys_si0(ifilenumber,5)))
    
    ! Attach the objects as tracers
    call storage_getbase_int(rbloodflow%h_IobjectCoordsIdx,&
        p_IobjectCoordsIdx)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords,&
        p_DobjectCoords)
    
    ! Set tracers
    call ucd_setTracers (rexport, p_DobjectCoords)

    ! Set tracer data
    allocate(DtracerData(size(p_DobjectCoords,2)))
    
    ! Loop over points of all objects
    do iobj = 1, size(p_IobjectCoordsIdx)-1
      do i = p_IobjectCoordsIdx(iobj), p_IobjectCoordsIdx(iobj+1)-1
        DtracerData(i) = real(iobj, DP)
      end do
    end do
    
    call ucd_addTracerVariable (rexport, 'object', DtracerData)
    deallocate(DtracerData)
    
    ! Attach simulation time
    if (present(dtime))&
        call ucd_setSimulationTime (rexport,dtime)

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
    integer, dimension(:), pointer :: p_IobjectCoordsIdx
    real(DP) :: a,b,c,L,t,w
    integer, dimension(2) :: Isize
    integer :: icase,ipoint,npoints,i,n

    ! Release object(s) from previous evaluation
    if (rbloodflow%h_IobjectCoordsIdx .ne. ST_NOHANDLE) then
      call storage_free(rbloodflow%h_IobjectCoordsIdx)
    end if
    
    if (rbloodflow%h_DobjectCoords .ne. ST_NOHANDLE) then
      call storage_free(rbloodflow%h_DobjectCoords)
    end if
    
    ! Get value from parameter list
    call parlst_getvalue_int(rbloodflow%rparlist, 'Object', 'icase', icase)
    select case(icase)

    case (1)

      !-------------------------------------------------------------------------
      ! Thin object
      !-------------------------------------------------------------------------
      
      call parlst_getvalue_int(rbloodflow%rparlist, 'Object', 'npoints', npoints)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'L', L)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'c', c)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'w', w)

      ! Allocate storage for the coordinates
      Isize = (/2, npoints/)
      call storage_new('bloodflow_evalObject', 'DobjectCoords', Isize,&
          ST_DOUBLE, rbloodflow%h_DobjectCoords, ST_NEWBLOCK_ZERO)
      call storage_getbase_double2d(rbloodflow%h_DobjectCoords,&
          p_DobjectCoords)

      ! Allocate storage for the index
      call storage_new('bloodflow_evalObject', 'IobjectCoordsIdx', 2,&
          ST_INT, rbloodflow%h_IobjectCoordsIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(rbloodflow%h_IobjectCoordsIdx,&
          p_IobjectCoordsIdx)
      
      ! Initialize index
      p_IobjectCoordsIdx = (/1, npoints+1/)
      
      if (abs(w) .le. SYS_EPSREAL) then
        
        ! Create a thin (quasi 1D) object
        do ipoint = 0, npoints-1
          
          ! Compute x-coordinate
          p_DobjectCoords(1,ipoint+1) = L*ipoint/real(npoints-1, DP)
          
          ! Compute y-coordinate
          p_DobjectCoords(2,ipoint+1) = c*sin(dtime)*&
              (3*L*p_DobjectCoords(1,ipoint+1)**2 -&
                   p_DobjectCoords(1,ipoint+1)**3)
        end do
        
      else

        ! Create first part of the 2D object
        do ipoint = 0, int(npoints/2.0)-1
          
          ! Compute x-coordinate
          p_DobjectCoords(1,ipoint+1) = L*ipoint/real(int(npoints/2.0)-1, DP)
         
          ! Compute y-coordinate
          p_DobjectCoords(2,ipoint+1) = c*sin(dtime)*&
              (3*L*p_DobjectCoords(1,ipoint+1)**2 -&
                   p_DobjectCoords(1,ipoint+1)**3)
        end do
          
        ! Create second part of the 2D object
        do ipoint = int(npoints/2.0), npoints-1
          
          ! Compute x-coordinate
          p_DobjectCoords(1,ipoint+1) = L*(npoints-ipoint-1)&
              /real(npoints-int(npoints/2.0)-1, DP)

          ! Compute y-coordinate
          p_DobjectCoords(2,ipoint+1) = c*sin(dtime)*&
              (3*L*p_DobjectCoords(1,ipoint+1)**2 -&
                   p_DobjectCoords(1,ipoint+1)**3) + w
        end do

      end if
      

    case (2)

      !-------------------------------------------------------------------------
      ! Rotating ellipse
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(rbloodflow%rparlist, 'Object', 'npoints', npoints)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'a', a)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'b', b)

      ! Allocate storage for the coordinates
      Isize = (/2, npoints/)
      call storage_new('bloodflow_evalObject', 'DobjectCoords', Isize,&
          ST_DOUBLE, rbloodflow%h_DobjectCoords, ST_NEWBLOCK_ZERO)
      call storage_getbase_double2d(rbloodflow%h_DobjectCoords,&
          p_DobjectCoords)
      
      ! Allocate storage for the index
      call storage_new('bloodflow_evalObject', 'IobjectCoordsIdx', 2,&
          ST_INT, rbloodflow%h_IobjectCoordsIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(rbloodflow%h_IobjectCoordsIdx,&
          p_IobjectCoordsIdx)

      ! Initialize index
      p_IobjectCoordsIdx = (/1, npoints+1/)
      
      ! Create ellipse
      do ipoint = 0, npoints-1
        
        ! Compute parameter value
        t = 2*SYS_PI*ipoint/real(npoints-1, DP)
        
        ! Compute x-coordinate
        p_DobjectCoords(1,ipoint+1) = cos(dtime)*a*cos(t)-sin(dtime)*b*sin(t)
        
        ! Compute y-coordinate
        p_DobjectCoords(2,ipoint+1) = sin(dtime)*a*cos(t)+cos(dtime)*b*sin(t)
      end do
      
    case (3)
  
      !-------------------------------------------------------------------------
      ! Rotating rotor
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(rbloodflow%rparlist, 'Object', 'n', n)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'L', L)
      call parlst_getvalue_double(rbloodflow%rparlist, 'Object', 'w', w)

      ! Allocate storage for the coordinates
      Isize = (/2, 4*n/)
      call storage_new('bloodflow_evalObject', 'DobjectCoords', Isize,&
          ST_DOUBLE, rbloodflow%h_DobjectCoords, ST_NEWBLOCK_ZERO)
      call storage_getbase_double2d(rbloodflow%h_DobjectCoords,&
          p_DobjectCoords)
      
      ! Allocate storage for the index
      call storage_new('bloodflow_evalObject', 'IobjectCoordsIdx', n+1,&
          ST_INT, rbloodflow%h_IobjectCoordsIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(rbloodflow%h_IobjectCoordsIdx,&
          p_IobjectCoordsIdx)

      ! Create each blade separately
      do i = 1, n
        
        ! Set position of first point of the object
        p_IobjectCoordsIdx(i) = (i-1)*4+1
        
        ! Compute parameter value
        t = 2*SYS_PI*(i-1)/real(n, DP)
        
        ! Compute x-coordinates
        p_DobjectCoords(1,(i-1)*4+1) = -cos(dtime+t)*w
        p_DobjectCoords(1,(i-1)*4+2) = -cos(dtime+t)*w-sin(dtime+t)*L
        p_DobjectCoords(1,(i-1)*4+3) =  cos(dtime+t)*w-sin(dtime+t)*L
        p_DobjectCoords(1,(i-1)*4+4) =  cos(dtime+t)*w
        
        ! Compute y-coordinate
        p_DobjectCoords(2,(i-1)*4+1) = -sin(dtime+t)*w
        p_DobjectCoords(2,(i-1)*4+2) = -sin(dtime+t)*w+cos(dtime+t)*L
        p_DobjectCoords(2,(i-1)*4+3) =  sin(dtime+t)*w+cos(dtime+t)*L
        p_DobjectCoords(2,(i-1)*4+4) =  sin(dtime+t)*w
        
      end do

      p_IobjectCoordsIdx(n+1) = n*4+1


    case default
      call output_line('Invalid test case!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'bloodflow_evalObject')
      call sys_halt()
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
    integer :: iadapt,npreadapt,npostadapt
    
    
    ! Initialize adaptation structure from triangulation
    call hadapt_initFromParameterlist(rbloodflow%rhadapt, rbloodflow%rparlist, 'Adaptation')
    call hadapt_initFromTriangulation(rbloodflow%rhadapt, rbloodflow%rtriangulation)

    
    ! Perform prescribed number of standard h-adaptation steps
    call parlst_getvalue_int(rbloodflow%rparlist, 'Adaptation', 'npreadapt', npreadapt)

    do iadapt = 1, npreadapt
      
      ! Evaluate the indicator
      call bloodflow_evalIndicator(rbloodflow)
      
      ! Convert it to standard h-adaptation indicator
      call bloodflow_convertRefIndicator(rbloodflow)
      
      ! Adapt the computational mesh
      call hadapt_performAdaptation(rbloodflow%rhadapt, rbloodflow%rindicator)
      
      ! Generate raw mesh from adaptation structure
      call hadapt_generateRawMesh(rbloodflow%rhadapt, rbloodflow%rtriangulation)
    end do


    ! Perform prescribed number of h-adaptation steps, required to
    ! adapt the computational mesh to the interface boundary
    call parlst_getvalue_int(rbloodflow%rparlist, 'Adaptation', 'npostadapt', npostadapt)

    do iadapt = 1, npostadapt
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
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
    integer, dimension(:), pointer :: p_Iindicator, p_IobjectCoordsIdx
    integer, dimension(:), pointer :: IelementPatch 
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Dpoint0Coords, Dpoint1Coords,Daux
    real(DP) :: dscale, dscaleMax
    integer :: ive,jve,iel,jel,i,i1,i2,i3,istatus,idx
    integer :: ipoint,ipatch,npatch,ipos,iobj
  
    
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
        rbloodflow%rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    call storage_getbase_int(&
        rbloodflow%rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    call storage_getbase_int(&
        rbloodflow%rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
    call storage_getbase_double2d(&
        rbloodflow%rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_double2d(&
        rbloodflow%h_DobjectCoords,&
        p_DobjectCoords)
    call storage_getbase_int(&
        rbloodflow%h_IobjectCoordsIdx,&
        p_IobjectCoordsIdx)
    
    
    ! Loop over all objects    
    do iobj = 1, size(p_IobjectCoordsIdx)-1
      
      ! Initialize point number
      ipoint = p_IobjectCoordsIdx(iobj)
      
      !-------------------------------------------------------------------------
      ! (1) Find element surrounding/meeting at first point of the object:
      !
      !     The algorithm is really simply. An extensive search over all
      !     element of the triangulation is performed and the first
      !     element which either surrounds the first point of the thin
      !     object or is connected to this point via a corner certex or
      !     an edge is selected.
      !-------------------------------------------------------------------------
      
      ! Get global coordinates of first point
      Dpoint0Coords =  p_DobjectCoords(:,ipoint)
      
      ! Loop over all elements in triangulation
      findElementOfFirstVertex: do&
          iel = 1, rbloodflow%rtriangulation%NEL
        
        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        ! Check if the point is 'inside' or on the boundary of the element
        call PointInTriangleTest(DtriaCoords, Dpoint0Coords,&
            POINT_EQUAL_TOLERANCE, istatus)
        if (istatus .ge. 0) exit findElementOfFirstVertex
        
      end do findElementOfFirstVertex
      
      ! Mark first element for potential refinement
      p_Iindicator(iel) = ibset(p_Iindicator(iel), istatus)
      
      ! Update point number
      ipoint = ipoint+1
      
      !-------------------------------------------------------------------------
      ! (2) Find elements surrounding/meeting at all points of the object:
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
      !-------------------------------------------------------------------------

      ! Loop over all remaining points
      findElementsOfOtherVertices: do&
          while (ipoint .le. p_IobjectCoordsIdx(iobj+1)-1)
        
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
            
            ! Check if we are at the boundary
            if (jel .ne. 0) then
              
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
              
            end if
          end do
          npatch = ipatch
          
        case (4:6)               
          ! Previous point was located on the edge of the element
          jel = p_IneighboursAtElement(istatus-3, iel)
          
          ! Create two element patch if there is an adjacent element
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
            call LinesIntersectionTest(Dpoint0Coords, Dpoint1Coords,&
                DtriaCoords(:,ive), DtriaCoords(:,mod(ive,3)+1),&
                istatus, dscale)
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
            p_Iindicator(iel) = ior(p_Iindicator(iel),&
            BITFIELD_MULTI_INTERSECTION)
        
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
        
        ! Mark elements surrounding this point
        do idx = p_IelementsAtVertexIdx(i),&
                 p_IelementsAtVertexIdx(i+1)-1
          
          ! Get global element number
          jel = p_IelementsAtVertex(idx)
          
          ! Check if we are at the boundary
          if (jel .ne. 0) then
            
            ! Append element to the list of elements adjacent to the object
            if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne.&
                BITFIELD_INLIST) then
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
            
          end if
        end do
        
      case (4:6)
        ! Previous point was located on the edge of the element, thus,
        ! mark adjacent element and create two element patch
        jel = p_IneighboursAtElement(istatus-3, iel)
        
        ! Mark adjacent element if it exists
        if (jel .ne. 0) then
          
          ! Append element to the list of elements adjacent to the object
          if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne.&
              BITFIELD_INLIST) then
            p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_INLIST)
            call list_appendToList(rbloodflow%relementList, jel, ipos)
          end if
          
          ! Mark element for potential refinement
          do jve = 1, TRIA_NVETRI2D
            if (p_IneighboursAtElement(jve, jel) .eq. iel) then
              
              ! Check if edge has been intersected previously
              if (btest(p_Iindicator(jel), jve+3))&
                  p_Iindicator(jel) = ior(p_Iindicator(jel),&
                  BITFIELD_MULTI_INTERSECTION)
              
              ! Mark edge as intersected
              p_Iindicator(jel) = ibset(p_Iindicator(jel), jve+3)
              exit
            end if
          end do
          
        end if
        
      end select
      
    end do

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
    real(DP), dimension(:,:), pointer :: p_DobjectCoords,&
        p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement,&
        p_IneighboursAtElement
    integer, dimension(:), pointer :: p_Iindicator, p_Imarker,&
        p_IvertexAge, p_InodalProperty
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Daux
    real(DP) :: dscale
    integer, dimension(2) :: Isize
    integer :: ive,jve,nve,iel,nel,ivt,jvt,nvt
    integer :: i1,i2,i3,istatus,ipoint,ipos,iresult
    

    ! Initialize marker structure (if required)
    if (rbloodflow%rhadapt%h_Imarker .ne. ST_NOHANDLE)&
        call storage_free(rbloodflow%rhadapt%h_Imarker)
    call storage_new('bloodflow_convertMoveRefIndicator', &
        'Imarker', rbloodflow%rindicator%NEQ, ST_INT,&
        rbloodflow%rhadapt%h_Imarker, ST_NEWBLOCK_ZERO)

    ! Set pointers
    call storage_getbase_int(&
        rbloodflow%rhadapt%h_Imarker, p_Imarker)
    call storage_getbase_int(&
        rbloodflow%rhadapt%h_IvertexAge, p_IvertexAge)
    call storage_getbase_int(&
        rbloodflow%rhadapt%h_InodalProperty, p_InodalProperty)
    call storage_getbase_int2d(&
        rbloodflow%rhadapt%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d(&
        rbloodflow%rhadapt%h_IneighboursAtElement, p_IneighboursAtElement)
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
    
    
    ! Set state of all vertices to "free". Note that vertices of the
    ! initial triangulation are always "locked", i.e. have no positive age.
    do ivt = 1, size(p_IvertexAge, 1)
      p_IvertexAge(ivt) = abs(p_IvertexAge(ivt))
    end do

    !---------------------------------------------------------------------------
    ! (0) Loop over all elements of the triangulation and lock those
    !     vertices which must not be repositioned, that is, vertices
    !     which have a younger edge neighbor and/or belong to the
    !     initial triangulation.
    !---------------------------------------------------------------------------
    
    locking: do iel = 1, rbloodflow%rhadapt%NEL
      
      ! Get number of vertices per elements
      nve = hadapt_getNVE(rbloodflow%rhadapt, iel)

      ! Loop over all edges of the element
      do ive = 1, nve
        
        ! Get global vertex numbers
        ivt = p_IverticesAtElement(ive, iel)
        jvt = p_IverticesAtElement(mod(ive, nve)+1, iel)
        
        ! Check if endpoints have different age and "lock" the older one
        if (abs(p_IvertexAge(ivt)) .lt.&
            abs(p_IvertexAge(jvt))) then
          p_IvertexAge(ivt) = -abs(p_IvertexAge(ivt))
        elseif (abs(p_IvertexAge(jvt)) .lt.&
                abs(p_IvertexAge(ivt))) then
          p_IvertexAge(jvt) = -abs(p_IvertexAge(jvt))
        end if
      end do
    end do locking
    

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

        ! One or more edges of the current element are intersected
        ! multiple times so that repositioning vertices does not
        ! suffice. Therefore, set the refinement indicator to
        ! enforce regular subdivision of the elements.
        nve = hadapt_getNVE(rbloodflow%rhadapt, iel)

        select case(nve)
        case (TRIA_NVETRI2D)
          p_Imarker(iel) = MARK_REF_TRIA4TRIA
          
        case (TRIA_NVEQUAD2D)
          p_Imarker(iel) = MARK_REF_QUAD4QUAD
          
          ! Increase number of vertices to be created by one to
          ! account for the new vertex in the center of the quad
          rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1

        case default
          call output_line('Unsupported type of element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bloodflow_performAdaptation')
          call sys_halt()
        end select
        
        ! Compute number of new vertices to be created at edge midpoints
        do ive = 1, nve
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
        ! Set the refinement indicator to regular refinement
        p_Imarker(iel) = MARK_REF_TRIA4TRIA

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

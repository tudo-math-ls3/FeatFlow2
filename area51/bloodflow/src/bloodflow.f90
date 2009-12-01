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
!# 7.) bloodflow_convertRefIndicator
!#     -> Converts the indicator function into refinement indicator
!#
!# 8.) bloodflow_performAdaptation
!#     -> Performs (r)h-adaptation
!# 
!# The following auxiliary routines are available:
!#
!# 1.) TestPointInTriangle2D
!#     -> Tests if a point is located inside a triangle
!#
!# 2.) TestTriangleLineIntersection2D
!#     -> Tests if a line segment intersects with an edge of the triangle
!#
!# 3.) getBarycentricCoords2D
!#     -> Calculates the barycentric coordinates of a point w.r.t. triangle
!#
!# 4.) signedArea2D
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
  real(DP), parameter :: POINT_EQUAL_TOLERANCE = 1e-8

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
  integer, parameter :: BITFIELD_POINT4 = ibset(0,4)

  ! Bitfield to identify the edges
  integer, parameter :: BITFIELD_EDGE1 = ibset(0,5)
  integer, parameter :: BITFIELD_EDGE2 = ibset(0,6)
  integer, parameter :: BITFIELD_EDGE3 = ibset(0,7)
  integer, parameter :: BITFIELD_EDGE4 = ibset(0,8)

  ! Bitfield to identify elements in list
  integer, parameter :: BITFIELD_INLIST = ibset(0,9)

  ! Bitfield to identify multi-intersected edges
  integer, parameter :: BITFIELD_MULTI_INTERSECTION = ibset(0,10)

  ! Bitfield used to check point intersection
  integer, parameter :: BITFIELD_POINT_INTERSECTION = BITFIELD_POINT1 +&
                                                      BITFIELD_POINT2 +&
                                                      BITFIELD_POINT3 +&
                                                      BITFIELD_POINT4

  ! Bitfield used to check edge intersection
  integer, parameter :: BITFIELD_EDGE_INTERSECTION = BITFIELD_EDGE1 +&
                                                     BITFIELD_EDGE2 +&
                                                     BITFIELD_EDGE3 +&
                                                     BITFIELD_EDGE4

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

    ! Allocate auxiliary memory
    allocate(DtracerData(size(p_DobjectCoords,2)))
    
    ! Loop over points of all objects
    do iobj = 1, size(p_IobjectCoordsIdx)-1
      do i = p_IobjectCoordsIdx(iobj), p_IobjectCoordsIdx(iobj+1)-1
        DtracerData(i) = real(iobj, DP)
      end do
    end do

    ! Set tracer data
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
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(TRIA_NVETRI2D) :: DlineParam
    real(DP) :: dparamMax
    integer, dimension(TRIA_NVETRI2D) :: Iedgestatus
    integer :: ipoint,ipos,iobj,ive,jve,iel,jel,ivt,istatus,idx,iel0,i1,i2,i3

    
    ! Release the indicator vector (if any)
    call lsyssc_releaseVector(rbloodflow%rindicator)

    ! Release the list of elements (if any)
    call list_releaseList(rbloodflow%relementList)
    
    ! Create new indicator vector as integer array
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true., ST_INT)
    call lsyssc_getbase_int(rbloodflow%rindicator, p_Iindicator)
    
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
            
      ! Loop over all elements in the triangulation
      do iel = 1, rbloodflow%rtriangulation%NEL
        
        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        ! Check if the point is 'inside' or on the boundary of the element
        call TestPointInTriangle2D(DtriaCoords,&
            p_DobjectCoords(:,ipoint), 0.0_DP, istatus)
        if (istatus .eq. 1) exit
      end do
      
      ! Loop over all remaining points
      do while (ipoint .le. p_IobjectCoordsIdx(iobj+1)-2)
        
        !-------------------------------------------------------------------------
        ! (2) Compute segment-edge intersection for all edges of the
        !     element patch surrounding the current element.
        !     Moreover, find the element number for which the distance
        !     between the intersection point and the current point is
        !     largest. Once that element is found, we know that the
        !     endpoint of the segment must be located in that element
        !     because otherwise there would be another intersection
        !     point for which the parameter would be even larger.
        !-------------------------------------------------------------------------

        ! Initialize the parameter value by zero, the next element
        ! number by the number of the current element and the local
        ! edge where the intersection takes place by zero. This
        ! settings ensures that the next point will be automatically
        ! considered in case the line segment [ipoint,ipoint+1] is
        ! completely located inside the current element.
        dparamMax = 0.0_DP;   iel0 = iel
                
        ! Loop over all vertices of the current element
        do ive = 1, tria_getNVE(p_IverticesAtElement, iel)
          
          ! Get global vertex number
          ivt = p_IverticesAtElement(ive, iel)
          
          ! Loop over all elements surrounding this point
          do idx = p_IelementsAtVertexIdx(ivt),&
                   p_IelementsAtVertexIdx(ivt+1)-1
            
            ! Get global element number
            jel = p_IelementsAtVertex(idx)
            
            ! If we are neither at the boundary not edge-adjacent to
            ! the current elment iel then update the patch
            if ((jel .ne. 0) .and. (jel .ne. iel) .and.&
                (jel .ne. p_IneighboursAtElement(ive,iel))) then
           
              ! Get vertices at element
              i1 = p_IverticesAtElement(1, jel)
              i2 = p_IverticesAtElement(2, jel)
              i3 = p_IverticesAtElement(3, jel)
              
              ! Get global coordinates of corner vertices
              DtriaCoords(:,1) = p_DvertexCoords(:,i1)
              DtriaCoords(:,2) = p_DvertexCoords(:,i2)
              DtriaCoords(:,3) = p_DvertexCoords(:,i3)
              
              ! Test if the endpoint is located in this element. If
              ! this is the case we indicate the new starting element
              ! by a negative sign since no additional point in
              ! triangle tests are required once this patch is done
              call TestPointInTriangle2D(DtriaCoords,&
                  p_DobjectCoords(:,ipoint+1), 0.0_DP, istatus)
              if (istatus .eq. 1) iel0 = -jel

              ! Test if line segment intersects with triangle edges
              call TestTriangleLineIntersection2D(DtriaCoords,&
                  p_DobjectCoords(:,ipoint:ipoint+1), 0.0_DP, .false.,&
                  Iedgestatus, DlineParam=DlineParam)

              ! Loop over all edges of the triangle
              do jve = 1, TRIA_NVETRI2D
                
                ! Check if the distance between the new intersection
                ! point and the starting point of the segment is
                ! larger than the maximum distance of previous points
                if (DlineParam(jve) .gt. dparamMax) then
                  dparamMax = DlineParam(jve)
                  if (iel0 .ge. 0) iel0 = jel
                end if
                
                if (Iedgestatus(jve) .eq. 1) then
                  ! Check if edge has been intersected previously,
                  ! then mark this element as multiply intersected one
                  if (btest(p_Iindicator(jel), jve+TRIA_MAXNVE2D))&
                      p_Iindicator(jel) = ior(p_Iindicator(jel),&
                                              BITFIELD_MULTI_INTERSECTION)
                  
                  ! Mark element for potential refinement
                  p_Iindicator(jel) = ibset(p_Iindicator(jel), jve+TRIA_MAXNVE2D)
                end if
              end do
              
              ! Append element to the list of elements adjacent to the object
              if (iand(p_Iindicator(jel), BITFIELD_INLIST) .ne. BITFIELD_INLIST) then
                p_Iindicator(jel) = ior(p_Iindicator(jel), BITFIELD_INLIST)          
                call list_appendToList(rbloodflow%relementList, jel, ipos)
              end if
            end if
          end do
        end do

        ! If the new starting element has been determined already then
        !  we do not have to perform additional point in triangle
        !  tests, and hence, we can directly proceed to the next line
        !  segment of the polygon
        if (iel0 .lt. 0) then
          iel = -iel0
          ipoint = ipoint+1
        else
          
          ! Either iel0 or its neighbor is the element where the
          ! endpoint of the line segment is located. We are on the
          ! safe side if we take iel0 as new starting element and
          ! check if the endpoint is located inside the element so
          ! that we can proceed to the next line segment.
          iel = iel0
          
          ! Get vertices at element
          i1 = p_IverticesAtElement(1, iel)
          i2 = p_IverticesAtElement(2, iel)
          i3 = p_IverticesAtElement(3, iel)
          
          ! Get global coordinates of corner vertices
          DtriaCoords(:,1) = p_DvertexCoords(:,i1)
          DtriaCoords(:,2) = p_DvertexCoords(:,i2)
          DtriaCoords(:,3) = p_DvertexCoords(:,i3)
          
          ! Check if the end point is located inside the element
          call TestPointInTriangle2D(DtriaCoords,&
              p_DobjectCoords(:,ipoint+1), 0.0_DP, istatus)
          
          ! If the end point is located in the element then we increase
          ! the point number and proceed with the next line segment
          if (istatus .eq. 1) ipoint = ipoint+1
        end if
      end do
    end do

!!$    >>> PROOF OF CONCEPT < CODE <<<
!!$    do iobj = 1, size(p_IobjectCoordsIdx)-1
!!$
!!$      ! Initialize point number
!!$      do ipoint = p_IobjectCoordsIdx(iobj),&
!!$                  p_IobjectCoordsIdx(iobj+1)-2
!!$
!!$        do iel = 1, rbloodflow%rtriangulation%NEL
!!$          
!!$          ! Get vertices at element
!!$          i1 = p_IverticesAtElement(1, iel)
!!$          i2 = p_IverticesAtElement(2, iel)
!!$          i3 = p_IverticesAtElement(3, iel)
!!$          
!!$          ! Get global coordinates of corner vertices
!!$          DtriaCoords(:,1) = p_DvertexCoords(:,i1)
!!$          DtriaCoords(:,2) = p_DvertexCoords(:,i2)
!!$          DtriaCoords(:,3) = p_DvertexCoords(:,i3)
!!$          
!!$          ! Test if line segment intersects with triangle edges
!!$          call TestTriangleLineIntersection2D(DtriaCoords,&
!!$              p_DobjectCoords(:,ipoint:ipoint+1), 0.0_DP, .false., Iedgestatus)
!!$          
!!$          do ive = 1, TRIA_NVETRI2D
!!$            if (Iedgestatus(ive) .ne. 0 ) then
!!$              
!!$              ! Check if edge has been intersected previously
!!$              if (btest(p_Iindicator(iel), ive+TRIA_MAXNVE2D))&
!!$                  p_Iindicator(iel) = ior(p_Iindicator(iel),&
!!$                                          BITFIELD_MULTI_INTERSECTION)
!!$              
!!$              ! Append element to the list of elements adjacent to the object
!!$              if (iand(p_Iindicator(iel), BITFIELD_INLIST)&
!!$                  .ne. BITFIELD_INLIST) then
!!$                p_Iindicator(iel) = ior(p_Iindicator(iel), BITFIELD_INLIST)
!!$                call list_appendToList(rbloodflow%relementList, iel, ipos)
!!$              end if
!!$              
!!$              ! Mark element for potential refinement
!!$              p_Iindicator(iel) = ibset(p_Iindicator(iel), ive+TRIA_MAXNVE2D)
!!$            end if
!!$          end do
!!$        end do
!!$      end do
!!$    end do
!!$    >>> END OF PROOF OF CONCEPT CODE <<<<

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
    type(t_vectorScalar) :: rvector
    real(DP), dimension(:), pointer :: p_Dindicator
    integer, dimension(:), pointer :: p_Iindicator
    integer :: ipos, iel

    ! Create new indicator vector in double precision
    call lsyssc_createVector(rvector,&
        rbloodflow%rtriangulation%NEL, .true., ST_DOUBLE)
    
    ! Convert the (integer) indicator vector to double precision
    call lsyssc_getbase_double(rvector, p_Dindicator)
    call lsyssc_getbase_int(rbloodflow%rindicator, p_Iindicator)
    
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
      
      if (iand(p_Iindicator(iel), BITFIELD_EDGE_INTERSECTION) .ne. 0)&
          p_Dindicator(iel) = 1.0
    end do list

    ! Replace the indicator vector by the double precision data and
    ! release the auxiliary vector
    call lsyssc_swapVectors(rvector, rbloodflow%rindicator)
    call lsyssc_releaseVector(rvector)

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
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords,&
        DpointCoords, DpointCoordsAux
    real(DP), dimension(TRIA_NVETRI2D) :: DedgeParam, DedgeParamAux
    real(DP) :: dscale
    integer, dimension(TRIA_NVETRI2D) :: Iedgestatus, IedgestatusAux
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
    ! initial triangulation are always "locked", i.e. have no positive
    ! age. Moreover, boundary vertices are also locked at the moment
    do ivt = 1, rbloodflow%rhadapt%NVT
      p_IvertexAge(ivt) = abs(p_IvertexAge(ivt))*&
          merge(1, -1, p_InodalProperty(ivt) .eq. 0)
    end do

    !---------------------------------------------------------------------------
    ! (0) Loop over all elements of the triangulation and lock those
    !     vertices which must not be repositioned, that is, vertices
    !     which have a younger edge neighbor and/or belong to the
    !     initial triangulation.
    !---------------------------------------------------------------------------
    
!!$    locking: do iel = 1, rbloodflow%rhadapt%NEL
!!$      
!!$      ! Get number of vertices per elements
!!$      nve = hadapt_getNVE(rbloodflow%rhadapt, iel)
!!$
!!$      ! Loop over all edges of the element
!!$      do ive = 1, nve
!!$        
!!$        ! Get global vertex numbers
!!$        ivt = p_IverticesAtElement(ive, iel)
!!$        jvt = p_IverticesAtElement(mod(ive, nve)+1, iel)
!!$        
!!$        ! Check if endpoints have different age and "lock" the older one
!!$        if (abs(p_IvertexAge(ivt)) .lt.&
!!$            abs(p_IvertexAge(jvt))) then
!!$          p_IvertexAge(ivt) = -abs(p_IvertexAge(ivt))
!!$        elseif (abs(p_IvertexAge(jvt)) .lt.&
!!$                abs(p_IvertexAge(ivt))) then
!!$          p_IvertexAge(jvt) = -abs(p_IvertexAge(jvt))
!!$        end if
!!$      end do
!!$    end do locking
    

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
    
!!$      ! Check if multi-intersections are present
!!$      if (iand(p_Iindicator(iel), BITFIELD_MULTI_INTERSECTION) .eq.&
!!$                                  BITFIELD_MULTI_INTERSECTION) then
!!$
!!$        ! Clear indicator
!!$        p_Iindicator(iel) = 0
!!$
!!$        ! One or more edges of the current element are intersected
!!$        ! multiple times so that repositioning vertices does not
!!$        ! suffice. Therefore, set the refinement indicator to
!!$        ! enforce regular subdivision of the elements.
!!$        nve = hadapt_getNVE(rbloodflow%rhadapt, iel)
!!$
!!$        select case(nve)
!!$        case (TRIA_NVETRI2D)
!!$          p_Imarker(iel) = MARK_REF_TRIA4TRIA
!!$          
!!$        case (TRIA_NVEQUAD2D)
!!$          p_Imarker(iel) = MARK_REF_QUAD4QUAD
!!$          
!!$          ! Increase number of vertices to be created by one to
!!$          ! account for the new vertex in the center of the quad
!!$          rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1
!!$          
!!$        case default
!!$          call output_line('Unsupported type of element!',&
!!$              OU_CLASS_ERROR,OU_MODE_STD,'bloodflow_performAdaptation')
!!$          call sys_halt()
!!$        end select
!!$        
!!$        ! Compute number of new vertices to be created at edge midpoints
!!$        do ive = 1, nve
!!$          if (p_IneighboursAtElement(ive, iel) .eq. 0) then
!!$            
!!$            ! Edge is adjacent to boundary
!!$            rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1
!!$            
!!$          elseif(p_Imarker(p_IneighboursAtElement(ive, iel)) .eq. 0) then
!!$
!!$            ! Edge has not been marked in previous steps
!!$            rbloodflow%rhadapt%increaseNVT = rbloodflow%rhadapt%increaseNVT+1
!!$
!!$          end if
!!$        end do
!!$        
!!$        ! That's it for multiply intersected edges
!!$        cycle list
!!$      end if
          
!!$      ! Remove "multi-intersection" flag from indicator
!!$      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_MULTI_INTERSECTION))
!!$
!!$      ! Remove "in list" flag from indicator
!!$      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_INLIST))
!!$      
!!$      ! Remove "corner vertices" flags from indicator (if any).
!!$      ! This may be used in future versions of this code.
!!$      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_POINT_INTERSECTION))
!!$      
!!$      ! Remove "interior vertex" flag from indicator (if any).
!!$      ! This may be used in future versions of this code.
!!$      p_Iindicator(iel) = iand(p_Iindicator(iel), not(BITFIELD_INNER))
      

      ! Check status of intersected element edges
      select case (iand(p_Iindicator(iel), BITFIELD_EDGE_INTERSECTION))

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
        
        
      case (BITFIELD_EDGE1 + BITFIELD_EDGE2,&
            BITFIELD_EDGE2 + BITFIELD_EDGE3,&
            BITFIELD_EDGE1 + BITFIELD_EDGE3,&
            BITFIELD_EDGE1, BITFIELD_EDGE2, BITFIELD_EDGE3)
        ! One or two edges are intersected!
        
        ! Get vertices at element
        i1 = p_IverticesAtElement(1, iel)
        i2 = p_IverticesAtElement(2, iel)
        i3 = p_IverticesAtElement(3, iel)
        
        ! Get global coordinates of corner vertices
        DtriaCoords(:,1) = p_DvertexCoords(:,i1)
        DtriaCoords(:,2) = p_DvertexCoords(:,i2)
        DtriaCoords(:,3) = p_DvertexCoords(:,i3)
        
        IedgeStatus  = 0
        DedgeParam   = 0.0_DP
        DpointCoords = 0.0_DP

        ! Loop over all segments and check intersection point
        do ipoint = 1, size(p_DobjectCoords,2)-1
          
          ! Test if line segment intersects with triangle edges
          call TestTriangleLineIntersection2D(DtriaCoords,&
              p_DobjectCoords(:,ipoint:ipoint+1), 0.0_DP,&
              .false., IedgestatusAux, DpointCoordsAux,&
              DedgeParam=DedgeParamAux)
          
          where (IedgestatusAux .eq. 1)
            Iedgestatus = IedgestatusAux
            DedgeParam  = DedgeParamAux
            DpointCoords(1,:) = DpointCoordsAux(1,:)
            DpointCoords(2,:) = DpointCoordsAux(2,:)
          end where
        end do
        
        ! Loop over all edges of the current element
        nve = tria_getNVE(p_IverticesAtElement, iel)
        do ive = 1, nve
          
          ! Skip edges which are not intersected
          if (Iedgestatus(ive) .ne. 1) cycle
          
          if ((DedgeParam(ive) .lt. POINT_COLLAPSE_TOLERANCE) .and.&
              (DedgeParam(ive) .gt. POINT_EQUAL_TOLERANCE)) then
            ! The "1/3"-rule applies
            ivt = p_IverticesAtElement(ive,iel)
            if (p_IvertexAge(ivt) .gt. 0) then
              ! Update coordinate of vertex in quadtree
              iresult = qtree_moveInQuadtree(&
                  rbloodflow%rhadapt%rVertexCoordinates2D,&
                  p_DvertexCoords(:, ivt), DpointCoords(:,ive))
              
              ! Adjust vertex coordinates and lock it
              p_DvertexCoords(:, ivt) = DpointCoords(:,ive)
              p_IvertexAge(ivt) = -abs(p_IvertexAge(ivt))
            end if
            
            ! Remove the refinement indicator for current edge
            p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+TRIA_MAXNVE2D)

          elseif((DedgeParam(ive) .gt. 1-POINT_COLLAPSE_TOLERANCE) .and.&
                 (DedgeParam(ive) .lt. 1-POINT_EQUAL_TOLERANCE)) then
            ! The "2/3"-rule applies
            ivt = p_IverticesAtElement(mod(ive, nve)+1,iel)
            if (p_IvertexAge(ivt) .gt. 0) then
              ! Update coordinate of vertex in quadtree
              iresult = qtree_moveInQuadtree(&
                  rbloodflow%rhadapt%rVertexCoordinates2D,&
                  p_DvertexCoords(:, ivt), DpointCoords(:,ive))
              ! Adjust vertex coordinates and lock it
              p_DvertexCoords(:,ivt) = DpointCoords(:,ive)
              p_IvertexAge(ivt) = -abs(p_IvertexAge(ivt))
            end if

            ! Remove the refinement indicator for current edge
            p_Iindicator(iel) = ibclr(p_Iindicator(iel), ive+TRIA_MAXNVE2D)
            
          else
            ! Mark element for regular refinement
            call markRefine(iel)
          end if
        end do
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

  contains

    subroutine markRefine(iel)
      integer, intent(in) :: iel

      ! local variables
      integer :: ive,ivt

      ! Check if the element has reached maximum refinement level
      if (maxval(abs(p_IvertexAge(p_IverticesAtElement(1:TRIA_NVETRI2D&
          ,iel)))) .lt. rbloodflow%rhadapt%nsubdividemax) then
        
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
      end if
    end subroutine markRefine

  end subroutine bloodflow_performAdaptation
  
  !*****************************************************************************

!<subroutine>
  
  pure subroutine TestTriangleLineIntersection2D(DtriaCoords,&
      DlineCoords, dtolerance, bquicktest, Istatus,&
      DpointCoords, DlineParam, DedgeParam)

!<description>
    
    ! This subroutine test if the given line intersects with the edges
    ! of the given triangle. If this is the case, than Istatus(ive)=1
    ! is returned, where ive=1,2,3 is the local edge number. If the
    ! edge is parallel to the line or if there is no common point then
    ! Istatus(ive)=0 is returned. If present the parameter value
    ! DlineParam(ive) is calculated as follows:
    !
    ! $ X = P + t*(Q-P) $
    !
    ! where X is the intersection point and P and Q denote the start
    ! and endpoint of the line, respectively. If present then the
    ! coordinates of the intersection point DpointCoords(ive)=X are
    ! returned. If present, the parameter value p for the edge
    !
    ! $ X = A + p*(B-A) $
    !
    ! is returned, where A and B denote two consecutive vertices if
    ! the triangle, i.e., first-second, second-third, third-first.

!</description>
    
!<input>
    
    ! Coordinates of the triangle
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(in) :: DtriaCoords

    ! Coordinates of the line
    real(DP), dimension(NDIM2D,TRIA_NVELINE1D), intent(in) :: DlineCoords

    ! Tolerance parameter
    real(DP), intent(in) :: dtolerance

    ! Flag: If .TRUE. then all edges are tested
    logical, intent(in) :: bquicktest
!</input>

!<output>

    ! Status of the intersection test
    integer, dimension(TRIA_NVETRI2D), intent(out) :: Istatus
    
    ! OPTIONAL: Coordinates of the intersection point
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(out), optional :: DpointCoords

    ! OPTIONAL: Parameter value of the intersection point
    real(DP), dimension(TRIA_NVETRI2D), intent(out), optional :: DlineParam
    real(DP), dimension(TRIA_NVETRI2D), intent(out), optional :: DedgeParam
       
!</output>

!</subroutine>

    ! Description of the algorith: Given the three vertices A,B and C
    ! of the triangle, each point X in the 2D-plane can be uniquely
    ! described by the barycentric coordinates (u,v):
    !
    ! $ X = u*A+v*B+(1-u-v)*C $
    !
    ! The line between two points P and Q is parametrized as follows:
    !
    ! $ X = P + t*(Q-P) $
    ! 
    ! with parameter value $0\le t\le 1$. A necessary condition for
    ! intersection of the line with the triangle is
    !
    ! $(A-C)*u + (B-C)*v + (Q-O)*t = P-C $
    !
    ! which yields two equations (for x- and y-components) for the
    ! three unknowns u,v and t. In addition to the inequalities
    ! 0\le t\le 1$, the third equality which ensures that the
    ! intersection point is located on an edge of the triangle reads
    !
    ! $ u=0$   or   $v=0$   or   $1-u-v=0$, that is, $u+v=1$
    !
    ! In other word, the three linear systems to be solved are
    !
    ! |A-C  B-C  P-Q| |u| |P-C|
    ! |             |*|v|=|   |   (1)
    ! | 1    0    0 | |t| | 0 |
    !
    ! |A-C  B-C  P-Q| |u| |P-C|
    ! |             |*|v|=|   |   (2)
    ! | 0    1    0 | |t| | 0 |
    !
    ! |A-C  B-C  P-Q| |u| |P-C|
    ! |             |*|v|=|   |   (3)
    ! | 1    1    0 | |t| | 1 |
    !
    
    ! local variables
    real(DP), dimension(NDIM2D) :: Daux
    real(DP) :: detu,detv,detuv,detaux,dpar,dpar1,dpar2,dbaryc,dist
    
    ! Initialization
    Istatus = 0
    if (present(DlineParam)) DlineParam = 0.0_DP
    if (present(DedgeParam)) DedgeParam = 0.0_DP
    if (present(DpointCoords)) DpointCoords = 0.0_DP

    !---------------------------------------------------------------------------
    ! Second edge BC
    !---------------------------------------------------------------------------

    ! Compute determinate for the second edge BC
    !  (B.x-C.x)*(P.y-Q.y)-(B.y-C.y)*(P.x-Q.x)
    detu = (DtriaCoords(1,2)-DtriaCoords(1,3))&
         * (DlineCoords(2,1)-DlineCoords(2,2))&
         - (DtriaCoords(2,2)-DtriaCoords(2,3))&
         * (DlineCoords(1,1)-DlineCoords(1,2))

    ! Intersection test for the second edge BC
    if (abs(detu) .gt. SYS_EPSREAL) then
      
      ! Compute auxiliary determinant with third column replaced by
      ! the right-hand side (B.x-C.x)*(P.y-C.y)-(B.y-C.y)*(P.x-C.x)
      detaux = (DtriaCoords(1,2)-DtriaCoords(1,3))&
             * (DlineCoords(2,1)-DtriaCoords(2,3))&
             - (DtriaCoords(2,2)-DtriaCoords(2,3))&
             * (DlineCoords(1,1)-DtriaCoords(1,3))
      
      ! Compute parameter value of the intersection point 
      ! located on the line through points P and Q
      dpar = detaux/detu

      ! Check if parameter value is bounded by 0 and 1
      if ((dpar .ge. dtolerance) .and. (dpar .le. 1-dtolerance)) then
        
        ! Compute auxiliary determinant with second column replaced by
        ! the right-hand side (P.x-C.x)*(P.y-Q.y)-(P.y-C.y)*(P.x-Q.x)
        detaux = (DlineCoords(1,1)-DtriaCoords(1,3))&
               * (DlineCoords(2,1)-DlineCoords(2,2))&
               - (DlineCoords(2,1)-DtriaCoords(2,3))&
               * (DlineCoords(1,1)-DlineCoords(1,2))

        ! Compute barycentric coordinate of the intersection point
        dbaryc = detaux/detu
        
        ! Check if intersection point is "inside" the triangle, that
        ! is, if the barycentric coordinate v satisfies 0 <= v <= 1
        if ((dbaryc .ge. dtolerance) .and. (dbaryc .le. 1-dtolerance)) then

          ! Segments PQ and BC intersect
          Istatus(2) = 1
          if (present(DlineParam)) DlineParam(2) = dpar
          if (present(DedgeParam)) DedgeParam(2) = 1-dbaryc
          if (present(DpointCoords)) DpointCoords(:,2)=&
              DlineCoords(:,1)+dpar*(DlineCoords(:,2)-DlineCoords(:,1))

          ! Return, if we only need to find any intersection point
          if (bquicktest) return
        end if
        
      end if

    else

      ! Second edge is parallel to the line segment. Test if they
      ! coincide, that is, if their distance equals zero
      dist = abs( (DtriaCoords(1,2)-DlineCoords(1,1))*&
                  (DlineCoords(2,1)-DlineCoords(2,2))+&
                  (DtriaCoords(2,2)-DlineCoords(2,1))*&
                  (DlineCoords(1,2)-DlineCoords(1,1)) )/&
            sqrt( (DlineCoords(2,1)-DlineCoords(2,2))**2+&
                  (DlineCoords(1,2)-DlineCoords(1,1))**2 )
      
      if (dist .le. SYS_EPSREAL) then
        ! Compute parameter value for third and second vertex of the
        ! triangle expressed in terms of the line parameter form
        if (abs(DlineCoords(1,2)-DlineCoords(1,1)) .gt.&
            abs(DlineCoords(2,2)-DlineCoords(2,1))) then
          dpar2 = (DtriaCoords(1,3)-DlineCoords(1,1))/&
                  (DlineCoords(1,2)-DlineCoords(1,1))
          dpar1 = (DtriaCoords(1,2)-DlineCoords(1,1))/&
                  (DlineCoords(1,2)-DlineCoords(1,1))
        else
          dpar2 = (DtriaCoords(2,3)-DlineCoords(2,1))/&
                  (DlineCoords(2,2)-DlineCoords(2,1))
          dpar1 = (DtriaCoords(2,2)-DlineCoords(2,1))/&
                  (DlineCoords(2,2)-DlineCoords(2,1))
        end if

        dpar = max(dpar1,dpar2)

        if ((min(dpar1,dpar2) .gt. 1-dtolerance) .or.&
            (dpar .lt. dtolerance)) then
          ! Segment PQ and BC have no common subinterval, and hence,
          ! no parameter values need to be computed
        else
          Istatus(2) = 2
          dpar = min(dpar,1.0_DP)
          Daux = DlineCoords(:,1)+dpar*(DlineCoords(:,2)-DlineCoords(:,1))
          if (present(DlineParam)) DlineParam(2) = dpar
          if (present(DpointCoords)) DpointCoords(:,2) = Daux
          if (present(DedgeParam)) then
            if (abs(DtriaCoords(1,3)-DtriaCoords(1,2)) .gt.&
                abs(DtriaCoords(2,3)-DtriaCoords(2,2))) then
              DedgeParam(2) = (Daux(1)-DtriaCoords(1,2))/&
                     (DtriaCoords(1,3)-DtriaCoords(1,2))
            else
              DedgeParam(2) = (Daux(2)-DtriaCoords(2,2))/&
                     (DtriaCoords(2,3)-DtriaCoords(2,2))
            end if
          end if
        end if
      end if
    end if
    
    !---------------------------------------------------------------------------
    ! Third edge CA
    !---------------------------------------------------------------------------
    
    ! Compute determinate for the third edge CA
    !  -(A.x-C.x)*(P.y-Q.y)+(A.y-C.y)*(P.x-Q.x)
    detv = -(DtriaCoords(1,1)-DtriaCoords(1,3))&
         *  (DlineCoords(2,1)-DlineCoords(2,2))&
         +  (DtriaCoords(2,1)-DtriaCoords(2,3))&
         *  (DlineCoords(1,1)-DlineCoords(1,2))

    ! Intersection test for the third edge CA
    if (abs(detv) .gt. SYS_EPSREAL) then
      
      ! Compute auxiliary determinant twith third column replaced by
      ! the right-hand side -(A.x-C.x)*(P.y-C.y)+(A.y-C.y)*(P.x-C.x)
      detaux = -(DtriaCoords(1,1)-DtriaCoords(1,3))&
             *  (DlineCoords(2,1)-DtriaCoords(2,3))&
             +  (DtriaCoords(2,1)-DtriaCoords(2,3))&
             *  (DlineCoords(1,1)-DtriaCoords(1,3))

      ! Compute parameter value of intersection point
      ! located on the line through points P and Q
      dpar = detaux/detv

      ! Check if parameter value is bounded by 0 and 1
      if ((dpar .ge. dtolerance) .and. (dpar .le. 1-dtolerance)) then
        
        ! Compute auxiliary determinant with first column replaced by
        ! the right-hand side -(P.x-C.x)*(P.y-Q.y)+(P.y-C.y)*(P.x-Q.x)
        detaux = -(DlineCoords(1,1)-DtriaCoords(1,3))&
               *  (DlineCoords(2,1)-DlineCoords(2,2))&
               +  (DlineCoords(2,1)-DtriaCoords(2,3))&
               *  (DlineCoords(1,1)-DlineCoords(1,2))
      
        ! Compute barycentric coordinate of the intersection point
        dbaryc = detaux/detv
        
        ! Check if intersection point is "inside" the triangle, that
        ! is, if the barycentric coordinate u satisfies 0 <= u <= 1
        if ((dbaryc .ge. dtolerance) .and. (dbaryc .le. 1-dtolerance)) then

          ! Segments PQ and AC intersect
          Istatus(3) = 1
          if (present(DlineParam)) DlineParam(3) = dpar
          if (present(DedgeParam)) DedgeParam(3) = dbaryc
          if (present(DpointCoords)) DpointCoords(:,3)=&
              DlineCoords(:,1)+dpar*(DlineCoords(:,2)-DlineCoords(:,1))

          ! Return, if we only need to find any intersection point
          if (bquicktest) return
        end if
        
      end if

    else

      ! Third edge is parallel to the line segment. Test if they
      ! coincide, that is, if their distance equals zero
       dist = abs( (DtriaCoords(1,3)-DlineCoords(1,1))*&
                  (DlineCoords(2,1)-DlineCoords(2,2))+&
                  (DtriaCoords(2,3)-DlineCoords(2,1))*&
                  (DlineCoords(1,2)-DlineCoords(1,1)) )/&
            sqrt( (DlineCoords(2,1)-DlineCoords(2,2))**2+&
                  (DlineCoords(1,2)-DlineCoords(1,1))**2 )
      
      if (dist .le. SYS_EPSREAL) then
        ! Compute parameter values for first and third vertex of the
        ! triangle expressed in terms of the line parameter form
        if (abs(DlineCoords(1,2)-DlineCoords(1,1)) .gt.&
            abs(DlineCoords(2,2)-DlineCoords(2,1))) then
          dpar2 = (DtriaCoords(1,1)-DlineCoords(1,1))/&
                  (DlineCoords(1,2)-DlineCoords(1,1))
          dpar1 = (DtriaCoords(1,3)-DlineCoords(1,1))/&
                  (DlineCoords(1,2)-DlineCoords(1,1))
        else
          dpar2 = (DtriaCoords(2,1)-DlineCoords(2,1))/&
                  (DlineCoords(2,2)-DlineCoords(2,1))
          dpar1 = (DtriaCoords(2,3)-DlineCoords(2,1))/&
                  (DlineCoords(2,2)-DlineCoords(2,1))
        end if

        dpar = max(dpar1,dpar2)

        if ((min(dpar1,dpar2) .gt. 1-dtolerance) .or.&
            (dpar .lt. dtolerance)) then
          ! Segment PQ and CA have no common subinterval, and hence,
          ! no parameter values need to be computed
        else
          Istatus(3) = 2
          dpar = min(dpar,1.0_DP)
          Daux = DlineCoords(:,1)+dpar*(DlineCoords(:,2)-DlineCoords(:,1))
          if (present(DlineParam)) DlineParam(3) = dpar
          if (present(DpointCoords)) DpointCoords(:,3) = Daux
          if (present(DedgeParam)) then
            if (abs(DtriaCoords(1,1)-DtriaCoords(1,3)) .gt.&
                abs(DtriaCoords(2,1)-DtriaCoords(2,3))) then
              DedgeParam(3) = (Daux(1)-DtriaCoords(1,3))/&
                     (DtriaCoords(1,1)-DtriaCoords(1,3))
            else
              DedgeParam(3) = (Daux(2)-DtriaCoords(2,3))/&
                     (DtriaCoords(2,1)-DtriaCoords(2,3))
            end if
          end if
        end if
      end if
    end if
    
    !---------------------------------------------------------------------------
    ! First edge AB
    !---------------------------------------------------------------------------

    ! The treatment of the first edge AB requires most computational
    ! work, and hence, it is checked after the computationally less
    ! expensive tests for the second and third edge.
    
    ! Compute determinate for the first edge AB 
    detuv = detu + detv
    
    ! Intersection test for the first edge AB. 
    if (abs(detuv) .gt. SYS_EPSREAL) then
      
      ! Compute auxiliary determinant with third column
      ! replaced by the right-hand side
      !  (B.x-C.x)*(P.y-C.y)-(B.y-C.y)*(P.x-C.x)
      ! -(A.x-C.x)*(P.y-C.y)+(A.y-C.y)*(P.x-C.x)
      ! +(A.x-C.x)*(B.y-C.y)-(A.y-C.y)*(B.x-C.x)
      detaux = (DtriaCoords(1,2)-DtriaCoords(1,3))&
             * (DlineCoords(2,1)-DtriaCoords(2,3))&
             - (DtriaCoords(2,2)-DtriaCoords(2,3))&
             * (DlineCoords(1,1)-DtriaCoords(1,3))&
             - (DtriaCoords(1,1)-DtriaCoords(1,3))&
             * (DlineCoords(2,1)-DtriaCoords(2,3))&
             + (DtriaCoords(2,1)-DtriaCoords(2,3))&
             * (DlineCoords(1,1)-DtriaCoords(1,3))&
             + (DtriaCoords(1,1)-DtriaCoords(1,3))&
             * (DtriaCoords(2,2)-DtriaCoords(2,3))&
             - (DtriaCoords(2,1)-DtriaCoords(2,3))&
             * (DtriaCoords(1,2)-DtriaCoords(1,3))

      ! Compute parameter value of intersection point
      ! located on the line through points P and Q
      dpar = detaux/detuv

      ! Check if parameter value is bounded by 0 and 1
      if ((dpar .ge. dtolerance) .and. (dpar .le. 1-dtolerance)) then

        ! Compute auxiliary determinant with first column replaced by
        ! the right-hand side
        !  (P.x-C.x)*(P.y-Q.y)-(P.y-C.y)*(P.x-Q.x)
        ! -(A.x-C.x)*(P.y-Q.y)+(A.y-C.y)*(P.x-Q.x)
        detaux = (DlineCoords(1,1)-DtriaCoords(1,3))&
               * (DlineCoords(2,1)-DlineCoords(2,2))&
               - (DlineCoords(2,1)-DtriaCoords(2,3))&
               * (DlineCoords(1,1)-DlineCoords(1,2))&
               - (DtriaCoords(1,1)-DtriaCoords(1,3))&
               * (DlineCoords(2,1)-DlineCoords(2,2))&
               + (DtriaCoords(2,1)-DtriaCoords(2,3))&
               * (DlineCoords(1,1)-DlineCoords(1,2))
        
        ! Compute barycentric coordinate of the intersection point
        dbaryc = detaux/detuv
        
        ! Check if intersection point is "inside" the triangle, that
        ! is, if the barycentric coordinate u satisfies 0 <= u <= 1
        if ((dbaryc .ge. dtolerance) .and. (dbaryc .le. 1-dtolerance)) then
          
          ! Segments PQ and AB intersect
          Istatus(1) = 1
          if (present(DlineParam)) DlineParam(1) = dpar
          if (present(DedgeParam)) DedgeParam(1) = dbaryc
          if (present(DpointCoords)) DpointCoords(:,1)=&
              DlineCoords(:,1)+dpar*(DlineCoords(:,2)-DlineCoords(:,1))
        end if
        
      end if

    else
      
      ! First edge is parallel to the line segment. Test if they
      ! coincide, that is, if their distance equals zero
      dist = abs( (DtriaCoords(1,1)-DlineCoords(1,1))*&
                  (DlineCoords(2,1)-DlineCoords(2,2))+&
                  (DtriaCoords(2,1)-DlineCoords(2,1))*&
                  (DlineCoords(1,2)-DlineCoords(1,1)) )/&
            sqrt( (DlineCoords(2,1)-DlineCoords(2,2))**2+&
                  (DlineCoords(1,2)-DlineCoords(1,1))**2 )

      if (dist .le. SYS_EPSREAL) then
        ! Compute parameter value for first and second vertex of the
        ! triangle expressed in terms of the line parameter form
        if (abs(DlineCoords(1,2)-DlineCoords(1,1)) .gt.&
            abs(DlineCoords(2,2)-DlineCoords(2,1))) then
          dpar2 = (DtriaCoords(1,2)-DlineCoords(1,1))/&
                  (DlineCoords(1,2)-DlineCoords(1,1))
          dpar1 = (DtriaCoords(1,1)-DlineCoords(1,1))/&
                  (DlineCoords(1,2)-DlineCoords(1,1))
        else
          dpar2 = (DtriaCoords(2,2)-DlineCoords(2,1))/&
                  (DlineCoords(2,2)-DlineCoords(2,1))
          dpar1 = (DtriaCoords(2,1)-DlineCoords(2,1))/&
                  (DlineCoords(2,2)-DlineCoords(2,1))
        end if

        dpar = max(dpar1,dpar2)

        if ((min(dpar1,dpar2) .gt. 1-dtolerance) .or.&
            (dpar .lt. dtolerance)) then
          ! Segment PQ and AB have no common subinterval, and hence,
          ! no parameter values need to be computed
        else
          Istatus(1) = 2
          dpar = min(dpar,1.0_DP)
          Daux = DlineCoords(:,1)+dpar*(DlineCoords(:,2)-DlineCoords(:,1))
          if (present(DlineParam)) DlineParam(1) = dpar
          if (present(DpointCoords)) DpointCoords(:,1) = Daux
          if (present(DedgeParam)) then
            if (abs(DtriaCoords(1,2)-DtriaCoords(1,1)) .gt.&
                abs(DtriaCoords(2,2)-DtriaCoords(2,1))) then
              DedgeParam(1) = (Daux(1)-DtriaCoords(1,1))/&
                     (DtriaCoords(1,2)-DtriaCoords(1,1))
            else
              DedgeParam(1) = (Daux(2)-DtriaCoords(2,1))/&
                     (DtriaCoords(2,2)-DtriaCoords(2,1))
            end if
          end if
        end if
      end if
    end if
    
  end subroutine TestTriangleLineIntersection2D

  !*****************************************************************************

!<subroutine>
  
  pure subroutine TestPointInTriangle2D(DtriaCoords, DpointCoords,&
      dtolerance, istatus)

!<description>
    
    ! This subroutine calculates the relation of the given point P
    ! compare to the triangle which is defined by its three corners
    ! The meaning of the resulting istatus is as follows:
    !
    ! istatus:
    !  = 0 : point P is outside of the triangle
    !  = 1 : point P is located inside the triangle

!</description>
    
!<input>
    
    ! Coordinates of the triangle
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(in) :: DtriaCoords

    ! Coordinates of the point
    real(DP), dimension(NDIM2D), intent(in) :: DpointCoords

    ! Tolerance for point collapse
    real(DP), intent(in) :: dtolerance

!</input>

!<output>

    ! Status of the test
    integer, intent(out) :: istatus

!</output>

!</subroutine>
    
    ! local variables
    real(DP) :: area,area1,area2,area3
    real(DP) :: baryc1,baryc2,baryc3
    
    ! Compute area of global and sub-triangles
    area  = signedArea2D(DtriaCoords(:,1), DtriaCoords(:,2), DtriaCoords(:,3))

    area1 = signedArea2D(DpointCoords, DtriaCoords(:,2), DtriaCoords(:,3))
    area2 = signedArea2D(DpointCoords, DtriaCoords(:,3), DtriaCoords(:,1))
    area3 = signedArea2D(DpointCoords, DtriaCoords(:,1), DtriaCoords(:,2))
    
    ! Compute barycentric coordinates
    baryc1 = area1/area
    baryc2 = area2/area
    baryc3 = area3/area

    ! Determine status
    if ((min(baryc1,baryc2,baryc3) .ge. dtolerance) .and.&
        (max(baryc1,baryc2,baryc3) .le. 1-dtolerance)) then
      istatus = 1
    else
      istatus = 0
    end if
    
  end subroutine TestPointInTriangle2D
  
  !*****************************************************************************

!<subroutine>
  
  pure subroutine getBarycentricCoords2D(TriaCoords, P, BarycentricCoords)

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
    area  = signedArea2D(TriaCoords(:,1), TriaCoords(:,2), TriaCoords(:,3))
    area1 = signedArea2D(P, TriaCoords(:,2), TriaCoords(:,3))
    area2 = signedArea2D(P, TriaCoords(:,3), TriaCoords(:,1))
    area3 = signedArea2D(P, TriaCoords(:,1), TriaCoords(:,2))
    
    ! Compute barycentric coordinates
    BarycentricCoords(1) = area1/area                             
    BarycentricCoords(2) = area2/area
    BarycentricCoords(3) = area3/area
    
  end subroutine getBarycentricCoords2D
  
  !*****************************************************************************

!<function>
  
  pure function signedArea2D(P1,P2,P3) result(area)
    
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
    
  end function signedArea2D

end module bloodflow

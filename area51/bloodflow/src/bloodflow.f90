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
    real(DP), dimension(:), pointer :: p_Ddata
    character(len=SYS_STRLEN) :: sucdfilename

    ! Get filename and start GMV output
    call parlst_getvalue_string(rbloodflow%rparlist, 'Output', 'ucdfilename', sucdfilename)
    call ucd_startGMV(rexport, UCD_FLAG_STANDARD, rbloodflow%rtriangulation, sucdfilename)

    ! Attach the indicator vector
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Ddata)
    call ucd_addVariableElementBased(rexport, 'Indicator', UCD_VAR_STANDARD, p_Ddata)
    
    ! Write UCD exporter to file
    call ucd_write(rexport)
    call ucd_release(rexport)

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
    call storage_new('bloodflow_evalObject', 'DobjectCoords', (/2, npoints/),&
                     ST_DOUBLE, rbloodflow%h_DobjectCoords, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)
    
    do ipoint = 0, npoints-1
      
      ! Compute x-coordinate
      p_DobjectCoords(1,ipoint+1) = L*ipoint/real(npoints-1, DP)

      ! Compute y-coordinate
      p_DobjectCoords(2,ipoint+1) = c*SIN(dtime)*(3*L*p_DobjectCoords(1,ipoint+1)**2 -&
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
    real(DP), dimension(NDIM2D,TRIA_NVETRI2D) :: DtriaCoords
    real(DP), dimension(NDIM2D) :: Dpoint
    integer, dimension(:,:), pointer :: P_IverticesAtElement
    integer :: iel,ipoint,i1,i2,i3

    
    ! Release the indicator vector
    call lsyssc_releaseVector(rbloodflow%rindicator)

    ! Create new indicator vector
    call lsyssc_createVector(rbloodflow%rindicator,&
        rbloodflow%rtriangulation%NEL, .true.)
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Dindicator)

    
    ! Set pointers
    call storage_getbase_int2d(rbloodflow%rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
    call storage_getbase_double2d(rbloodflow%rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)

    ! Loop over all elements in triangulation
    do iel = 1, rbloodflow%rtriangulation%NEL

      ! Get vertices at element
      i1 = p_IverticesAtElement(1, iel)
      i2 = p_IverticesAtElement(2, iel)
      i3 = p_IverticesAtElement(3, iel)

      ! Get global coordinates
      DtriaCoords(:,1) = p_DvertexCoords(:,i1)
      DtriaCoords(:,2) = p_DvertexCoords(:,i2)
      DtriaCoords(:,3) = p_DvertexCoords(:,i3)

      ! Loop over all points of the thin object
      do ipoint = 1, size(p_DobjectCoords,2)

        ! Get point coordinates
        Dpoint = p_DobjectCoords(:,ipoint)
        
        ! Check if point is 'inside' the element
        if (isInTriangle(DtriaCoords, Dpoint)) then
          p_Dindicator(iel) = 1
          cycle
        end if
      end do
    end do
    
  contains
    
    !***************************************************************************
    
    pure function isInTriangle(TriaCoords, P)
      
      real(DP), dimension(NDIM2D,TRIA_NVETRI2D), intent(IN) :: TriaCoords
      real(DP), dimension(NDIM2D), intent(IN) :: P
      logical :: isInTriangle

      ! local variables
      real(DP), dimension(NDIM2D) :: v0,v1,v2
      real(DP) :: dot00, dot01, dot02, dot11, dot12, invDenom, u, v

      ! Compute vectors
      v0 = TriaCoords(:,3) - TriaCoords(:,1)
      v1 = TriaCoords(:,2) - TriaCoords(:,1)
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
      
      ! Check if point is in triangle
      isInTriangle = (u > 0) .and. (v > 0) .and. (u + v < 1)

    end function isInTriangle

  end subroutine bloodflow_evalIndicator

end module bloodflow

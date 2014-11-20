!##############################################################################
!# Tutorial 007c: Identify points in the mesh, write a VTK file.
!##############################################################################

module tutorial007c

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  use meshregion
  use fparser

  use element
  use spatialdiscretisation
  use linearsystemscalar
  use ucd

  implicit none
  private
  
  public :: start_tutorial007c

contains

  ! ***************************************************************************

  subroutine start_tutorial007c

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_ucdExport) :: rexport
    type(t_vectorScalar) :: rx

    type(t_meshRegion) :: rmeshRegion
    integer :: i
    integer, dimension(:), pointer :: p_Idata
    real(DP), dimension(:), pointer :: p_Ddata
    character(len=SYS_STRLEN) :: scondition
    character(LEN=SYS_STRLEN) :: spostdir
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 007c")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 65x65-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 64, 64)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Initialise the parser.
    ! =================================

    call fparser_init()

    ! =================================
    ! Identify all points on the boundary.
    ! =================================
    
    ! The condition fetches all points in a circle around (0.5,0.5)
    scondition = "(X-0.5)^2 + (Y-0.5)^2 < 0.25^2"
    
    ! Get inner points in the circle
    call mshreg_createFromExpression(rmeshRegion, rtriangulation, &
        MSHREG_IDX_VERTEX, .false., scondition)

    ! =================================
    ! Discretise with Q1.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create a scalar vector
    ! =================================

    call lsyssc_createVector (rdiscretisation,rx)
    
    ! Clear the vector
    call lsyssc_clearVector (rx)
    
    ! Get the vector data
    call lsyssc_getbase_double (rx,p_Ddata)
    
    ! Get the indices of the points in the mesh region
    call storage_getbase_int (rmeshRegion%h_IvertexIdx,p_Idata)
    
    ! Set all points in the circle to value =1.0.
    do i=1,rmeshRegion%NVT
      p_Ddata ( p_Idata(i) ) = 1.0_DP
    end do

    ! =================================
    ! Write a VTK file with the mesh
    ! and this solution vector.
    ! =================================
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call output_line ("Writing file '"//trim(spostdir)//"/tutorial007c.vtk'.")
      
      ! Open
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial007c.vtk")
    else
      call output_line ("Writing file './post/tutorial007c.vtk'.")
      
      ! Open
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "./post/tutorial007c.vtk")
    end if
                       
    ! Pass rx as solution.
    call ucd_addVectorByVertex (rexport, "circle", UCD_VAR_STANDARD, rx)
          
    ! Write / close             
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsyssc_releaseVector (rx)
    
    ! Release the mesh region
    call mshreg_done(rmeshRegion)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Clean up the parser.
    call fparser_done()
    
    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module

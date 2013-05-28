!##############################################################################
!# Tutorial 006m: Read a greymap picture, translate into Q1 solution
!##############################################################################

module tutorial006m

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  
  use vectorio
  use pprocsolution
  use ucd

  implicit none
  private
  
  public :: start_tutorial006m

contains

  ! ***************************************************************************

  subroutine start_tutorial006m

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_vectorScalar) :: rx
    type(t_pgm) :: rpgm
    type(t_ucdExport) :: rexport
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dx
    
    ! Number of points in x- and y-direction
    integer, parameter :: npointsX = 64
    integer, parameter :: npointsY = 64
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006m")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 64x64-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, &
        0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, npointsX-1, npointsY-1)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rdiscretisation
    ! which describes the discretisation.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create a scalar vector
    ! =================================

    call lsyssc_createVector (rdiscretisation,rx)
    
    ! =================================
    ! Read the picture from the file.
    ! =================================
    
    call ppsol_readPGM(0, "./pre/cfdlogo.pgm", rpgm)

    ! =================================
    ! Translate into 2D array. As points,
    ! specify the coordinates of the
    ! points in our mesh.
    ! =================================
    
    ! Get the point coordinates
    call storage_getbase_double2D (rtriangulation%h_DvertexCoords,p_DvertexCoords)

    ! Get the DOF array
    call lsyssc_getbase_double (rx,p_Dx)
    
    ! Translate the picture into our Q_1 DOF array.   
    call ppsol_initArrayPGMDP(rpgm, p_DvertexCoords, p_Dx)

    ! =================================
    ! Release the greymap image
    ! =================================

    call ppsol_releasePGM (rpgm)
    
    ! =================================
    ! Output to files. 
    ! =================================
    
    call output_line ("Writing file 'post/tutorial006m.vtk'")

    ! Open / write / close; write the solution to a VTK file.
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,"post/tutorial006m.vtk")
    call ucd_addVectorByVertex (rexport, "solution", &
        UCD_VAR_STANDARD, rx)
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsyssc_releaseVector (rx)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module

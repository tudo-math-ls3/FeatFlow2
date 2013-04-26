!##############################################################################
!# Tutorial 013a: The analytical boundary
!##############################################################################

module tutorial013a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary

  implicit none
  private
  
  public :: start_tutorial013a

contains

  ! ***************************************************************************

  subroutine start_tutorial013a

    ! Declare some variables.
    type(t_boundary) :: rboundary
    
    integer :: nbct,nseg
    real(DP) :: dmaxpar01, dmaxparLen
    real(DP), dimension(2,4) :: Dcoords
    real(DP), dimension(2,4) :: Dnormals

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 013a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Read the underlying domain
    ! =================================

    call boundary_read_prm(rboundary, "pre/QUAD.prm")

    ! =================================
    ! Gather statistics
    ! =================================

    ! Number of boundary components
    nbct = boundary_igetNBoundComp(rboundary)
    
    ! Number of boundary segments on 1st boundary component
    nseg = boundary_igetNsegments(rboundary,1)
    
    ! Maximum parameter value of 1st boundary component,
    ! 0-1 and length parametrisation
    dmaxpar01 = boundary_dgetMaxParVal(rboundary, 1)
    dmaxparLen = boundary_dgetMaxParVal(rboundary, 1, BDR_PAR_LENGTH)
    
    ! The four corners, at parmeter values 0.0, 1.0, 2.0 and 3.0
    call boundary_getCoords(rboundary, 1, 0.0_DP, Dcoords(1,1), Dcoords(2,1))
    call boundary_getCoords(rboundary, 1, 1.0_DP, Dcoords(1,2), Dcoords(2,2))
    call boundary_getCoords(rboundary, 1, 2.0_DP, Dcoords(1,3), Dcoords(2,3))
    call boundary_getCoords(rboundary, 1, 3.0_DP, Dcoords(1,4), Dcoords(2,4))
    
    ! Normal vectors at the midpoints of the edges, parameter values 0.5, 1.5, 2.5 and 3.5
    call boundary_getNormalVec2D(rboundary, 1, 0.5_DP, Dnormals(1,1), Dnormals(2,1))
    call boundary_getNormalVec2D(rboundary, 1, 1.5_DP, Dnormals(1,2), Dnormals(2,2))
    call boundary_getNormalVec2D(rboundary, 1, 2.5_DP, Dnormals(1,3), Dnormals(2,3))
    call boundary_getNormalVec2D(rboundary, 1, 3.5_DP, Dnormals(1,4), Dnormals(2,4))

    ! =================================
    ! Print statistics
    ! =================================

    call output_line ("Domain 'QUAD.prm'. Statistics.")
    call output_line ("------------------------------")
    
    call output_line ("Number of boundary components     : "//&
        trim(sys_siL(nbct,10)))
        
    call output_lbrk()

    call output_line ("Number of boundary segments, BC1  : "//&
        trim(sys_siL(nseg,10)))

    call output_line ("Max. parameter value (0-1)        : "//&
        trim(sys_sdL(dmaxpar01,3)))

    call output_line ("Max. parameter value (length)     : "//&
        trim(sys_sdL(dmaxparLen,3)))
        
    call output_lbrk()

    call output_line ("Coordinates corner 1              : "//&
        trim(sys_sdL(Dcoords(1,1),3)) // ", " // trim(sys_sdL(Dcoords(2,1),3)))
    
    call output_line ("Coordinates corner 2              : "//&
        trim(sys_sdL(Dcoords(1,2),3)) // ", " // trim(sys_sdL(Dcoords(2,2),3)))
    
    call output_line ("Coordinates corner 3              : "//&
        trim(sys_sdL(Dcoords(1,3),3)) // ", " // trim(sys_sdL(Dcoords(2,3),3)))
    
    call output_line ("Coordinates corner 4              : "//&
        trim(sys_sdL(Dcoords(1,4),3)) // ", " // trim(sys_sdL(Dcoords(2,4),3)))
        
    call output_lbrk()

    call output_line ("Normal vector edge 1              : "//&
        trim(sys_sdL(Dnormals(1,1),3)) // ", " // trim(sys_sdL(Dnormals(2,1),3)))

    call output_line ("Normal vector edge 2              : "//&
        trim(sys_sdL(Dnormals(1,2),3)) // ", " // trim(sys_sdL(Dnormals(2,2),3)))

    call output_line ("Normal vector edge 3              : "//&
        trim(sys_sdL(Dnormals(1,3),3)) // ", " // trim(sys_sdL(Dnormals(2,3),3)))

    call output_line ("Normal vector edge 4              : "//&
        trim(sys_sdL(Dnormals(1,4),3)) // ", " // trim(sys_sdL(Dnormals(2,4),3)))

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the boundary defintiion
    call boundary_release(rboundary)
    
  end subroutine

end module

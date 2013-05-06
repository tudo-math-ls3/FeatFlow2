!##############################################################################
!# Tutorial 013b: The analytical boundary, extended version
!##############################################################################

module tutorial013b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary

  implicit none
  private
  
  public :: start_tutorial013b

contains

  ! ***************************************************************************

  subroutine start_tutorial013b

    ! Declare some variables.
    type(t_boundary) :: rboundary
    
    integer :: nbct,nsega,nsegb
    real(DP) :: dpar, dparlen
    real(DP) :: dmaxpar01a, dmaxparLena, dmaxpar01b, dmaxparLenb
    real(DP), dimension(2,4) :: Dcoords
    real(DP), dimension(2,4) :: Dnormals

    real(DP), dimension(2,8) :: DcoordsCircle
    real(DP), dimension(2,8) :: DnormalsCircle
    
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 013b")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Read the underlying domain
    ! =================================

    call boundary_read_prm(rboundary, "pre/bench1.prm")

    ! =================================
    ! Gather statistics
    ! =================================

    ! Number of boundary components
    nbct = boundary_igetNBoundComp(rboundary)
    
    ! Number of boundary segments on 1st boundary component
    nsega = boundary_igetNsegments(rboundary,1)
    
    ! Maximum parameter value of 1st boundary component,
    ! 0-1 and length parametrisation
    dmaxpar01a = boundary_dgetMaxParVal(rboundary, 1)
    dmaxparLena = boundary_dgetMaxParVal(rboundary, 1, BDR_PAR_LENGTH)

    ! Number of boundary segments on 2nd boundary component (circle)
    nsegb = boundary_igetNsegments(rboundary,2)
    
    ! Maximum parameter value of 2nd boundary component,
    ! 0-1 and length parametrisation
    dmaxpar01b = boundary_dgetMaxParVal(rboundary, 2)
    dmaxparLenb = boundary_dgetMaxParVal(rboundary, 2, BDR_PAR_LENGTH)
    
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

    ! Points and corresponding normal vectors on the circle
    do i=1,8
      ! Parameter value of a point
      ! (0-1 parmetrisation, circle is one segment in the range 0-1)
      dpar = real(i-1,DP)/8.0_DP
      
      call boundary_getCoords(rboundary, 2, dpar, &
          DcoordsCircle(1,i), DcoordsCircle(2,i))
    
      call boundary_getNormalVec2D(rboundary, 2, dpar, &
          DnormalsCircle(1,i), DnormalsCircle(2,i))
    end do

    ! =================================
    ! Print statistics
    ! =================================

    call output_line ("Domain 'bench1.prm'. Statistics.")
    call output_line ("--------------------------------")
    
    call output_line ("Number of boundary components     : "//&
        trim(sys_siL(nbct,10)))
        
    call output_lbrk()

    call output_line ("Number of boundary segments, BC1  : "//&
        trim(sys_siL(nsega,10)))

    call output_line ("Max. parameter value (0-1), BC1   : "//&
        trim(sys_sdL(dmaxpar01a,3)))

    call output_line ("Max. parameter value (length), BC1: "//&
        trim(sys_sdL(dmaxparLena,3)))
        
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

    call output_lbrk()

    call output_line ("Number of boundary segments, BC2  : "//&
        trim(sys_siL(nsegb,10)))

    call output_line ("Max. parameter value (0-1), BC2   : "//&
        trim(sys_sdL(dmaxpar01b,3)))

    call output_line ("Max. parameter value (length), BC2: "//&
        trim(sys_sdL(dmaxparLenb,3)))
        
    call output_lbrk()
    call output_line ("Points on the circle:")

    do i=1,8
      ! Parameter value, 0-1 parametrisation
      dpar = real(i-1,DP)/8.0_DP
      
      ! Corresponding parameter value in length parametrisation
      dparlen = boundary_convertParameter(rboundary, 2, dpar, BDR_PAR_01, BDR_PAR_LENGTH)
      
      call output_line (" Parameter value "//&
          trim(sys_sdL(dpar,3))//"="//trim(sys_sdL(dparlen,3))//"      : "//&
          trim(sys_sdL(DcoordsCircle(1,i),3)) // ", " // trim(sys_sdL(DcoordsCircle(2,i),3)))
    end do

    call output_lbrk()
    call output_line ("Corresponding normal vectors:")

    do i=1,8
      ! Parameter value, 0-1 parametrisation
      dpar = real(i-1,DP)/8.0_DP
      
      ! Corresponding parameter value in length parametrisation
      dparlen = boundary_convertParameter(rboundary, 2, dpar, BDR_PAR_01, BDR_PAR_LENGTH)
      
      call output_line (" Parameter value "//&
          trim(sys_sdL(dpar,3))//"="//trim(sys_sdL(dparlen,3))//"      : "//&
          trim(sys_sdL(DnormalsCircle(1,i),3)) // ", " // trim(sys_sdL(DnormalsCircle(2,i),3)))
    end do

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the boundary definition
    call boundary_release(rboundary)
    
  end subroutine

end module

!##############################################################################
!# Tutorial 013c: Complete analysis of a boundary
!##############################################################################

module tutorial013c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use boundary

  implicit none
  private
  
  public :: start_tutorial013c

contains

  ! ***************************************************************************

  subroutine start_tutorial013c

    ! Declare some variables.
    type(t_boundary) :: rboundary
    
    integer :: nbct,nseg,ibct,iseg
    real(DP) :: dpar1, dparlen1,dpar2,dparlen2,dlen,dlenlen
    real(DP) :: dx,dy,dnx,dny
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 013c")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Read the underlying domain
    ! =================================

    call boundary_read_prm(rboundary, "pre/bench1.prm")

    ! =================================
    ! Gather / print statistics
    ! =================================

    ! ---------------------------------
    ! Global statistics
    ! ---------------------------------
    nbct = boundary_igetNBoundComp(rboundary)
    
    call output_line ("Domain 'bench1.prm'. Statistics.")
    call output_line ("--------------------------------")
    
    call output_line ("Number of boundary components     : "//&
        trim(sys_siL(nbct,10)))
        
    call output_lbrk()

    ! -------------------------------------------
    ! Loop through the boundary components    
    ! -------------------------------------------
    do ibct = 1,nbct

      ! ---------------------------------
      ! Boundary component statistics
      ! ---------------------------------
      
      ! Length, in 0-1 and length parametrisation
      dlen = boundary_dgetMaxParVal(rboundary, ibct)
      dlenlen = boundary_dgetMaxParVal(rboundary, ibct, BDR_PAR_LENGTH)
      
      ! Number of boundary segments
      nseg = boundary_igetNsegments (rboundary,ibct)

      call output_lbrk()
      call output_line ("  BC " // trim(sys_siL(ibct,10)) // &
          " - Length (0-1 par.)               : "//trim(sys_sdL(dlen,10)))
      call output_line ("  BC " // trim(sys_siL(ibct,10)) // &
          " - Length (length par.)            : "//trim(sys_sdL(dlenlen,10)))
      call output_line ("  BC " // trim(sys_siL(ibct,10)) // &
          " - Number of boundary segments     : "//trim(sys_siL(nseg,10)))
          
      ! -------------------------------------------
      ! Loop through the boundary segments    
      ! -------------------------------------------
      do iseg = 1,nseg

        call output_lbrk()
        call output_line ("  BC " // trim(sys_siL(ibct,10)) // &
            ", segment " // trim(sys_siL(iseg,10)) // " - Statistics:")
      
        ! -------------------------------------------
        ! Create a boundary region that describes this segment
        call boundary_createRegion (rboundary, ibct, iseg, rboundaryRegion)
        
        ! -------------------------------------------
        ! Print statistics
        ! -------------------------------------------
        
        ! -------------------------------------------
        ! Min/Max parameter value
        dpar1 = rboundaryRegion%dminParam
        dpar2 = rboundaryRegion%dmaxParam
        
        call output_line ("    Start/End parameter (0-1 par.)      : "//&
            trim(sys_sdL(dpar1,3)) // " / " // trim(sys_sdL(dpar2,3)) )
            
        ! -------------------------------------------
        ! Min/max par value, length parametrisation
        dparlen1 = boundary_convertParameter(rboundary, ibct, dpar1, BDR_PAR_01, BDR_PAR_LENGTH)
        dparlen2 = boundary_convertParameter(rboundary, ibct, dpar2, BDR_PAR_01, BDR_PAR_LENGTH)
        
        call output_line ("    Start/End parameter (length par.)   : "//&
            trim(sys_sdL(dparlen1,3)) // " / " // trim(sys_sdL(dparlen2,3)) )
        
        ! -------------------------------------------
        ! Length of the segment
        dlen = boundary_getRegionLength (rboundary,rboundaryRegion)
        call output_line ("    Length                              : "//&
            trim(sys_sdL(dlen,3)))
        
        ! -------------------------------------------
        ! Start point
        call boundary_getCoords(rboundary,ibct,dpar1,dx,dy)

        call output_line ("    Start point                         : "//&
            trim(sys_sdL(dx,3)) // " / " // trim(sys_sdL(dy,3)) )
  
        ! -------------------------------------------
        ! End point
        call boundary_getCoords(rboundary,ibct,dpar2,dx,dy)

        call output_line ("    End point                           : "//&
            trim(sys_sdL(dx,3)) // " / " // trim(sys_sdL(dy,3)) )
        
        ! -------------------------------------------
        ! Midpoint
        call boundary_getCoords(rboundary,ibct,0.5_DP*(dpar1+dpar2),dx,dy)

        call output_line ("    Midpoint                            : "//&
            trim(sys_sdL(dx,3)) // " / " // trim(sys_sdL(dy,3)) )

        ! -------------------------------------------
        ! Normal vector in the midpoint
        call boundary_getNormalVec2D(rboundary,ibct,0.5_DP*(dpar1+dpar2),dnx,dny)

        call output_line ("    Normal vector in the midpoint       : "//&
            trim(sys_sdL(dnx,3)) // " / " // trim(sys_sdL(dny,3)) )

      end do ! iseg
      
    end do ! iccd
    
    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the boundary definition
    call boundary_release(rboundary)
    
  end subroutine

end module

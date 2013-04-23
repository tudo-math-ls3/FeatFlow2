!##############################################################################
!# Tutorial 005c: Basic linear algebra
!##############################################################################

module tutorial005c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use linearalgebra

  implicit none
  private
  
  public :: start_tutorial005c

contains

  ! ***************************************************************************

  subroutine start_tutorial005c
  
    ! Declare some variables
    real(DP), dimension(4) :: Dvec1, Dvec2, Dvec3, Dvec4, Dvec
    integer, dimension(4) :: Iperm
    real(DP) :: d
    integer :: i

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 005c")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Initialise vectors
    ! =================================
    
    Dvec1 = (/  1.0_DP, 2.0_DP,  3.0_DP, 4.0_DP /)
    Dvec2 = (/ -1.0_DP, 1.0_DP, -1.0_DP, 1.0_DP /)
    Dvec3 = (/  0.0_DP, 0.0_DP,  1.0_DP, 0.0_DP /)
    Dvec4 = (/  2.0_DP,-1.0_DP,  0.0_DP, 0.0_DP /)
    
    ! =================================
    ! Clear
    !   Dvec = 0
    ! =================================
    
    call lalg_clearVector (Dvec)

    do i=1,4
      call output_line (" " // sys_sdL(Dvec(i),2), bnolinebreak=(i .ne. 4) )
    end do

    ! =================================
    ! Set to constant
    !   Dvec = (1,1,1,1)
    ! =================================
    
    call lalg_setVector (Dvec,1.0_DP)

    do i=1,4
      call output_line (" " // sys_sdL(Dvec(i),2), bnolinebreak=(i .ne. 4) )
    end do

    ! =================================
    ! Copy and scale:
    !   Dvec = 2*Dvec1
    ! =================================
    
    call lalg_copyVector (Dvec1, Dvec)
    call lalg_scaleVector (Dvec, 2.0_DP)

    do i=1,4
      call output_line (" " // sys_sdL(Dvec(i),2), bnolinebreak=(i .ne. 4) )
    end do
    
    ! =================================
    ! Linear combination.
    !   Dvec = Dvec1 - (1,1,1,1)
    ! =================================
    
    call lalg_setVector (Dvec,1.0_DP)
    call lalg_vectorLinearComb (Dvec1,Dvec,1.0_DP,-1.0_DP)

    do i=1,4
      call output_line (" " // sys_sdL(Dvec(i),2), bnolinebreak=(i .ne. 4) )
    end do

    ! =================================
    ! Add constant
    !   Dvec = Dvec1 - (1,1,1,1)
    ! =================================
    
    call lalg_copyVector (Dvec1,Dvec)
    call lalg_vectorAddScalar (Dvec,-1.0_DP)

    do i=1,4
      call output_line (" " // sys_sdL(Dvec(i),2), bnolinebreak=(i .ne. 4) )
    end do

    ! =================================
    ! Sort vector backwards,
    !   Dvec = back(Dvec4)
    ! =================================
    
    Iperm = (/ 4, 3, 2, 1 /)
    call lalg_vectorSort (Dvec4, Dvec, Iperm)

    do i=1,4
      call output_line (" " // sys_sdL(Dvec(i),2), bnolinebreak=(i .ne. 4) )
    end do

    ! =================================
    ! Scalar product
    !   result = < Dvec2, Dvec3 >
    ! =================================
    
    d = lalg_scalarProduct(Dvec2,Dvec3)
    
    call output_line ("<vec2,vec3> = " // sys_sdL(d,2) )

    ! =================================
    ! Different norms
    !   result = || Dvec4 || 
    ! =================================
    
    ! SUM
    d = lalg_norm(Dvec4,LINALG_NORMSUM)
    call output_line ("||vec4||_SUM   = " // sys_sdL(d,2) )

    ! Euclid-norm
    d = lalg_norm(Dvec4,LINALG_NORMEUCLID)
    call output_line ("||vec4||_2     = " // sys_sdL(d,2) )

    ! L1-norm
    d = lalg_norm(Dvec4,LINALG_NORML1)
    call output_line ("||vec4||_l1    = " // sys_sdL(d,2) )

    ! L2-norm
    d = lalg_norm(Dvec4,LINALG_NORML2)
    call output_line ("||vec4||_l2    = " // sys_sdL(d,2) )

    ! MAX-norm
    d = lalg_norm(Dvec4,LINALG_NORMMAX)
    call output_line ("||vec4||_MAX   = " // sys_sdL(d,2) )

  end subroutine

end module

!##############################################################################
!# Tutorial 005a: Basic math
!##############################################################################

module tutorial005a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use mprimitives

  implicit none
  private
  
  public :: start_tutorial005a

contains

  ! ***************************************************************************

  subroutine start_tutorial005a
  
    ! Declare some variables
    real(DP) :: dval1, dval2, dval3

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("Hello world. This is FEAT-2. Tutorial 005a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Kronecker Symbol
    ! =================================
    
    call output_line ("kronecker(3,5) = " // trim(sys_siL(mprim_kronecker(3,5),10)) )
    call output_line ("kronecker(5,5) = " // trim(sys_siL(mprim_kronecker(5,5),10)) )
    
    ! =================================
    ! Signum
    ! =================================
    
    call output_line ("signum(-5) = " // trim(sys_siL(mprim_signum(-5),10)) )
    call output_line ("signum( 0) = " // trim(sys_siL(mprim_signum( 0),10)) )
    call output_line ("signum( 2) = " // trim(sys_siL(mprim_signum( 2),10)) )

    ! =================================
    ! Linear rescaling
    ! =================================
    
    ! Linear interpolation of [0,1] to [5,7]
    call mprim_linearRescale(0.00_DP,0.0_DP,1.0_DP,5.0_DP,7.0_DP,dval1)
    call mprim_linearRescale(0.25_DP,0.0_DP,1.0_DP,5.0_DP,7.0_DP,dval2)
    call mprim_linearRescale(1.00_DP,0.0_DP,1.0_DP,5.0_DP,7.0_DP,dval3)
    
    call output_line ("linrescale(0.00,[0,1],[5,7]) = " // trim(sys_sdL(dval1,10)) )
    call output_line ("linrescale(0.25,[0,1],[5,7]) = " // trim(sys_sdL(dval2,10)) )
    call output_line ("linrescale(1.00,[0,1],[5,7]) = " // trim(sys_sdL(dval3,10)) )

    ! =================================
    ! Quadratic interpolation
    ! =================================
    
    ! Quadratic interpolation, i.e. evaluate a quadratic
    ! polynomial p with p(-1)=2, p(0)=0, p(1)=2
    call mprim_quadraticInterpolation (-1.0_DP,2.0_DP,0.0_DP,2.0_DP,dval1)
    call mprim_quadraticInterpolation ( 0.0_DP,2.0_DP,0.0_DP,2.0_DP,dval2)
    call mprim_quadraticInterpolation ( 0.5_DP,2.0_DP,0.0_DP,2.0_DP,dval3)
    
    call output_line ("quadinterpol(-1.0,[2,0,2]) = " // trim(sys_sdL(dval1,10)) )
    call output_line ("quadinterpol( 0.0,[2,0,2]) = " // trim(sys_sdL(dval2,10)) )
    call output_line ("quadinterpol( 0.5,[2,0,2]) = " // trim(sys_sdL(dval3,10)) )
    
    ! =================================
    ! Parabolic profile
    ! =================================
    
    ! This is a quadratic polynomial p with p(0)=p(len)=0, p(mid)=max.
    dval1 = mprim_getParabolicProfile (0.5_DP,2.0_DP,10.0_DP)
    dval2 = mprim_getParabolicProfile (1.0_DP,2.0_DP,10.0_DP)
    dval3 = mprim_getParabolicProfile (1.5_DP,2.0_DP,10.0_DP)
    
    call output_line ("parprofile(0.5,len=2,max=10) = " // trim(sys_sdL(dval1,10)) )
    call output_line ("parprofile(1.0,len=2,max=10) = " // trim(sys_sdL(dval2,10)) )
    call output_line ("parprofile(1.5,len=2,max=10) = " // trim(sys_sdL(dval3,10)) )

  end subroutine

end module

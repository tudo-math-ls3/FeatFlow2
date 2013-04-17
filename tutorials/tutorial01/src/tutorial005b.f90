!##############################################################################
!# Tutorial 005b: Working with small, dense matrices
!##############################################################################

module tutorial005b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use mprimitives

  implicit none
  private
  
  public :: start_tutorial005b

contains

  ! ***************************************************************************

  subroutine start_tutorial005b
  
    ! Declare some variables
    real(DP), dimension(4,4) :: DaTrans, Da, DaInv
    integer :: i,j
    logical :: bsuccess

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 005b")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Transpose a matrix
    ! =================================
    
    ! Initialise a matrix - transposed, as Fortran works column-wise
    DaTrans = reshape ( (/ 2.0_DP, -1.0_DP, 0.0_DP, 0.0_DP, &
                           0.0_DP,  2.0_DP,-1.0_DP, 0.0_DP, &
                           0.0_DP,  0.0_DP, 2.0_DP,-1.0_DP, &
                           0.0_DP,  0.0_DP, 0.0_DP, 2.0_DP /), (/ 4, 4 /) )
                       
    ! Transpose the matrix.
    ! This is an optimised variant of the TRANSPOSE command of Fortran.
    call mprim_transposeMatrix(DaTrans,Da)
    
    ! =================================
    ! Invert the matrix
    ! =================================
    
    ! Invert the matrix
    DaInv(:,:) = Da(:,:)
    call mprim_invertMatrixPivot(DaInv,4,bsuccess)
    
    if (.not. bsuccess) then
      call sys_halt()
    end if
    
    ! =================================
    ! Output
    ! =================================

    call output_line ("Transposed matrix:")
    do i=1,4
      do j=1,4
        call output_line (sys_adjustr(sys_sdL(DaTrans(i,j),2),6), bnolinebreak=(j .ne. 4) )
      end do
    end do

    call output_lbrk()
    call output_line ("Matrix:")
    do i=1,4
      do j=1,4
        call output_line (sys_adjustr(sys_sdL(Da(i,j),2),6), bnolinebreak=(j .ne. 4) )
      end do
    end do

    call output_lbrk()
    call output_line ("Its inverse:")

    do i=1,4
      do j=1,4
        call output_line (" " // sys_sdL(DaInv(i,j),2), bnolinebreak=(j .ne. 4) )
      end do
    end do
    
  end subroutine

end module

! This file describes the interfaces of the callback routines
! of the discrete operators for scalar conservation laws

! ***************************************************************************

interface

! ***************************************************************************

!<subroutine>

pure subroutine fcb_calcConvection(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the matrix coefficients $k_{ij}$ and $k_{ji}$
    ! resulting from the discretisation of the convective operator.
!</description>

use fsystem
use dofmapping

!<input>
  real(DP),               intent(IN)  :: u_i,u_j
  real(DP), dimension(:), intent(IN)  :: C_ij,C_ji
  integer(PREC_DOFIDX),   intent(IN)  :: i,j
!</input>

!<output>
  real(DP), intent(OUT) :: k_ij,k_ji
!</output>

!</subroutine>

end subroutine fcb_calcConvection

end interface

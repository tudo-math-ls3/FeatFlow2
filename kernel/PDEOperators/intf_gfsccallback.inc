  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! of the discrete operators for scalar conservation laws
  
  ! ***************************************************************************
  
  interface
  
  ! ***************************************************************************

!<subroutine>

    pure subroutine fcb_calcMatrix(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

    use fsystem

!<description>
    ! This subroutine computes the local matrices for the edge
    ! $ij$. Given the conservative variables at the two states $i$ and
    ! $j$ and the multi-component weighting coefficients $C_{ij}$ and
    ! $C_{ji}$, this subroutine computes the matrix coefficients
    ! $k_{ij} and $k_{ji} which are assembled into the global system
    ! matrix. Moreover, the diffusion coefficient $d_{ij}$ is computed.
!</description>

!<input>
  real(DP), intent(IN) :: u_i, u_j
  real(DP), dimension(:), intent(IN)  :: C_ij, C_ji
  integer, intent(IN) :: i,j
!</input>

!<outptu>
  real(DP), intent(OUT) :: k_ij, k_ji, d_ij
!</output>

!</subroutine>

  end subroutine fcb_calcMatrix

end interface

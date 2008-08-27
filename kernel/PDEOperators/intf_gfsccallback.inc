! This file describes the interfaces of the callback routines
! of the discrete operators for scalar conservation laws

! ***************************************************************************

INTERFACE

! ***************************************************************************

!<subroutine>

PURE SUBROUTINE fcb_calcConvection(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, istatus)

!<description>
    ! This subroutine computes the matrix coefficients $k_{ij}$ and $k_{ji}$
    ! resulting from the discretisation of the convective operator.
!</description>

USE fsystem
USE dofmapping

!<input>
  REAL(DP),               INTENT(IN)  :: u_i,u_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: C_ij,C_ji
  INTEGER(PREC_DOFIDX),   INTENT(IN)  :: i,j
!</input>

!<inputoutput>
  INTEGER, INTENT(INOUT), OPTIONAL    :: istatus
!</inputoutput>

!<output>
  REAL(DP), INTENT(OUT)               :: k_ij,k_ji
!</output>

!</subroutine>

END SUBROUTINE fcb_calcConvection

END INTERFACE
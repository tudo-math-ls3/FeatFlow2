INTERFACE

!<subroutine>

SUBROUTINE fcb_calcVelocity(u_i,u_j,i,j,v_ij,v_ji)

!<description>
    ! This subroutine computes the (characteristic) velocity
    ! for the edge that connects nodes I and J.
!</description>

USE fsystem
USE dofmapping

!<input>
  REAL(DP), INTENT(IN)                :: u_i,u_j
  INTEGER(PREC_DOFIDX),  INTENT(IN)   :: i,j
!</input>

!<output>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: v_ij,v_ji
!</output>

!</subroutine>

END SUBROUTINE fcb_calcVelocity

END INTERFACE

	INTERFACE

	SUBROUTINE fcb_getVelocity(u_i,u_j,i,j,v_ij,v_ji)

	USE fsystem

	REAL(DP), INTENT(IN)                :: u_i,u_j
	INTEGER, INTENT(IN)                 :: i,j
	REAL(DP), DIMENSION(:), INTENT(OUT) :: v_ij,v_ji

	END SUBROUTINE fcb_getVelocity

	END INTERFACE
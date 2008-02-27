! This file describes the interfaces of the callback routines
! of the discrete operators for hyperbolic systems

! ***************************************************************************

INTERFACE

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcCharacteristics(U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

USE fsystem

!<description>
    ! This subroutine computes the characteristic transformation variables.
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dweight
!</input>

!<output>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: W_ij
  REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: Lbd_ij
  REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: R_ij
  REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: L_ij
!</output>

!</subroutine>

END SUBROUTINE fcb_calcCharacteristics

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcFlux(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji)

USE fsystem

!<description>
    ! This subroutine computes the internodal fluxes.
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: C_ij,C_ji
  REAL(DP),               INTENT(IN)  :: dscale
!</input>

!<output>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: F_ij,F_ji
!</output>

!</subroutine>

END SUBROUTINE fcb_calcFlux

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcMatrix(U_i, U_j, C_ij, C_ji, A_ij, S_ij)

USE fsystem

!<description>
    ! This subroutine computes the local matrices.
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: C_ij,C_ji
!</input>

!<outptu>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: A_ij,S_ij
!</output>

!</subroutine>

END SUBROUTINE fcb_calcMatrix

END INTERFACE
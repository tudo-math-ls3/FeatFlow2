! This file describes the interfaces of the callback routines
! of the discrete operators for hyperbolic systems

! ***************************************************************************
INTERFACE

SUBROUTINE fcb_getRoeMatrix(U_i, U_j, Dweight, A_ij)
USE fsystem

REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
REAL(DP), DIMENSION(:), INTENT(IN)  :: Dweight
REAL(DP), DIMENSION(:), INTENT(OUT) :: A_ij

END SUBROUTINE fcb_getRoeMatrix

SUBROUTINE fcb_getDissipation(U_i, U_j, Dweight, D_ij)

USE fsystem

REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
REAL(DP), DIMENSION(:), INTENT(IN)  :: Dweight
REAL(DP), DIMENSION(:), INTENT(OUT) :: D_ij

END SUBROUTINE fcb_getDissipation


SUBROUTINE fcb_getCharacteristics(U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

USE fsystem

REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
REAL(DP), DIMENSION(:), INTENT(IN)  :: Dweight
REAL(DP), DIMENSION(:), INTENT(OUT) :: W_ij
REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: Lbd_ij
REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: R_ij
REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: L_ij
END SUBROUTINE fcb_getCharacteristics

END INTERFACE
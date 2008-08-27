! This file describes the interfaces of the callback routines
! of the discrete operators for hyperbolic systems

! ***************************************************************************

INTERFACE

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcCharacteristics(U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij, istatus)

USE fsystem

!<description>
    ! This subroutine computes the characteristic transformation variables.
    ! Given the conservative variables at the two states $i$ and $j$ and the
    ! multi-component weighting coefficient, this routine computes:
    !
    ! 1.) The characteristic solution difference $W_{ij}=R_{ij}(U_j-U_i)$
    !
    ! 2.) The diagonal matrix of real eigenvalues $\Lambda_{ij}$
    !
    ! 3.) The matrix of right eigenvectors $R_{ij}$
    !
    ! 4.) The matrix of left eigenvectors $L_{ij}$
    !
    ! For hyperbolic systems, the following characteristic decomposition holds
    !
    ! $$A_{ij} = L_{ij} \Lambda_{ij} R_{ij}$$
    !
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Dweight
!</input>

!<inputoutput>
  INTEGER, INTENT(INOUT), OPTIONAL    :: istatus
!</inputoutput>

!<output>
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: W_ij
  REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: Lbd_ij
  REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: R_ij
  REAL(DP), DIMENSION(:), INTENT(OUT), OPTIONAL :: L_ij
!</output>

!</subroutine>

END SUBROUTINE fcb_calcCharacteristics

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcFlux(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji, istatus)

USE fsystem

!<description>
    ! This subroutine computes the internodal fluxes.
    ! Given the conservative variables at the two states $i$ and $j$ and the
    ! multi-component weighting coefficients, this routine computes:
    !
    ! 1.) The internodal flux from node $j$ into node $i$: $F_{ij}$
    !
    ! 2.) The internodal flux from node $i$ into node $j$: $F_{ji}$
    !
    ! In some cases it is useful to have skew-symmetry, i.e. $F_{ji}=-F_{ij}$
    ! but this is not mandatory, and hence, this routine computes both.
    !
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: C_ij,C_ji
  REAL(DP),               INTENT(IN)  :: dscale
!</input>

!<inputoutput>
  INTEGER, INTENT(INOUT), OPTIONAL    :: istatus
!</inputoutput>

!<output>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: F_ij,F_ji
!</output>

!</subroutine>

END SUBROUTINE fcb_calcFlux

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcRawFlux(U_i, U_j, C_ij, C_ji, dscale, F_ij, istatus)

USE fsystem

!<description>
    ! This subroutine computes the internodal fluxes.
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: C_ij,C_ji
  REAL(DP),               INTENT(IN)  :: dscale
!</input>

!<inputoutput>
  INTEGER, INTENT(INOUT), OPTIONAL    :: istatus
!</inputoutput>

!<output>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: F_ij
!</output>

!</subroutine>

END SUBROUTINE fcb_calcRawFlux

! ***************************************************************************

!<subroutine>

SUBROUTINE fcb_calcMatrix(U_i, U_j, C_ij, C_ji, A_ij, S_ij, istatus)

USE fsystem

!<description>
    ! This subroutine computes the local matrices.
    ! Given the conservative variables at the two states $i$ and $j$ and the
    ! multi-component weighting coefficients, this routine computes:
    !
    ! 1.) The skew-symmetric part of the local system matrix, i.e. $A_{ji}=-A_{ij}$
    !
    ! 2.) The symmetric part of the local system matrix, i.e. $S_{ji}=S_{ij}$
    !
!</description>

!<input>
  REAL(DP), DIMENSION(:), INTENT(IN)  :: U_i,U_j
  REAL(DP), DIMENSION(:), INTENT(IN)  :: C_ij,C_ji
!</input>

!<inputoutput>
  INTEGER, INTENT(INOUT), OPTIONAL    :: istatus
!</inputoutput>

!<outptu>
  REAL(DP), DIMENSION(:), INTENT(OUT) :: A_ij,S_ij
!</output>

!</subroutine>

END SUBROUTINE fcb_calcMatrix

END INTERFACE
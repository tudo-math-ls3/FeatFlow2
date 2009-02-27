  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! of the discrete operators for hyperbolic systems

  ! ***************************************************************************

  interface

    ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcCharacteristics(U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij, istatus)

    use fsystem

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
  real(DP), dimension(:), intent(IN)  :: U_i,U_j
  real(DP), dimension(:), intent(IN)  :: Dweight
!</input>

!<inputoutput>
  integer, intent(INOUT), optional    :: istatus
!</inputoutput>

!<output>
  real(DP), dimension(:), intent(OUT)           :: W_ij
  real(DP), dimension(:), intent(OUT), optional :: Lbd_ij
  real(DP), dimension(:), intent(OUT), optional :: R_ij
  real(DP), dimension(:), intent(OUT), optional :: L_ij
!</output>

!</subroutine>

  end subroutine fcb_calcCharacteristics

! ***************************************************************************

!<subroutine>

  subroutine fcb_calcFlux(U_i, U_j, C_ij, C_ji, dscale, F_ij, F_ji, istatus)

  use fsystem

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
  real(DP), dimension(:), intent(IN)  :: U_i,U_j
  real(DP), dimension(:), intent(IN)  :: C_ij,C_ji
  real(DP),               intent(IN)  :: dscale
!</input>

!<inputoutput>
  integer, intent(INOUT), optional    :: istatus
!</inputoutput>

!<output>
  real(DP), dimension(:), intent(OUT) :: F_ij,F_ji
!</output>

!</subroutine>

  end subroutine fcb_calcFlux

! ***************************************************************************

!<subroutine>

  subroutine fcb_calcRawFlux(U_i, U_j, C_ij, C_ji, dscale, F_ij, istatus)

  use fsystem

!<description>
    ! This subroutine computes the internodal fluxes.
!</description>

!<input>
  real(DP), dimension(:), intent(IN)  :: U_i,U_j
  real(DP), dimension(:), intent(IN)  :: C_ij,C_ji
  real(DP),               intent(IN)  :: dscale
!</input>

!<inputoutput>
  integer, intent(INOUT), optional    :: istatus
!</inputoutput>

!<output>
  real(DP), dimension(:), intent(OUT) :: F_ij
!</output>

!</subroutine>

  end subroutine fcb_calcRawFlux

! ***************************************************************************

!<subroutine>

  subroutine fcb_calcMatrix(U_i, U_j, C_ij, C_ji, A_ij, S_ij, istatus)

  use fsystem

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
  real(DP), dimension(:), intent(IN)  :: U_i,U_j
  real(DP), dimension(:), intent(IN)  :: C_ij,C_ji
!</input>

!<inputoutput>
  integer, intent(INOUT), optional    :: istatus
!</inputoutput>

!<outptu>
  real(DP), dimension(:), intent(OUT) :: A_ij,S_ij
!</output>

!</subroutine>

  end subroutine fcb_calcMatrix

end interface

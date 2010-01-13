  ! -*- mode: f90; -*-
  
  ! This file describes the interfaces of the callback routines
  ! of the discrete operators for hyperbolic systems

  ! ***************************************************************************

  interface

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcMatrixDiagonal(U_i, C_ii, i, dscale, K_ii)

    use fsystem

!<description>
    ! Given the variables at node $i$ and the weighting coefficient
    ! $C_{ii}$, this subroutine computes the matrix coefficient
    ! $K_{ii} which is assembled into the diagonal of the global
    ! system matrix.
!</description>

!<input>
  ! Variable at state $i$
  real(DP), dimension(:), intent(in) :: U_i

  ! Multi-component weighting coefficient
  real(DP), dimension(:), intent(in) :: C_ii

  ! Scaling parameter
  real(DP), intent(in) :: dscale

  ! State $i$
  integer, intent(in) :: i
!</input>

!<outptu>
  ! Matrix coefficient
  real(DP), dimension(:), intent(out) :: K_ii
!</output>

!</subroutine>

  end subroutine fcb_calcMatrixDiagonal

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcMatrix(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)

    use fsystem

!<description>
    ! This subroutine computes the local matrices for the edge
    ! $ij$. Given the variables at the two states $i$ and $j$ and the
    ! multi-component weighting coefficients $C_{ij}$ and $C_{ji}$,
    ! this subroutine computes the matrix coefficients $K_{ij} and
    ! $K_{ji} which are assembled into the global system
    ! matrix. Moreover, the diffusion coefficient $D_{ij}$ is
    ! computed.
!</description>

!<input>
  ! Variables at states $i$ and $j$
  real(DP), dimension(:), intent(in) :: U_i, U_j

  ! Multi-component weighting coefficients
  real(DP), dimension(:), intent(in) :: C_ij, C_ji

  ! Scaling parameter
  real(DP), intent(in) :: dscale

  ! States $i$ and $j$
  integer, intent(in) :: i,j
!</input>

!<outptu>
  ! Matrix coefficients
  real(DP), dimension(:), intent(out) :: K_ij, K_ji, D_ij
!</output>

!</subroutine>

  end subroutine fcb_calcMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcFlux(&
      U_i, U_j, C_ij, C_ji, i, j, dscale, F_ij, F_ji)

    use fsystem

!<description>
    ! This subroutine computes the internodal fluxes.  Given the
    ! variables at the two states $i$ and $j$ and the multi-component
    ! weighting coefficients, this routine computes:
    !
    ! 1.) The internodal flux from node $j$ into node $i$: $F_{ij}$
    !
    ! 2.) The internodal flux from node $i$ into node $j$: $F_{ji}$
    !
    ! In some cases it is useful to have skew-symmetry,
    ! i.e. $F_{ji}=-F_{ij}$ but this is not mandatory, and hence, this
    ! routine computes both.
!</description>

!<input>
  ! Variables at states $i$ and $j$
  real(DP), dimension(:), intent(in) :: U_i, U_j

  ! Multi-component weighting coefficients
  real(DP), dimension(:), intent(in) :: C_ij, C_ji

  ! Scaling parameter
  real(DP), intent(in) :: dscale

  ! States $i$ and $j$
  integer, intent(in) :: i,j
!</input>

!<output>
  ! Fluxes
  real(DP), dimension(:), intent(out) :: F_ij, F_ji
!</output>

!</subroutine>

  end subroutine fcb_calcFlux

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcFluxFCT(U1_i, U1_j, U2_i, U2_j,&
      C_ij, C_ji, i, j, dscale1, dscale2, F_ij)

    use fsystem

!<description>
    ! This subroutine computes the internodal fluxes for the FCT algorithm.
    ! Given the variables at the two states $i$ and $j$ this routine computes:
    !
    ! The internodal flux from node $j$ into node $i$: $F_{ij}=-F_{ji}$
!</description>

!<input>
  ! Variables at states $i$ and $j$
  real(DP), dimension(:), intent(in) :: U1_i, U1_j
  real(DP), dimension(:), intent(in) :: U2_i, U2_j

  ! Multi-component weighting coefficients
  real(DP), dimension(:), intent(in) :: C_ij, C_ji

  ! Scaling parameters
  real(DP), intent(in) :: dscale1, dscale2

  ! States $i$ and $j$
  integer, intent(in) :: i,j
!</input>

!<output>
  ! raw antidiffusive flux from node $j$ into node $i$
  real(DP), dimension(:), intent(out) :: F_ij
!</output>

!</subroutine>

  end subroutine fcb_calcFluxFCT

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcTransformation(U_i, U_j, F_ij, G_ij, G_ji)

    use fsystem

!<description>
    ! This subroutine computes the transformed fluxes.
!</description>

!<input>
  ! Variables at states $i$ and $j$
  real(DP), dimension(:), intent(in) :: U_i, U_j

  ! Raw antidiffusive flux from node $j$ into node $i$
  real(DP), dimension(:), intent(in) :: F_ij
!</input>

!<output>
  ! Transformed raw antidiffusive fluxes
  real(DP), dimension(:), intent(out) :: G_ij,G_ji
  !</output>

!</subroutine>

  end subroutine fcb_calcTransformation

  ! ***************************************************************************

!<subroutine>

  subroutine fcb_calcCharacteristics(U_i, U_j, Dweight, W_ij, Lbd_ij, R_ij, L_ij)

    use fsystem

!<description>
    ! This subroutine computes the characteristic transformation
    ! variables.  Given the conservative variables at the two states
    ! $i$ and $j$ and the multi-component weighting coefficient, this
    ! routine computes:
    !
    ! 1.) The characteristic solution difference
    !     $W_{ij}=R_{ij}(U_j-U_i)$
    !
    ! 2.) The diagonal matrix of real eigenvalues $\Lambda_{ij}$
    !
    ! 3.) The matrix of right eigenvectors $R_{ij}$
    !
    ! 4.) The matrix of left eigenvectors $L_{ij}$
    !
    ! For hyperbolic systems, the following characteristic
    ! decomposition holds
    !
    ! $$A_{ij} = L_{ij} \Lambda_{ij} R_{ij}$$
!</description>

!<input>
  real(DP), dimension(:), intent(in)  :: U_i,U_j
  real(DP), dimension(:), intent(in)  :: Dweight
!</input>

!<output>
  real(DP), dimension(:), intent(out), optional :: W_ij
  real(DP), dimension(:), intent(out), optional :: Lbd_ij
  real(DP), dimension(:), intent(out), optional :: R_ij
  real(DP), dimension(:), intent(out), optional :: L_ij
!</output>

!</subroutine>

  end subroutine fcb_calcCharacteristics

  end interface

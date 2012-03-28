!##############################################################################
!# ****************************************************************************
!# <name> coarsegridcorrection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains different methods to calculate the coarse grid
!# parameter $\alpha$ in the multigrid solver. There are general methods
!# provided here as well as special calculation methods for specific
!# systems.
!#
!# The following routines can be found here:
!#
!# 1.) cgcor_init
!#     -> Initialise a coarse grid correction structure
!#
!# 2.) cgcor_release
!#     -> Releases a coarse grid correction structure
!#
!# 3.) cgcor_calcOptimalCorrection
!#     -> Calculate the optimal alpha for the correction.
!#
!#  FAQ - Frequently asked Questions  \\
!# ---------------------------------- \\
!# 1.) What is this p_DequationWeights-parameter in the t_coarseGridCorrection
!#     structure for? I cannot understand?
!#
!#   This parameter allows to multiply an equation by a given factor
!#   while calculating optimal correction factors. This is used when
!#   calculating residuals. Imagine you want to calculate the residual
!#   of a 2D Navier-Stokes system:
!#
!# <verb>
!#    (r1)     (r1)   (A11  A12  B1)  (x1)
!#    (r2)  =  (r2) - (A21  A22  B2)  (x2)
!#    (rp)     (rp)   (D1   D2     )  (xp)
!# </verb>
!#
!#   Now imagine that you set the array p_DequationWeights to (1,1,-1); then,
!#   residuals are calculated of the following system:
!#
!# <verb>
!#    (r1)     ( f1)   ( A11    A12  B1)  (x1)
!#    (r2)  =  ( f2) - ( A21    A22  B2)  (x2)
!#    (rp)     (-fp)   (-D1    -D2     )  (xp)
!# </verb>
!#
!#   i.e. the 3rd equation is multiplied by -1. This is just a small trick
!#   to get better convergence rates by enforcing/destroying symmetry in
!#   an operator. The modified residual gives then another correction
!#   factor for the coarse grid correction.
!#
!# 2.) How to use that array?
!#
!#   Set up a coarse grid correction structure with
!#
!#    call cgcor_init (rcoarseGridCorrection,3)
!#
!#   where "3" is in this example the number of equations in your PDE.
!#   Then modify the rcoarseGridCorrection%p_DequationWeights array manually;
!#   the array is predefined with 1.0.
!#
!# </purpose>
!##############################################################################

module coarsegridcorrection

!$use omp_lib
  use fsystem
  use linearsystemscalar
  use linearsystemblock
  use filtersupport

  implicit none

  private

!<constants>

!<constantblock description="Method identifier for coarse grid correction">

  ! Standard constant damping. A correction factor of 1.0 is used for
  ! the coarse grid correction if dalphaMin < 1.0 < dalphaMax.
  ! Otherwise, the constant dalphaMin is chosen.
  ! Remark: This is the optimal setting for a scalar equation with conformal
  !  elements; however, using nonconformal elements or nonscalar equations
  !  might make it necessary to switch to another method.
  integer, parameter, public :: CGCOR_STANDARD       = 0

  ! Damping by energy minimisation.
  ! Remark: This is the optimal setting for a scalar Laplace equation
  !  with nonconformal Rannacher-Turek element Ex30/Ex31.
  integer, parameter, public :: CGCOR_SCALARENERGYMIN = 1

  ! Damping by defect minimisation.
  ! Remark: Easy/cheap to calculate, but no functional analytic background.
  integer, parameter, public :: CGCOR_SCALARDEFMIN    = 2

!</constantblock>

!</constants>


!<types>

!<typeblock>

  ! Coarse grid correction structure; defines the behaviour of
  ! how to calculate the optimal coarse grid correction.
  type t_coarseGridCorrection

    ! Method identifier. This is one of the CGCOR_xxxx constants and specifies
    ! the algorithm that is used to calculate the optimal coarse grid
    ! correction.
    integer :: ccorrectionType = CGCOR_STANDARD

    ! Minimum damping parameter
    real(DP) :: dalphaMin = -10.0_DP

    ! Maximum damping parameter
    real(DP) :: dalphaMax = 10.0_DP

    ! A list of weights for the different equations.
    ! If not associated, a standard value of 1.0 is assumed for every equation.
    ! If associated, the calculation routine for the optimal
    ! coarse grid correction will multiply residuals of the corresponding
    ! equation with this factor.
    ! Example: 2D Navier-Stokes: (1,1,-1)
    ! -> The pressure equation is multiplied by -1 before taking the norm
    !    of the residual.
    ! Can be used to symmetrise equations.
    real(DP), dimension(:), pointer :: p_DequationWeights => null()

  end type

  public :: t_coarseGridCorrection

!</typeblock>

!</types>

  public :: cgcor_init
  public :: cgcor_release
  public :: cgcor_calcOptimalCorrection

  ! ***************************************************************************

contains

  ! ***************************************************************************

!<subroutine>

  subroutine cgcor_init (rcoarseGridCorrection,nequations)

!<description>
  ! Initialises a coarse grid correction structure.
!</description>

!<input>
  ! OPTIONAL: Number of equations.
  ! If specified, the p_DequationWeights array in the coarse grid correction
  ! structure is allocated and initially filled with 1.0-values.
  integer, intent(in), optional :: nequations
!</input>

!<output>
  ! A coarse grid correction structure to be initialised
  type(t_coarseGridCorrection), intent(out) :: rcoarseGridCorrection
!</output>

    ! Initialise by default initialisation
    !
    ! Exception: The p_DequationWeights array.
    if (present(nequations)) then
      allocate(rcoarseGridCorrection%p_DequationWeights(nequations))
      rcoarseGridCorrection%p_DequationWeights = 1.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cgcor_release (rcoarseGridCorrection)

!<description>
  ! Releases a coarse grid correction structure.
!</description>

!<inputoutput>
  ! A coarse grid correction structure to be released
  type(t_coarseGridCorrection), intent(inout) :: rcoarseGridCorrection
!</inputoutput>

    ! Release the pointer to the equation weights if associated
    if (associated(rcoarseGridCorrection%p_DequationWeights)) then
      deallocate(rcoarseGridCorrection%p_DequationWeights)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cgcor_calcOptimalCorrection (rcoarseGridCorrection,&
                                          rmatrix,rvector,rrhs,rcorrVector,&
                                          rtempVector,p_RfilterChain,dalpha)

!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively. The parameter ccorrectionType in the
  ! rcoarseGridCorrection structure defines the method how to calculate this
  ! parameter.
!</description>

!<input>
  ! A coarse grid correction structure specifying the algorithm to use.
  type(t_coarseGridCorrection), intent(in) :: rcoarseGridCorrection

  ! The block matrix of the system Ax=b which is solved by multigrid
  type(t_matrixBlock), intent(in) :: rmatrix

  ! The current (uncorrected) solution vector
  type(t_vectorBlock), intent(in) :: rvector

  ! The current RHS vector
  type(t_vectorBlock), intent(in) :: rrhs

  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=<tex>$P^{-1}(b-Ax)$</tex> and <tex>$P^{-1}=$</tex> multigrid on the
  ! coarse level.
  type(t_vectorBlock), intent(in) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
!</input>

!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  type(t_vectorBlock), intent(inout) :: rtempVector
!</inputoutput>

!<output>
  ! The optimal correction parameter <tex>$\alpha$</tex> for the coarse grid correction.
  real(DP), intent(out) :: dalpha
!</output>

!</subroutine>

  if (rcoarseGridCorrection%dalphaMax .le. rcoarseGridCorrection%dalphaMin) then
    dalpha = 1.0_DP
  else
    ! Which method to use?

    select case (rcoarseGridCorrection%ccorrectionType)
    case (CGCOR_SCALARENERGYMIN)
      if (.not. associated(rcoarseGridCorrection%p_DequationWeights)) then
        ! Standard energy minimisation
        call cgcor_calcCorrEnergyMin (rmatrix,rvector,rrhs,rcorrVector,&
                                      rtempVector,p_RfilterChain,dalpha)
      else
        ! Special energy minimisation with different weights for the
        ! different equations.
        call cgcor_calcCorrEnergyMinWeighted (rmatrix,rvector,rrhs,rcorrVector,&
            rtempVector,p_RfilterChain,dalpha,rcoarseGridCorrection%p_DequationWeights)
      end if
    case (CGCOR_SCALARDEFMIN)
      call cgcor_calcCorrDefMin (rmatrix,rvector,rrhs,rcorrVector,&
                                rtempVector,p_RfilterChain,dalpha)
    case DEFAULT !(=CGCOR_STANDARD)
      dalpha = 1.0_DP   ! Standard setting
    end select
  end if

  ! Make sure it is in the interval given by the dalphaMin/dalphaMax
  dalpha = max(min(dalpha,rcoarseGridCorrection%dalphaMax), &
                          rcoarseGridCorrection%dalphaMin)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cgcor_calcCorrEnergyMin (rmatrix,rvector,rrhs,rcorrVector,&
                                      rtempVecBlock,p_RfilterChain,dalpha)

!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
!</description>

!<input>
  ! The block matrix of the system Ax=b which is solved by multigrid
  type(t_matrixBlock), intent(in) :: rmatrix

  ! The current (uncorrected) solution vector
  type(t_vectorBlock), intent(in) :: rvector

  ! The current RHS vector
  type(t_vectorBlock), intent(in) :: rrhs

  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=<tex>$P^{-1}(b-Ax)$</tex> and <tex>$P^{-1}=$</tex> multigrid on the
  ! coarse level.
  type(t_vectorBlock), intent(in) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
!</input>

!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  type(t_vectorBlock), intent(inout) :: rtempVecBlock
!</inputoutput>

!<output>
  ! The optimal correction parameter <tex>$\alpha$</tex> for the coarse grid correction.
  real(DP), intent(out) :: dalpha
!</output>

!</subroutine>

    ! local variables

    real(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek`s book):
    !
    !             ( f_k - A_k x_k  ,  corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction

    call lsysbl_copyVector(rrhs,rtempVecBlock)
    call lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    if (associated(p_RfilterChain)) then
      call filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    end if

    a = lsysbl_scalarProduct(rtempVecBlock,rcorrVector)

    ! Calculate the demoninator of the fraction
    call lsysbl_blockMatVec(rmatrix, rcorrVector, rtempVecBlock, 1.0_DP,0.0_DP)
    ! Apply the filter
    if (associated(p_RfilterChain)) then
      call filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    end if

    b = lsysbl_scalarProduct(rtempVecBlock,rcorrVector)

    ! Return the alpha.
    if (b .ne. 0.0_DP) then
      dalpha = a/b
    else
      dalpha = 1.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cgcor_calcCorrEnergyMinWeighted (rmatrix,rvector,rrhs,rcorrVector,&
      rtempVecBlock,p_RfilterChain,dalpha,DequationWeights)

!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
  ! For every equation (line of the matrix), DequationWeights allows to
  ! specify a weight which is multiplied to residuum norms before
  ! summing up the residuum norms of all subvectors to a common one.
!</description>

!<input>
  ! The block matrix of the system Ax=b which is solved by multigrid
  type(t_matrixBlock), intent(in) :: rmatrix

  ! The current (uncorrected) solution vector
  type(t_vectorBlock), intent(in) :: rvector

  ! The current RHS vector
  type(t_vectorBlock), intent(in) :: rrhs

  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=<tex>$P^{-1}(b-Ax)$</tex> and <tex>$P^{-1}=$</tex> multigrid on the
  ! coarse level.
  type(t_vectorBlock), intent(in) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain

  ! List of weights for every equation (=row of the matrix)
  real(DP), dimension(:), intent(in) :: DequationWeights
!</input>

!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  type(t_vectorBlock), intent(inout) :: rtempVecBlock
!</inputoutput>

!<output>
  ! The optimal correction parameter <tex>$\alpha$</tex> for the coarse grid correction.
  real(DP), intent(out) :: dalpha
!</output>

!</subroutine>

    ! local variables

    real(DP) :: a,b
    integer :: irow

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek`s book):
    !
    !             ( f_k - A_k x_k  ,  corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction

    call lsysbl_copyVector(rrhs,rtempVecBlock)
    call lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    if (associated(p_RfilterChain)) then
      call filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    end if

    ! Calculate a weighted scalar product using the equation weights.
    a = 0.0_DP
    do irow = 1,rtempVecBlock%nblocks
      a = a + DequationWeights(irow) * &
          lsyssc_scalarProduct(rtempVecBlock%RvectorBlock(irow),&
                               rcorrVector%RvectorBlock(irow))
    end do

    ! Calculate the demoninator of the fraction
    call lsysbl_blockMatVec(rmatrix, rcorrVector, rtempVecBlock, 1.0_DP,0.0_DP)
    ! Apply the filter
    if (associated(p_RfilterChain)) then
      call filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    end if

    b = 0.0_DP
    do irow = 1,rtempVecBlock%nblocks
      b = b + DequationWeights(irow) * &
          lsyssc_scalarProduct(rtempVecBlock%RvectorBlock(irow),&
                               rcorrVector%RvectorBlock(irow))
    end do

    ! Return the alpha.
    if (b .ne. 0.0_DP) then
      dalpha = a/b
    else
      dalpha = 1.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cgcor_calcCorrDefMin (rmatrix,rvector,rrhs,rcorrVector,&
                                   rtempVecBlock,p_RfilterChain,dalpha)

!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
!</description>

!<input>
  ! The block matrix of the system Ax=b which is solved by multigrid
  type(t_matrixBlock), intent(in) :: rmatrix

  ! The current (uncorrected) solution vector
  type(t_vectorBlock), intent(in) :: rvector

  ! The current RHS vector
  type(t_vectorBlock), intent(in) :: rrhs

  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=<tex>$P^{-1}(b-Ax)$</tex> and <tex>$P^{-1}=$</tex> multigrid on the
  ! coarse level.
  type(t_vectorBlock), intent(in) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  type(t_filterChain), dimension(:), pointer :: p_RfilterChain
!</input>

!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  type(t_vectorBlock), intent(inout) :: rtempVecBlock
!</inputoutput>

!<output>
  ! The optimal correction parameter <tex>$\alpha$</tex> for the coarse grid correction.
  real(DP), intent(out) :: dalpha
!</output>

!</subroutine>

    ! local variables

    type(t_vectorBlock) :: rtempBlock2
    real(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek`s book):
    !
    !             ( f_k - A_k x_k  ,  A_k corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  A_k corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction

    call lsysbl_copyVector(rrhs,rtempVecBlock)
    call lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    if (associated(p_RfilterChain)) then
      call filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    end if

    ! We need a second temp vector which is to be released afterwards.
    ! (Perhaps in a later implementation, we could describe rtempVecBlock
    ! to have double the size?)
    call lsysbl_createVecBlockIndirect (rtempVecBlock,rtempBlock2,.false.)

    call lsysbl_blockMatVec(rmatrix, rcorrVector, rtempBlock2, 1.0_DP,0.0_DP)
    ! Apply the filter
    if (associated(p_RfilterChain)) then
      call filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    end if

    a = lsysbl_scalarProduct(rtempVecBlock,rtempBlock2)

    ! Calculate the demoninator of the fraction

    b = lsysbl_scalarProduct(rtempBlock2,rtempBlock2)

    ! Return the alpha.
    if (b .ne. 0.0_DP) then
      dalpha = a/b
    else
      dalpha = 1.0_DP
    end if

    ! Release the 2nd temp vector
    call lsysbl_releaseVector (rtempBlock2)

  end subroutine

end module

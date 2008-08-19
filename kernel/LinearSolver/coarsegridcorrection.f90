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
!# 3.) cgcor_calcOptimalCorrection  - Calculate the optimal alpha
!#                                    for the correction.
!#
!#  FAQ - Frequently asked Questions
!# ----------------------------------
!# 1.) What is this p_DequationWeights-parameter in the t_coarseGridCorrection
!#     structure for? I cannot understand?
!#
!#   This parameter allows to multiply an equation by a given factor
!#   while calculating optimal correction factors. This is used when
!#   calculating residuals. Imagine you want to calculate the residual
!#   of a 2D Navier-Stokes system:
!#
!#    (r1)     (r1)   (A11  A12  B1)  (x1)
!#    (r2)  =  (r2) - (A21  A22  B2)  (x2)
!#    (rp)     (rp)   (B1^T B2^T   )  (xp)
!#
!#   Now imagine that you set the array p_DequationWeights to (1,1,-1); then,
!#   residuals are calculated of the following system:
!# 
!#    (r1)     ( f1)   ( A11    A12  B1)  (x1)
!#    (r2)  =  ( f2) - ( A21    A22  B2)  (x2)
!#    (rp)     (-fp)   (-B1^T  -B2^T   )  (xp)
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

MODULE coarsegridcorrection

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE filtersupport
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Method identifier for coarse grid correction">

  ! Standard constant damping. A correction factor of 1.0 is used for 
  ! the coarse grid correction if dalphaMin < 1.0 < dalphaMax.
  ! Otherwise, the constant dalphaMin is chosen.
  ! Remark: This is the optimal setting for a scalar equation with conformal
  !  elements; however, using nonconformal elements or nonscalar equations
  !  might make it necessary to switch to another method.
  INTEGER, PARAMETER :: CGCOR_STANDARD       = 0
  
  ! Damping by energy minimisation.
  ! Remark: This is the optimal setting for a scalar Laplace equation
  !  with nonconformal Rannacher-Turek element Ex30/Ex31.
  INTEGER, PARAMETER :: CGCOR_SCALARENERGYMIN = 1
  
  ! Damping by defect minimisation.
  ! Remark: Easy/cheap to calculate, but no functional analytic background.
  INTEGER, PARAMETER :: CGCOR_SCALARDEFMIN    = 2
  
!</constantblock>

!</constants>


!<types>

!<typeblock>

  ! Coarse grid correction structure; defines the behaviour of 
  ! how to calculate the optimal coarse grid correction.
  TYPE t_coarseGridCorrection
    
    ! Method identifier. This is one of the CGCOR_xxxx constants and specifies
    ! the algorithm that is used to calculate the optimal coarse grid
    ! correction.
    INTEGER :: ccorrectionType = CGCOR_STANDARD
    
    ! Minimum damping parameter
    REAL(DP) :: dalphaMin = -10.0_DP
    
    ! Maximum damping parameter
    REAL(DP) :: dalphaMax = 10.0_DP
    
    ! A list of weights for the different equations.
    ! If not associated, a standard value of 1.0 is assumed for every equation.
    ! If associated, the calculation routine for the optimal
    ! coarse grid correction will multiply residuals of the corresponding
    ! equation with this factor.
    ! Example: 2D Navier-Stokes: (1,1,-1)
    ! -> The pressure equation is multiplied by -1 before taking the norm
    !    of the residual.
    ! Can be used to symmetrise equations.
    REAL(DP), DIMENSION(:), POINTER :: p_DequationWeights => NULL()
    
  END TYPE
  
!</typeblock>

!</types>

  ! ***************************************************************************

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cgcor_init (rcoarseGridCorrection,nequations)
                                          
!<description>
  ! Initialises a coarse grid correction structure.
!</description>

!<input>
  ! OPTIONAL: Number of equations.
  ! If specified, the p_DequationWeights array in the coarse grid correction
  ! structure is allocated and initially filled with 1.0-values.
  INTEGER, INTENT(IN), OPTIONAL :: nequations
!</input>

!<output>
  ! A coarse grid correction structure to be initialised
  TYPE(t_coarseGridCorrection), INTENT(OUT) :: rcoarseGridCorrection
!</output>

    ! Initialise by default initialisation
    !
    ! Exception: The p_DequationWeights array.
    IF (PRESENT(nequations)) THEN
      ALLOCATE(rcoarseGridCorrection%p_DequationWeights(nequations))
      rcoarseGridCorrection%p_DequationWeights = 1.0_DP
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cgcor_release (rcoarseGridCorrection)
                                          
!<description>
  ! Releases a coarse grid correction structure.
!</description>

!<inputoutput>
  ! A coarse grid correction structure to be released
  TYPE(t_coarseGridCorrection), INTENT(INOUT) :: rcoarseGridCorrection
!</inputoutput>

    ! Release the pointer to the equation weights if associated
    IF (ASSOCIATED(rcoarseGridCorrection%p_DequationWeights)) THEN
      DEALLOCATE(rcoarseGridCorrection%p_DequationWeights)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cgcor_calcOptimalCorrection (rcoarseGridCorrection,&
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
  TYPE(t_coarseGridCorrection), INTENT(IN) :: rcoarseGridCorrection
  
  ! The block matrix of the system Ax=b which is solved by multigrid
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
  
  ! The current (uncorrected) solution vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
  
  ! The current RHS vector
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
  
  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=$P^{-1}(b-Ax)$ and $P^{-1}=$multigrid on the coarse level.
  TYPE(t_vectorBlock), INTENT(IN) :: rcorrVector
  
  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!</input>
  
!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! The optimal correction parameter $\alpha$ for the coarse grid correction.
  REAL(DP), INTENT(OUT) :: dalpha
!</output>

!</subroutine>

  IF (rcoarseGridCorrection%dalphaMax .LE. rcoarseGridCorrection%dalphaMin) THEN
    dalpha = 1.0_DP
  ELSE
    ! Which method to use?
    
    SELECT CASE (rcoarseGridCorrection%ccorrectionType)
    CASE (CGCOR_SCALARENERGYMIN)
      IF (.NOT. ASSOCIATED(rcoarseGridCorrection%p_DequationWeights)) THEN
        ! Standard energy minimisation
        CALL cgcor_calcCorrEnergyMin (rmatrix,rvector,rrhs,rcorrVector,&
                                      rtempVector,p_RfilterChain,dalpha)
      ELSE
        ! Special energy minimisation with different weights for the
        ! different equations.
        CALL cgcor_calcCorrEnergyMinWeighted (rmatrix,rvector,rrhs,rcorrVector,&
            rtempVector,p_RfilterChain,dalpha,rcoarseGridCorrection%p_DequationWeights)
      END IF
    CASE (CGCOR_SCALARDEFMIN)
      CALL cgcor_calcCorrDefMin (rmatrix,rvector,rrhs,rcorrVector,&
                                rtempVector,p_RfilterChain,dalpha)
    CASE DEFAULT !(=CGCOR_STANDARD)
      dalpha = 1.0_DP   ! Standard setting
    END SELECT
  END IF

  ! Make sure it's in the interval given by the dalphaMin/dalphaMax
  dalpha = MAX(MIN(dalpha,rcoarseGridCorrection%dalphaMax), &
                          rcoarseGridCorrection%dalphaMin)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cgcor_calcCorrEnergyMin (rmatrix,rvector,rrhs,rcorrVector,&
                                      rtempVecBlock,p_RfilterChain,dalpha)
                                          
!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
!</description>

!<input>
  ! The block matrix of the system Ax=b which is solved by multigrid
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
  
  ! The current (uncorrected) solution vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
  
  ! The current RHS vector
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
  
  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=$P^{-1}(b-Ax)$ and $P^{-1}=$multigrid on the coarse level.
  TYPE(t_vectorBlock), INTENT(IN) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!</input>
  
!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecBlock
!</inputoutput>

!<output>
  ! The optimal correction parameter $\alpha$ for the coarse grid correction.
  REAL(DP), INTENT(OUT) :: dalpha
!</output>

!</subroutine>

    ! local variables
    
    REAL(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek's book):
    !
    !             ( f_k - A_k x_k  ,  corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction
      
    CALL lsysbl_copyVector(rrhs,rtempVecBlock)
    CALL lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
    a = lsysbl_scalarProduct(rtempVecBlock,rcorrVector)
    
    ! Calculate the demoninator of the fraction
    CALL lsysbl_blockMatVec(rmatrix, rcorrVector, rtempVecBlock, 1.0_DP,0.0_DP)
    ! Apply the filter
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
    b = lsysbl_scalarProduct(rtempVecBlock,rcorrVector)
    
    ! Return the alpha.
    IF (b .NE. 0.0_DP) THEN
      dalpha = a/b
    ELSE
      dalpha = 1.0_DP
    END IF
    
  END SUBROUTINE  

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cgcor_calcCorrEnergyMinWeighted (rmatrix,rvector,rrhs,rcorrVector,&
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
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
  
  ! The current (uncorrected) solution vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
  
  ! The current RHS vector
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
  
  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=$P^{-1}(b-Ax)$ and $P^{-1}=$multigrid on the coarse level.
  TYPE(t_vectorBlock), INTENT(IN) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
  ! List of weights for every equation (=row of the matrix)
  REAL(DP), DIMENSION(:), INTENT(IN) :: DequationWeights
!</input>
  
!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecBlock
!</inputoutput>

!<output>
  ! The optimal correction parameter $\alpha$ for the coarse grid correction.
  REAL(DP), INTENT(OUT) :: dalpha
!</output>

!</subroutine>

    ! local variables
    
    REAL(DP) :: a,b
    INTEGER :: irow

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek's book):
    !
    !             ( f_k - A_k x_k  ,  corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction
      
    CALL lsysbl_copyVector(rrhs,rtempVecBlock)
    CALL lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
    ! Calculate a weighted scalar product using the equation weights.
    a = 0.0_DP
    DO irow = 1,rtempVecBlock%nblocks
      a = a + DequationWeights(irow) * &
          lsyssc_scalarProduct(rtempVecBlock%RvectorBlock(irow),&
                               rcorrVector%RvectorBlock(irow))
    END DO
    
    ! Calculate the demoninator of the fraction
    CALL lsysbl_blockMatVec(rmatrix, rcorrVector, rtempVecBlock, 1.0_DP,0.0_DP)
    ! Apply the filter
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
    b = 0.0_DP
    DO irow = 1,rtempVecBlock%nblocks
      b = b + DequationWeights(irow) * &
          lsyssc_scalarProduct(rtempVecBlock%RvectorBlock(irow),&
                               rcorrVector%RvectorBlock(irow))
    END DO
    
    ! Return the alpha.
    IF (b .NE. 0.0_DP) THEN
      dalpha = a/b
    ELSE
      dalpha = 1.0_DP
    END IF
    
  END SUBROUTINE  

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cgcor_calcCorrDefMin (rmatrix,rvector,rrhs,rcorrVector,&
                                   rtempVecBlock,p_RfilterChain,dalpha)
                                          
!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
!</description>

!<input>
  ! The block matrix of the system Ax=b which is solved by multigrid
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
  
  ! The current (uncorrected) solution vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
  
  ! The current RHS vector
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
  
  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=$P^{-1}(b-Ax)$ and $P^{-1}=$multigrid on the coarse level.
  TYPE(t_vectorBlock), INTENT(IN) :: rcorrVector

  ! Either NULL() or a pointer to a filter chain which must be applied
  ! to every defect vector (b-Ax).
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
!</input>
  
!<inputoutput>
  ! A block temporary vector. Must have the same size and structure as
  ! the RHS and the solution vector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVecBlock
!</inputoutput>

!<output>
  ! The optimal correction parameter $\alpha$ for the coarse grid correction.
  REAL(DP), INTENT(OUT) :: dalpha
!</output>

!</subroutine>

    ! local variables
    
    TYPE(t_vectorBlock) :: rtempBlock2
    REAL(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek's book):
    !
    !             ( f_k - A_k x_k  ,  A_k corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  A_k corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !    
    ! Calculate nominator of the fraction
      
    CALL lsysbl_copyVector(rrhs,rtempVecBlock)
    CALL lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
    ! We need a second temp vector which is to be released afterwards.
    ! (Perhaps in a later implementation, we could describe rtempVecBlock
    ! to have double the size?)
    CALL lsysbl_createVecBlockIndirect (rtempVecBlock,rtempBlock2,.FALSE.)
    
    CALL lsysbl_blockMatVec(rmatrix, rcorrVector, rtempBlock2, 1.0_DP,0.0_DP)
    ! Apply the filter
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
    a = lsysbl_scalarProduct(rtempVecBlock,rtempBlock2)
    
    ! Calculate the demoninator of the fraction
    
    b = lsysbl_scalarProduct(rtempBlock2,rtempBlock2)
    
    ! Return the alpha.
    IF (b .NE. 0.0_DP) THEN
      dalpha = a/b
    ELSE
      dalpha = 1.0_DP
    END IF
    
    ! Release the 2nd temp vector
    CALL lsysbl_releaseVector (rtempBlock2)
    
  END SUBROUTINE  

END MODULE

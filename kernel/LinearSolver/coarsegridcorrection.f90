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
!# 1.) cgcor_calcOptimalCorrection  - Calculate the optimal alpha
!#                                    for the correction.
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

!<constantblocks description="Method identifier for coarse grid correction">

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
    
  END TYPE
  
!</typeblock>

!</types>

  ! ***************************************************************************

CONTAINS

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
  ! A scalar temporary vector. Must be at least as large as rvector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtempVector
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
      CALL cgcor_calcCorrEnergyMin (rmatrix,rvector,rrhs,rcorrVector,&
                                    rtempVector,p_RfilterChain,dalpha)
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
                                      rtempVector,p_RfilterChain,dalpha)
                                          
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
  ! A scalar temporary vector. Must be at least as large as rvector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! The optimal correction parameter $\alpha$ for the coarse grid correction.
  REAL(DP), INTENT(OUT) :: dalpha
!</output>

!</subroutine>

    ! local variables
    
    TYPE(t_vectorBlock) :: rtempVecBlock
    REAL(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek's book):
    !
    !             ( f_k - A_k x_k  ,  corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  corr_k )
    !
    ! For this purpose, we need a temporary vector.
    ! Ok, the above vector is scalar, but long enough. For simplicity,
    ! we enforce the structure of rvector (=x_k) to the scalar vector and
    ! use its memory as block vector:
    
    CALL lsysbl_createVecFromScalar (rtempVector,rtempVecBlock)
    CALL lsysbl_enforceStructure (rvector,rtempVecBlock)
    
    ! Now, rtempVecBlock uses the memory of rtempVector but has the
    ! structure of rvector. Be careful not to release it :-)
    ! (as the caller is the owner of the vector)
    !
    ! Calculate nominator of the fraction
      
    CALL lsysbl_vectorCopy(rrhs,rtempVecBlock)
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

  SUBROUTINE cgcor_calcCorrDefMin (rmatrix,rvector,rrhs,rcorrVector,&
                                   rtempVector,p_RfilterChain,dalpha)
                                          
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
  ! A scalar temporary vector. Must be at least as large as rvector.
  ! The content is undefined on entry of this routine and will be
  ! undefined when the routine finishes.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! The optimal correction parameter $\alpha$ for the coarse grid correction.
  REAL(DP), INTENT(OUT) :: dalpha
!</output>

!</subroutine>

    ! local variables
    
    TYPE(t_vectorBlock) :: rtempVecBlock, rtempBlock2
    REAL(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek's book):
    !
    !             ( f_k - A_k x_k  ,  A_k corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  A_k corr_k )
    !
    ! For this purpose, we need a temporary vector.
    ! Ok, the above vector is scalar, but long enough. For simplicity,
    ! we enforce the structure of rvector (=x_k) to the scalar vector and
    ! use its memory as block vector:
    
    CALL lsysbl_createVecFromScalar (rtempVector,rtempVecBlock)
    CALL lsysbl_enforceStructure (rvector,rtempVecBlock)
    
    ! We need a second temp vector which is to be released afterwards.
    ! (Perhaps in a later implementation, we could describe rtempVecBlock
    ! to have double the size?)
    CALL lsysbl_createVecBlockIndirect (rtempVecBlock,rtempBlock2,.FALSE.)
    
    ! Now, rtempVecBlock uses the memory of rtempVector but has the
    ! structure of rvector. Be careful not to release it :-)
    ! (as the caller is the owner of the vector)
    !
    ! Calculate nominator of the fraction
      
    CALL lsysbl_vectorCopy(rrhs,rtempVecBlock)
    CALL lsysbl_blockMatVec(rmatrix, rvector, rtempVecBlock, -1.0_DP,1.0_DP)
    ! This is a defect vector - apply the filter chain.
    IF (ASSOCIATED(p_RfilterChain)) THEN
      CALL filter_applyFilterChainVec (rtempVecBlock,p_RfilterChain)
    END IF
    
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
    
    ! Release the 2nd vector
    CALL lsysbl_releaseVector(rtempBlock2)
    
  END SUBROUTINE  

END MODULE

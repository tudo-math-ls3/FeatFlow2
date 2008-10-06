!##############################################################################
!# ****************************************************************************
!# <name> prolrest2d_test2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module discreticises the "PDE"
!#
!#                                 u = f
!#
!# on two levels and solves the systems on both levels.
!# Afterwards, the fine mesh solution is compared with the prolongated
!# coarse mesh solution and the coarse mesh RHS is compared with the
!# restricted fine mesh RHS.
!#
!# The grid transfer is performed using true L2-projection based on the
!# 2-Level-Mass matrix assembled by the multileveloperators module.
!# </purpose>
!##############################################################################

MODULE prolrest2d_test2

  USE fsystem
  USE genoutput
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE triangulation
  USE spatialdiscretisation
  USE genoutput
  USE multileveloperators

  USE prolrest_aux
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE prolrest2d_2
  
!<description>
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriaC, rtriaF

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscrC, rdiscrF
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrixC, rmatrixF
    TYPE(t_vectorBlock) :: rsolC,rsolF,rrhsC,rrhsF,rtempC,rtempF
    
    TYPE(t_matrixScalar) :: rmatrix2Lvl,rmatrix2LvlT
    
    REAL(DP), DIMENSION(:), POINTER :: p_DsolF,p_DsolC,p_DrhsC,p_DrhsF,&
                                       p_DtmpF,p_DtmpC
    
    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverC,p_rsolverF

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    
    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMIN,i
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    

    NLMIN = 1
    
    CALL output_lbrk()
    CALL output_separator(OU_SEP_STAR)
    CALL output_line('Test 2: L2-Projection Operators VS Other Discretisation')
    CALL output_separator(OU_SEP_STAR)

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic coarse triangulation.
    CALL tria_readTriFile2D (rtriaC, './pre/QUAD.tri', rboundary)
     
    ! Create the coarse mesh
    CALL tria_quickRefine2LevelOrdering (NLMIN-1,rtriaC,rboundary)
    CALL tria_initStandardMeshFromRaw (rtriaC,rboundary)
    
    ! Create the fine mesh
    CALL tria_refine2LevelOrdering(rtriaC,rtriaF,rboundary)
    CALL tria_initStandardMeshFromRaw (rtriaF,rboundary)

    ! Create 2 discretisations for both meshes
    CALL spdiscr_initBlockDiscr2D (rdiscrC,1,rtriaC, rboundary)
    CALL spdiscr_initBlockDiscr2D (rdiscrF,1,rtriaF, rboundary)

    CALL spdiscr_initDiscr_simple (rdiscrC%RspatialDiscr(1), &
                                   EL_E037,CUB_G5X5,rtriaC, rboundary)
    CALL spdiscr_initDiscr_simple (rdiscrF%RspatialDiscr(1), &
                                   EL_E037,CUB_G5X5,rtriaF, rboundary)

    ! Allocate two block matrices
    CALL lsysbl_createMatBlockByDiscr (rdiscrC,rmatrixC)
    CALL lsysbl_createMatBlockByDiscr (rdiscrF,rmatrixF)
    
    ! And create the matrix structures
    CALL bilf_createMatrixStructure (rdiscrC%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixC%RmatrixBlock(1,1))
    CALL bilf_createMatrixStructure (rdiscrF%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixF%RmatrixBlock(1,1))
    
    ! Update the block matrices
    CALL lsysbl_updateMatStrucInfo (rmatrixC)
    CALL lsysbl_updateMatStrucInfo (rmatrixF)
    
    ! Set up the bilinearform for a mass matrix
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = 1.0 

    ! Build the two mass matrices
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrixC%RmatrixBlock(1,1))
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrixF%RmatrixBlock(1,1))
    
    ! Create the 2-Level matrix structure
    CALL mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
            rdiscrF%RspatialDiscr(1),LSYSSC_MATRIX9,rmatrix2Lvl)
    
    ! Build the 2-Level-Mass matrix
    CALL mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
                 rdiscrF%RspatialDiscr(1),.TRUE.,rmatrix2Lvl)
                 
    ! Transpose the 2-Level-Mass matrix
    CALL lsyssc_transposeMatrix(rmatrix2Lvl,rmatrix2LvlT,LSYSSC_TR_VIRTUAL)
    
    ! Allocate two block vectors for the RHS
    CALL lsysbl_createVecBlockIndMat (rmatrixC,rrhsC, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixF,rrhsF, .FALSE.)
    
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! Build the two RHS vectors
    CALL linf_buildVectorScalar (rdiscrC%RspatialDiscr(1),&
                 rlinform,.TRUE.,rrhsC%rvectorBlock(1),procRHS_Q2_2D)
    CALL linf_buildVectorScalar (rdiscrF%RspatialDiscr(1),&
                 rlinform,.TRUE.,rrhsF%rvectorBlock(1),procRHS_Q2_2D)
    
    ! Allocate solution and temporary vectors
    CALL lsysbl_createVecBlockIndMat (rmatrixC,rsolC, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixF,rsolF, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixC,rtempC, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixF,rtempF, .FALSE.)
    
    ! Clear solution vectors
    CALL lsysbl_clearVector(rsolC)
    CALL lsysbl_clearVector(rsolF)
    
    ! Create two UMFPACK solvers
    CALL linsol_initUMFPACK4(p_rsolverC)
    CALL linsol_initUMFPACK4(p_rsolverF)
    p_rsolverC%ioutputLevel = 2
    p_rsolverF%ioutputLevel = 2
    
    ! Initialise both solvers
    Rmatrices = (/rmatrixC/)
    CALL linsol_setMatrices(p_rsolverC,Rmatrices)
    CALL linsol_initStructure (p_rsolverC, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverC, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    Rmatrices = (/rmatrixF/)
    CALL linsol_setMatrices(p_rsolverF,Rmatrices)
    CALL linsol_initStructure (p_rsolverF, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverF, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Solve both systems
    CALL linsol_solveAdaptively (p_rsolverC,rsolC,rrhsC,rtempC)
    CALL linsol_solveAdaptively (p_rsolverF,rsolF,rrhsF,rtempF)
    
    ! Multiply 2-Level-Mass matrix with coarse mesh solution vector
    ! to get a rhs vector for the fine mesh.
    CALL lsyssc_scalarMatVec(rmatrix2Lvl, rsolC%RvectorBlock(1), &
                             rtempF%RvectorBlock(1), 1.0_DP, 0.0_DP)
    
    ! Multiply transposed 2-Level-Mass matrix with fine mesh solution
    ! vector to get a rhs vector for the coarse mesh.
    CALL lsyssc_scalarMatVec(rmatrix2LvlT, rsolF%RvectorBlock(1), &
                             rtempC%RvectorBlock(1), 1.0_DP, 0.0_DP)
    
    ! Get all the arrays
    CALL lsysbl_getbase_double(rsolF, p_DsolF)
    CALL lsysbl_getbase_double(rsolC, p_DsolC)
    CALL lsysbl_getbase_double(rrhsF, p_DrhsF)
    CALL lsysbl_getbase_double(rrhsC, p_DrhsC)
    CALL lsysbl_getbase_double(rtempF, p_DtmpF)
    CALL lsysbl_getbase_double(rtempC, p_DtmpC)
    
    ! -------------------------------------------------------------------------
    ! 2-Level-Mass Test
    ! -------------------------------------------------------------------------
    CALL output_separator(OU_SEP_MINUS)
    CALL output_line('2-Level-Mass Test')
    CALL output_separator(OU_SEP_MINUS)

    ! Calculate error of prolongation
    CALL lsysbl_vectorLinearComb(rrhsF,rtempF,1.0_DP,-1.0_DP,rsolF)
    CALL vec_filterByEps(p_DsolF)
    
    ! Print out all DOFs
    !                0         1         2         3         4         5
    !                 123456789-123456789-123456789-123456789-123456789-123456789
    !                  1234----12345678901234567890----12345678901234567890----
    CALL output_line(' DOF     fine rhs                prol. coarse solution   error')
    DO i=1, rsolF%NEQ
      CALL output_line(TRIM(sys_si(i,4)) // '    ' // &
        TRIM(sys_sdEP(p_DrhsF(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_DtmpF(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_DsolF(i),20,13)))
    END DO    
    
    ! -------------------------------------------------------------------------
    ! Transposed 2-Level-Mass Test
    ! -------------------------------------------------------------------------
    CALL output_separator(OU_SEP_MINUS)
    CALL output_line('Transposed 2-Level-Mass Test')
    CALL output_separator(OU_SEP_MINUS)

    ! Calculate error of restriction
    CALL lsysbl_vectorLinearComb(rrhsC,rtempC,1.0_DP,-1.0_DP,rsolC)
    CALL vec_filterByEps(p_DsolC)

    ! Print out all DOFs
    !                0         1         2         3         4         5
    !                 123456789-123456789-123456789-123456789-123456789-123456789
    !                  1234----12345678901234567890----12345678901234567890----
    CALL output_line(' DOF     coarse rhs              rest. fine solution     error')
    DO i=1, rrhsC%NEQ
      CALL output_line(TRIM(sys_si(i,4)) // '    ' // &
        TRIM(sys_sdEP(p_DrhsC(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_DtmpC(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_DsolC(i),20,13)))
    END DO    
    
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverF)
    CALL linsol_doneStructure (p_rsolverF)
    CALL linsol_doneData (p_rsolverC)
    CALL linsol_doneStructure (p_rsolverC)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverF)
    CALL linsol_releaseSolver (p_rsolverC)

    ! Release the 2-Level matrices
    CALL lsyssc_releaseMatrix (rmatrix2LvlT)
    CALL lsyssc_releaseMatrix (rmatrix2Lvl)

    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rtempF)
    CALL lsysbl_releaseVector (rtempC)
    CALL lsysbl_releaseVector (rsolF)
    CALL lsysbl_releaseVector (rsolC)
    CALL lsysbl_releaseVector (rrhsF)
    CALL lsysbl_releaseVector (rrhsC)
    CALL lsysbl_releaseMatrix (rmatrixF)
    CALL lsysbl_releaseMatrix (rmatrixC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscrF)
    CALL spdiscr_releaseBlockDiscr(rdiscrC)
    
    ! Release the triangulation. 
    CALL tria_done (rtriaF)
    CALL tria_done (rtriaC)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
  END SUBROUTINE

END MODULE

!##############################################################################
!# ****************************************************************************
!# <name> prolrest2d_test3 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module discreticises the "PDE"
!#
!#                                 u = f
!#
!# on two levels and solves the systems on both levels.
!# Afterwards, the prolongation / restriction results of the hard-coded
!# inter-level projection structure defined in the multilevelprojection
!# module are compared with the true L2-projection based on the
!# 2-Level-Mass matrix assembled by the multileveloperators module.
!# </purpose>
!##############################################################################

MODULE prolrest2d_test3

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
  USE ucd
  USE pprocerror
  USE genoutput
  USE multilevelprojection
  USE multileveloperators
    
  USE prolrest_aux
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE prolrest2d_3
  
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
    TYPE(t_vectorBlock) :: rsolC,rsolF,rrhsC,rrhsF,rtempC,rtempF,&
                           rvecProl,rvecRest,rvecL2Prol,rvecL2Rest
    
    REAL(DP), DIMENSION(:), POINTER :: p_DtmpF,p_DtmpC,p_Dprol,p_Drest,&
                                       p_DL2prol,p_DL2rest

    TYPE(t_matrixScalar) :: rmatrix2Lvl
    
    ! A temporary vector for the projection
    TYPE(t_vectorScalar) :: rprjTmpVec

    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rproj
    TYPE(t_interlevelProjectionBlock) :: rL2proj

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverC,p_rsolverF

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    
    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMIN,i
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    

    NLMIN = 2
    
    CALL output_lbrk()
    CALL output_separator(OU_SEP_STAR)
    CALL output_line('Test 3: L2-Projection Operators VS Hard-Coded Operators')
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
                                   EL_Q1,CUB_G5X5,rtriaC, rboundary)
    CALL spdiscr_initDiscr_simple (rdiscrF%RspatialDiscr(1), &
                                   EL_Q1,CUB_G5X5,rtriaF, rboundary)

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

    ! Create the 2-Level matrix structure
    CALL mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
            rdiscrF%RspatialDiscr(1),LSYSSC_MATRIX9,rmatrix2Lvl)
    
    ! Build the 2-Level-Mass matrix
    CALL mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
                 rdiscrF%RspatialDiscr(1),.TRUE.,rmatrix2Lvl)

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
    
    ! Create the vectors for prolongation, restriction and interpolation
    CALL lsysbl_createVecBlockIndMat (rmatrixC,rvecRest,.FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixF,rvecProl,.FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixC,rvecL2Rest,.FALSE.)
    CALL lsysbl_createVecBlockIndMat (rmatrixF,rvecL2Prol,.FALSE.)

    ! Initialise a two standard interlevel projection structure
    CALL mlprj_initProjectionMat (rproj,rmatrixF)
    CALL mlprj_initProjectionMat (rL2proj,rmatrixF)
    
    ! Set matrices for L2-projection
    CALL mlprj_initL2Proj (rL2proj%RscalarProjection(1,1), &
                           rmatrix2Lvl, rmatrixF%RmatrixBlock(1,1))
    !rL2proj%RscalarProjection(1,1)%depsAbsL2 = SYS_EPSREAL
    rL2proj%RscalarProjection(1,1)%imaxL2Iterations = 100
    
    ! Create a dummy scalar vector
    CALL lsyssc_createVector(rprjTmpVec, 0, .FALSE.)
    
    ! Prolongate coarse mesh solution vector
    CALL mlprj_performProlongation (rproj,rsolC,rvecProl,rprjTmpVec)
    CALL mlprj_performProlongation (rL2proj,rsolC,rvecL2Prol,rprjTmpVec)

    ! Restrict fine mesh right-hand-side vector
    CALL mlprj_performRestriction(rproj,rvecRest,rrhsF,rprjTmpVec)
    CALL mlprj_performRestriction(rL2proj,rvecL2Rest,rrhsF,rprjTmpVec)
    
    ! Get all the arrays
    CALL lsysbl_getbase_double(rvecProl, p_Dprol)
    CALL lsysbl_getbase_double(rvecRest, p_Drest)
    CALL lsysbl_getbase_double(rvecL2Prol, p_DL2prol)
    CALL lsysbl_getbase_double(rvecL2Rest, p_DL2rest)
    CALL lsysbl_getbase_double(rtempF, p_DtmpF)
    CALL lsysbl_getbase_double(rtempC, p_DtmpC)
    
    ! -------------------------------------------------------------------------
    ! Prolongation Test
    ! -------------------------------------------------------------------------
    CALL output_separator(OU_SEP_MINUS)
    CALL output_line('Prolongation Test')
    CALL output_separator(OU_SEP_MINUS)
    
    ! Calculate error of prolongation
    CALL lsysbl_vectorLinearComb(rvecL2Prol,rvecProl,1.0_DP,-1.0_DP,rtempF)
    CALL vec_filterByEps(p_DtmpF)

    ! Print out all DOFs
    !                0         1         2         3         4         5
    !                 123456789-123456789-123456789-123456789-123456789-123456789
    !                  1234----12345678901234567890----12345678901234567890----
    CALL output_line(' DOF     L2-projection           std. prolongation       error')
    DO i=1, rsolF%NEQ
      CALL output_line(TRIM(sys_si(i,4)) // '    ' // &
        TRIM(sys_sdEP(p_DL2prol(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_Dprol(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_DtmpF(i),20,13)))
    END DO    
    
    ! -------------------------------------------------------------------------
    ! Restriction Test
    ! -------------------------------------------------------------------------
    CALL output_separator(OU_SEP_MINUS)
    CALL output_line('Restriction Test')
    CALL output_separator(OU_SEP_MINUS)

    ! Calculate error of restriction
    CALL lsysbl_vectorLinearComb(rvecL2Rest,rvecRest,1.0_DP,-1.0_DP,rtempC)
    CALL vec_filterByEps(p_DtmpC)

    ! Print out all DOFs
    !                0         1         2         3         4         5
    !                 123456789-123456789-123456789-123456789-123456789-123456789
    !                  1234----12345678901234567890----12345678901234567890----
    CALL output_line(' DOF     L2-projection           std. restriction        error')
    DO i=1, rrhsC%NEQ
      CALL output_line(TRIM(sys_si(i,4)) // '    ' // &
        TRIM(sys_sdEP(p_DL2rest(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_Drest(i),20,13)) // '    ' // &
        TRIM(sys_sdEP(p_DtmpC(i),20,13)))
    END DO    
   
    ! Release the projection
    CALL mlprj_doneProjection (rL2proj)
    CALL mlprj_doneProjection (rproj)
    
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverF)
    CALL linsol_doneStructure (p_rsolverF)
    CALL linsol_doneData (p_rsolverC)
    CALL linsol_doneStructure (p_rsolverC)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverF)
    CALL linsol_releaseSolver (p_rsolverC)
    
    ! Release the 2-Level-Mass matrix
    CALL lsyssc_releaseMatrix (rmatrix2Lvl)

    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rvecL2Prol)
    CALL lsysbl_releaseVector (rvecL2Rest)
    CALL lsysbl_releaseVector (rvecProl)
    CALL lsysbl_releaseVector (rvecRest)
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

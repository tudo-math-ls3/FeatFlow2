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

module prolrest2d_test2

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use triangulation
  use spatialdiscretisation
  use genoutput
  use multileveloperators

  use prolrest_aux
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine prolrest2d_2
  
!<description>
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriaC, rtriaF

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscrC, rdiscrF
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixC, rmatrixF
    type(t_vectorBlock) :: rsolC,rsolF,rrhsC,rrhsF,rtempC,rtempF
    
    type(t_matrixScalar) :: rmatrix2Lvl,rmatrix2LvlT
    
    real(DP), dimension(:), pointer :: p_DsolF,p_DsolC,p_DrhsC,p_DrhsF,&
                                       p_DtmpF,p_DtmpC
    
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverC,p_rsolverF

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMIN,i
    
    ! Error indicator during initialisation of the solver
    integer :: ierror

    NLMIN = 1
    
    call output_lbrk()
    call output_separator(OU_SEP_STAR)
    call output_line('Test 2: L2-Projection Operators VS Other Discretisation')
    call output_separator(OU_SEP_STAR)

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic coarse triangulation.
    call tria_readTriFile2D (rtriaC, './pre/QUAD.tri', rboundary)
     
    ! Create the coarse mesh
    call tria_quickRefine2LevelOrdering (NLMIN-1,rtriaC,rboundary)
    call tria_initStandardMeshFromRaw (rtriaC,rboundary)
    
    ! Create the fine mesh
    call tria_refine2LevelOrdering(rtriaC,rtriaF,rboundary)
    call tria_initStandardMeshFromRaw (rtriaF,rboundary)

    ! Create 2 discretisations for both meshes
    call spdiscr_initBlockDiscr (rdiscrC,1,rtriaC, rboundary)
    call spdiscr_initBlockDiscr (rdiscrF,1,rtriaF, rboundary)

    call spdiscr_initDiscr_simple (rdiscrC%RspatialDiscr(1), &
                                   EL_Q2TB,CUB_G5X5,rtriaC, rboundary)
    call spdiscr_initDiscr_simple (rdiscrF%RspatialDiscr(1), &
                                   EL_Q2TB,CUB_G5X5,rtriaF, rboundary)

    ! Allocate two block matrices
    call lsysbl_createMatBlockByDiscr (rdiscrC,rmatrixC)
    call lsysbl_createMatBlockByDiscr (rdiscrF,rmatrixF)
    
    ! And create the matrix structures
    call bilf_createMatrixStructure (rdiscrC%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixC%RmatrixBlock(1,1))
    call bilf_createMatrixStructure (rdiscrF%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixF%RmatrixBlock(1,1))
    
    ! Update the block matrices
    call lsysbl_updateMatStrucInfo (rmatrixC)
    call lsysbl_updateMatStrucInfo (rmatrixF)
    
    ! Set up the bilinearform for a mass matrix
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = 1.0

    ! Build the two mass matrices
    call bilf_buildMatrixScalar (rform,.true.,rmatrixC%RmatrixBlock(1,1))
    call bilf_buildMatrixScalar (rform,.true.,rmatrixF%RmatrixBlock(1,1))
    
    ! Create the 2-Level matrix structure
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
            rdiscrF%RspatialDiscr(1),LSYSSC_MATRIX9,rmatrix2Lvl)
    
    ! Build the 2-Level-Mass matrix
    call mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
                 rdiscrF%RspatialDiscr(1),.true.,rmatrix2Lvl)
                 
    ! Transpose the 2-Level-Mass matrix
    call lsyssc_transposeMatrix(rmatrix2Lvl,rmatrix2LvlT,LSYSSC_TR_VIRTUAL)
    
    ! Allocate two block vectors for the RHS
    call lsysbl_createVecBlockIndMat (rmatrixC,rrhsC, .false.)
    call lsysbl_createVecBlockIndMat (rmatrixF,rrhsF, .false.)
    
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! Build the two RHS vectors
    call linf_buildVectorScalar (rdiscrC%RspatialDiscr(1),&
                 rlinform,.true.,rrhsC%rvectorBlock(1),procRHS_Q2_2D)
    call linf_buildVectorScalar (rdiscrF%RspatialDiscr(1),&
                 rlinform,.true.,rrhsF%rvectorBlock(1),procRHS_Q2_2D)
    
    ! Allocate solution and temporary vectors
    call lsysbl_createVecBlockIndMat (rmatrixC,rsolC, .false.)
    call lsysbl_createVecBlockIndMat (rmatrixF,rsolF, .false.)
    call lsysbl_createVecBlockIndMat (rmatrixC,rtempC, .false.)
    call lsysbl_createVecBlockIndMat (rmatrixF,rtempF, .false.)
    
    ! Clear solution vectors
    call lsysbl_clearVector(rsolC)
    call lsysbl_clearVector(rsolF)
    
    ! Create two UMFPACK solvers
    call linsol_initUMFPACK4(p_rsolverC)
    call linsol_initUMFPACK4(p_rsolverF)
    p_rsolverC%ioutputLevel = 2
    p_rsolverF%ioutputLevel = 2
    
    ! Initialise both solvers
    Rmatrices = (/rmatrixC/)
    call linsol_setMatrices(p_rsolverC,Rmatrices)
    call linsol_initStructure (p_rsolverC, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverC, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    Rmatrices = (/rmatrixF/)
    call linsol_setMatrices(p_rsolverF,Rmatrices)
    call linsol_initStructure (p_rsolverF, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverF, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    
    ! Solve both systems
    call linsol_solveAdaptively (p_rsolverC,rsolC,rrhsC,rtempC)
    call linsol_solveAdaptively (p_rsolverF,rsolF,rrhsF,rtempF)
    
    ! Multiply 2-Level-Mass matrix with coarse mesh solution vector
    ! to get a rhs vector for the fine mesh.
    call lsyssc_scalarMatVec(rmatrix2Lvl, rsolC%RvectorBlock(1), &
                             rtempF%RvectorBlock(1), 1.0_DP, 0.0_DP)
    
    ! Multiply transposed 2-Level-Mass matrix with fine mesh solution
    ! vector to get a rhs vector for the coarse mesh.
    call lsyssc_scalarMatVec(rmatrix2LvlT, rsolF%RvectorBlock(1), &
                             rtempC%RvectorBlock(1), 1.0_DP, 0.0_DP)
    
    ! Get all the arrays
    call lsysbl_getbase_double(rsolF, p_DsolF)
    call lsysbl_getbase_double(rsolC, p_DsolC)
    call lsysbl_getbase_double(rrhsF, p_DrhsF)
    call lsysbl_getbase_double(rrhsC, p_DrhsC)
    call lsysbl_getbase_double(rtempF, p_DtmpF)
    call lsysbl_getbase_double(rtempC, p_DtmpC)
    
    ! -------------------------------------------------------------------------
    ! 2-Level-Mass Test
    ! -------------------------------------------------------------------------
    call output_separator(OU_SEP_MINUS)
    call output_line('2-Level-Mass Test')
    call output_separator(OU_SEP_MINUS)

    ! Calculate error of prolongation
    call lsysbl_vectorLinearComb(rrhsF,rtempF,1.0_DP,-1.0_DP,rsolF)
    call vec_filterByEps(p_DsolF)
    
    ! Print out all DOFs
    !                0         1         2         3         4         5
    !                 123456789-123456789-123456789-123456789-123456789-123456789
    !                  1234----12345678901234567890----12345678901234567890----
    call output_line(' DOF     fine rhs                prol. coarse solution   error')
    do i=1, rsolF%NEQ
      call output_line(trim(sys_si(i,4)) // '    ' // &
        trim(sys_sdEP(p_DrhsF(i),20,13)) // '    ' // &
        trim(sys_sdEP(p_DtmpF(i),20,13)) // '    ' // &
        trim(sys_sdEP(p_DsolF(i),20,13)))
    end do
    
    ! -------------------------------------------------------------------------
    ! Transposed 2-Level-Mass Test
    ! -------------------------------------------------------------------------
    call output_separator(OU_SEP_MINUS)
    call output_line('Transposed 2-Level-Mass Test')
    call output_separator(OU_SEP_MINUS)

    ! Calculate error of restriction
    call lsysbl_vectorLinearComb(rrhsC,rtempC,1.0_DP,-1.0_DP,rsolC)
    call vec_filterByEps(p_DsolC)

    ! Print out all DOFs
    !                0         1         2         3         4         5
    !                 123456789-123456789-123456789-123456789-123456789-123456789
    !                  1234----12345678901234567890----12345678901234567890----
    call output_line(' DOF     coarse rhs              rest. fine solution     error')
    do i=1, rrhsC%NEQ
      call output_line(trim(sys_si(i,4)) // '    ' // &
        trim(sys_sdEP(p_DrhsC(i),20,13)) // '    ' // &
        trim(sys_sdEP(p_DtmpC(i),20,13)) // '    ' // &
        trim(sys_sdEP(p_DsolC(i),20,13)))
    end do
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverF)
    call linsol_doneStructure (p_rsolverF)
    call linsol_doneData (p_rsolverC)
    call linsol_doneStructure (p_rsolverC)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverF)
    call linsol_releaseSolver (p_rsolverC)

    ! Release the 2-Level matrices
    call lsyssc_releaseMatrix (rmatrix2LvlT)
    call lsyssc_releaseMatrix (rmatrix2Lvl)

    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempF)
    call lsysbl_releaseVector (rtempC)
    call lsysbl_releaseVector (rsolF)
    call lsysbl_releaseVector (rsolC)
    call lsysbl_releaseVector (rrhsF)
    call lsysbl_releaseVector (rrhsC)
    call lsysbl_releaseMatrix (rmatrixF)
    call lsysbl_releaseMatrix (rmatrixC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscrF)
    call spdiscr_releaseBlockDiscr(rdiscrC)
    
    ! Release the triangulation.
    call tria_done (rtriaF)
    call tria_done (rtriaC)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

end module

!##############################################################################
!# Tutorial 021a: Solve a simple nonlinear system using the Picard iteration
!##############################################################################
!# Solves the system
!#
!#    -Laplace(u) - u + u^3 = f
!#
!# with a Newton iteration.
!#
!# Applies the block assembly for the defect and the Picard matrix.
!# No usage of template matrices. Linear systems are solved using
!# the Gauss elimination.
!##############################################################################

module tutorial021a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use derivatives
  use element
  use cubature
  use spatialdiscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  
  use discretebc
  use bcassembly
  use meshregion
  use vectorfilters
  use matrixfilters
  
  use linearsolver
  use collection
  
  use ucd

  implicit none
  private
  
  public :: start_tutorial021a

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fgetBoundaryValuesMR (Icomponents,rdiscretisation,rmeshRegion,&
      cinfoNeeded,Iwhere,Dwhere,Dcoords,Dvalues,rcollection)

!<description>
  ! Auxiliary function that returns values on the boundary.
!</description>

!<input>
  ! Component specifier.
  integer, dimension(:), intent(in) :: Icomponents

  ! Discretisation structure of the underlying equation.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in) :: rmeshRegion

  ! The type of information, the routine should calculate.
  integer, intent(in) :: cinfoNeeded

  ! Information about the DOF to be processed
  integer, dimension(4), intent(in) :: Iwhere

  ! The coordinates of the point which is currently processed.
  real(DP), dimension(:), intent(in) :: Dwhere

  ! The coordinates of the point for which the boundary values are to be calculated.
  real(DP), dimension(:), intent(in) :: Dcoords
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional information.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! Dvalues(1) receives the calculated value on the boundary.
  real(DP), dimension(:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! Zero boundary values.
    Dvalues(1) = 0.0_DP

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fcalcDefect(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the nonlinear defect.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasIx, dbasIy, du, dux, duy, dx, dy, df
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    real(DP), dimension(:,:,:), pointer :: p_Du
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(1)
    p_DlocalVector => p_rvectorData%p_Dentry
    p_DbasTest => p_rvectorData%p_DbasTest
    
    ! Function u in the cubature points
    p_Du => revalVectors%p_RvectorData(1)%p_Ddata
  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Right-hand side f
        df = 100.0_DP*sin(2.0_DP*SYS_PI*dx)*sin(2.0_DP*SYS_PI*dy)
        
        ! Function value / derivative of u
        du  = p_Du(icubp,iel,DER_FUNC2D)
        dux = p_Du(icubp,iel,DER_DERIV2D_X)
        duy = p_Du(icubp,iel,DER_DERIV2D_Y)
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Psi_i:
        do idofe=1,p_rvectorData%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Psi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
          dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
          
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          !
          ! def  =  f - [-Laplace u - u + u^3]
          p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
              p_DcubWeight(icubp,iel) * &
              ( df * dbasI  -  dux * dbasIx  -    duy * dbasIy  &
                            +   du * dbasI   -  du**3 * dbasI )
          
        end do ! idofe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fcalcMatrix(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Picard iteration matrix of the operator 
    ! "F(u) = -Laplace(u) - u + u^3 - f".
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    real(DP) :: dbasI, dbasJ, du
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:,:), pointer :: p_Du

    ! =====================================================
    ! We have the operator
    !
    !    F(u) = -Laplace(u) - u + u^3 - f
    !         = -Laplace(u) - u + u^2 u - f
    !         = G(u)u
    ! 
    ! with the Picard operator
    !
    !    G(u)v = -Laplace(v) - v + u^2 v - f
    !
    ! being linear in v for fixed u. Thus, the Picard iteration
    ! matrix reads
    !
    !    A(u) = -Laplace - M + " u^2*M "
    !
    ! with the last term understood as the matrix from the bilinear
    ! form "(u^2 phi_j, phi_i)".
    ! =====================================================

    ! Laplace
    call bma_docalc_laplace (RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,1.0_DP,1,1)
    
    ! Mass
    call bma_docalc_mass (RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,-1.0_DP,1,1)
    
    ! (u^2*phi_j, phi_i)
    p_DcubWeight => rassemblyData%p_DcubWeight
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      
    p_Du => revalVectors%p_RvectorData(1)%p_Ddata

    do iel = 1,nelements
      do icubp = 1,npointsPerElement
        du  = p_Du(icubp,iel,DER_FUNC2D)
        
        do idofe=1,RmatrixData(1,1)%ndofTest
          dbasI  = p_DbasTest(idofe,DER_FUNC2D,icubp,iel)
          
          do jdofe=1,RmatrixData(1,1)%ndofTrial
            dbasJ  = p_DbasTrial(jdofe,DER_FUNC2D,icubp,iel)

            p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * du**2 * dbasJ*dbasI

          end do
        end do
      end do
    end do

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial021a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rdefect, rsolution, rtemp
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_discreteBC) :: rdiscreteBC
    type(t_meshRegion) :: rmeshRegion
    type(t_ucdExport) :: rexport
    
    type(t_fev2Vectors) :: rcoeffVectors
    type(t_linsolNode), pointer :: p_rsolverNode
    integer :: ierror
    
    real(DP) :: depsRel,dresInit,dres
    integer :: ite

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 021a")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 11x11-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 10, 10)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rspatialDiscr
    ! which describes the discretisation.
    ! We generate a 1x1 system with Q1.
    ! =================================

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr,EL_Q1_2D,rtriangulation)
    
    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr,rtriangulation)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_commitBlockDiscr (rblockDiscr)

    ! =================================
    ! Assemble a matrix, a RHS and create 
    ! an empty solution vector.
    ! =================================

    ! Use a Gauss 3x3 formula for the discretisation.
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! -----------------------------------------------------
    ! Initialise a block matrix and block vectors
    ! corresponding to the block discretisation.
    call lsysbl_createMatrix (rblockDiscr,rmatrix)
    call lsysbl_createVector (rblockDiscr,rdefect)
    call lsysbl_createVector (rblockDiscr,rsolution)
    call lsysbl_createVector (rblockDiscr,rtemp)

    ! Create an empty matrix in block (1,1) in CSR format
    call bilf_createMatrixStructure (rmatrix,1,1,LSYSSC_MATRIX9)
    call lsysbl_allocEmptyMatrix (rmatrix)

    ! =================================
    ! Discretise boundary conditions
    ! =================================
    
    ! Initialise a boundary condition structure
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Get a mesh region for the complete boundary
    call mshreg_createFromNodalProp(rmeshRegion, rtriangulation, MSHREG_IDX_ALL)
    
    ! Discretise Dirichlet boundary conditions
    call bcasm_newDirichletBConMR (rblockDiscr, 1, rdiscreteBC, &
        rmeshRegion, fgetBoundaryValuesMR)
        
    call mshreg_done(rmeshregion)

    ! =================================
    ! Initialise linear solver
    ! =================================
    
    ! Create an UMFPACK solver.
    call linsol_initUMFPACK4 (p_rsolverNode)
    
    ! Attach the system matrix
    call linsol_setMatrix (p_rsolverNode, rmatrix)
    
    ! Symbolic factorisation
    call linsol_initStructure (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line ("Error during symbolic factorisation.")
      call sys_halt()
    end if

    ! =================================
    ! Newton iteration
    ! =================================

    ! Nonlinear solution vector. Impose boundary conditions.
    call lsysbl_clearVector (rsolution)
    call vecfil_discreteBCsol (rsolution,rdiscreteBC)

    ! Relative stopping criterion
    depsrel = 1E-10_DP
    
    ! Current iteration counter
    ite = 0
    
    do
      ! =========================================
      ! Nonlinear defect
      ! =========================================
      call lsysbl_clearVector (rdefect)

      ! Evaluate solution + 1st derivative in the calculation of the RHS vector.
      ! The first derivative is necessary for -Laplace(u) in the RHS.
      call fev2_addVectorToEvalList (rcoeffVectors,rsolution%RvectorBlock(1),1)
      call bma_buildVector (rdefect,BMA_CALC_STANDARD,&
          fcalcDefect, revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)
      call fev2_releaseVectorList(rcoeffVectors)

      ! Boundary conditions in the defect vector
      call vecfil_discreteBCdef (rdefect,rdiscreteBC)
          
      ! =========================================
      ! Residual check
      ! =========================================
      
      ! Get the norm of the residual
      dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
      
      if (ite .eq. 0) then
        dresInit = dres
      end if
      
      call output_line ("Ite "//trim(sys_siL(ite,10))//", ||RES|| = "//sys_sdEL(dres,10))
      
      ! Stopping criterion
      if (dres .le. dresInit * depsrel) exit

      ! We start the next iteration.
      ite = ite + 1
      
      ! =========================================
      ! Newton matrix of the operator
      ! =========================================
      call lsysbl_clearMatrix (rmatrix)

      ! Evaluate the solution (without derivative) during the calculation of
      ! the matrix. The solution is the nonlinearity.      
      call fev2_addVectorToEvalList (rcoeffVectors,rsolution%RvectorBlock(1),0)
      call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
          fcalcMatrix, revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfo)
      call fev2_releaseVectorList(rcoeffVectors)

      ! Boundary conditions
      call matfil_discreteBC (rmatrix,rdiscreteBC)

      ! =========================================
      ! Solve the linear subproblem
      ! =========================================
      
      ! Initialise the data of the solver
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) then
        call output_line("Matrix singular!",OU_CLASS_ERROR)
        call sys_halt()
      end if
      
      ! Solve the system or die trying
      call linsol_precondDefect (p_rsolverNode, rdefect)
      
      ! Release solver data and structure
      call linsol_doneData (p_rsolverNode)
      
      ! =========================================
      ! Update the nonlinear solution
      ! =========================================
      
      ! Update the nonlinear solution
      call lsysbl_vectorLinearComb (rdefect,rsolution,1.0_DP,1.0_DP)
      
      ! Re-impose boundary conditions -- to be on the safe side.
      call vecfil_discreteBCsol (rsolution,rdiscreteBC)
    
    end do

    ! =========================================
    ! Cleanup linear solver
    ! =========================================
    
    ! Symbolic data
    call linsol_doneStructure (p_rsolverNode)
    
    ! Remaining solver data
    call linsol_releaseSolver (p_rsolverNode)

    ! =================================
    ! Create output files
    ! =================================

    call output_line ("Writing postprocessing files...")

    ! Open / write / close; write the solution to a VTK file.
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,"post/tutorial021a.vtk")
    call ucd_addVectorByVertex (rexport, "solution", &
        UCD_VAR_STANDARD, rsolution%RvectorBlock(1))
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the BC
    call bcasm_releaseDiscreteBC(rdiscreteBC)
    
    ! Release the matrix/vectors
    call lsysbl_releaseMatrix (rmatrix)
    call lsysbl_releaseVector (rtemp)
    call lsysbl_releaseVector (rdefect)
    call lsysbl_releaseVector (rsolution)
    
    ! Release the discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Cubature done.
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module

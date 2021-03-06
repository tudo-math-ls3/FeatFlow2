!##############################################################################
!# Tutorial 020b: Solve a pure Neumann Poisson problem with UMFPACK
!##############################################################################

module tutorial020b

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use triangulation
  use meshgeneration

  use derivatives
  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation

  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  use linearsolver
  use filtersupport
  use collection

  use lumping
  use matrixmodification
  use vectorfilters

  use ucd

  implicit none
  private

  public :: start_tutorial020b

contains

  !****************************************************************************

!<subroutine>

  subroutine frhs(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=dscale*32*y*(1-y)+32*x*(1-x) - 32/3.
    !
    ! This RHS fulfils the property having integral mean value =0 on [0,1]^2
    ! which is a necessary condition for a RHS in a pure Neumann problem.
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
    real(DP) :: dbasI, dval, dx, dy
    integer :: iel, icubp, idofe, idimfe, ndimfe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get the data arrays of subvector 1
    p_rvectorData => RvectorData(1)
    p_DlocalVector => p_rvectorData%p_Dentry
    p_DbasTest => p_rvectorData%p_DbasTest

    ! FE space dimension
    ndimfe = p_rvectorData%ndimfe

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Calculate the values of the RHS using the coordinates
        ! of the cubature points.
        dval = 32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx) - 32.0_DP/3.0_DP

        ! Loop over the dimensions of the FE space
        do idimfe = 0,ndimfe-1

          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Psi_i:
          do idofe=1,p_rvectorData%ndofTest

            ! Fetch the contributions of the (test) basis functions Psi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe+idimfe*p_rvectorData%ndofTest,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI

          end do ! idofe

        end do ! idimfe

      end do ! icubp

    end do ! iel

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial020b

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_matrixBlock) :: rmatrix,rmass
    type(t_matrixScalar) :: rtemplateFEM,rmassLumped
    type(t_vectorBlock) :: rrhs, rsolution, rtemp
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_ucdExport) :: rexport
    character(LEN=SYS_STRLEN) :: spostdir

    type(t_filterChain), dimension(1), target :: RfilterChain
    integer :: nfilters
    type(t_linsolNode), pointer :: p_rsolverNode
    integer :: ierror

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 020b")
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
    ! Create a template FEM matrix which provides the structure
    ! for all matrices.
    call bilf_createMatrixStructure (rspatialDiscr,LSYSSC_MATRIX9,rtemplateFEM)

    ! -----------------------------------------------------
    ! Enlarge the template FEM matrix such that the first row is full.
    ! Later, we will impose the lumped mass matrix into that row.
    call mmod_expandToFullRow (rtemplateFEM,1)

    ! -----------------------------------------------------
    ! Initialise a block matrix and block vectors
    ! corresponding to the block discretisation.
    call lsysbl_createVector (rblockDiscr,rrhs)
    call lsysbl_createVector (rblockDiscr,rsolution)
    call lsysbl_createVector (rblockDiscr,rtemp)

    ! -----------------------------------------------------
    ! Create an empty matrix in block (1,1) based on the template matrix
    call lsysbl_createMatrix (rblockDiscr,rmatrix)

    call lsysbl_duplicateMatrix (rtemplateFEM,rmatrix,1,1,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    ! Discretise a simple Laplace matrix into (1,1)
    call lsysbl_clearMatrix (rmatrix)
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
        bma_fcalc_laplace,rcubatureInfo=rcubatureInfo)

    ! -----------------------------------------------------
    ! Discretise a mass matrix. Take the structure of the
    ! Laplace matrix. Allocate new memory.
    call lsysbl_createMatrix (rblockDiscr,rmass)

    call lsysbl_duplicateMatrix (rtemplateFEM,rmass,1,1,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    ! Discretise into (1,1)
    call lsysbl_clearMatrix (rmass)
    call bma_buildMatrix (rmass,BMA_CALC_STANDARD,&
        bma_fcalc_mass,rcubatureInfo=rcubatureInfo)

    ! -----------------------------------------------------
    ! Create a lumped mass matrix by summing up the rows to the diagonal.
    call lsyssc_createDiagMatrixStruc (rspatialDiscr,LSYSSC_MATRIXD,rmassLumped)
    call lump_sumToDiagonal(rmass%RmatrixBlock(1,1),rmassLumped)

    ! -----------------------------------------------------
    ! Create a RHS vector into block (1)
    call lsysbl_clearVector (rrhs)
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
        frhs,rcubatureInfo=rcubatureInfo)

    ! =================================
    ! Prepare imposing the side condition int(u(1)) = 0
    ! =================================

    ! -----------------------------------------------------
    ! Insert the lumped mass matrix as first row to the system matrix.
    ! Overwrite the first row, which is anyway linearly dependent on the
    ! other rows. This realises the calculation of the int(u(1)).
    call mmod_replaceLineByLumpedMass (rmatrix%RmatrixBlock(1,1),1,rmassLumped)

    ! Replace the first row in the first block of in the RHS vector, f[1](1)
    ! by zero. In combination with the first row of the matrix realising the
    ! integral, this imposes the condition int(u(1))=0.
    call vecfil_oneEntryZero (rrhs,1,1)

    ! =================================
    ! Solve the system with Gauss elimination
    ! =================================

    call output_line ("Solving linear system...")

    ! ----------------
    ! Solver preparation
    ! ----------------

    ! Create a CG solver with the filter chain attached.
    call linsol_initUMFPACK4 (p_rsolverNode)

    ! Attach the system matrix
    call linsol_setMatrix (p_rsolverNode, rmatrix)

    ! Symbolic factorisation
    call linsol_initStructure (p_rsolverNode, ierror)

    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line ("Error during symbolic factorisation.")
      call sys_halt()
    end if

    ! Numeric factorisation
    call linsol_initData (p_rsolverNode, ierror)

    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line ("Error during numeric factorisation.")
      call sys_halt()
    end if

    ! ----------------
    ! Solve the system
    ! ----------------

    ! Clear the solution
    call lsysbl_clearVector (rsolution)

    ! Solve
    call linsol_solveAdaptively (p_rsolverNode,rsolution,rrhs,rtemp)

    ! ----------------
    ! Cleanup
    ! ----------------

    ! Numeric data
    call linsol_doneData (p_rsolverNode)

    ! Symbolic data
    call linsol_doneStructure (p_rsolverNode)

    ! Remaining solver data
    call linsol_releaseSolver (p_rsolverNode)

    ! Filter chain
    call filter_doneFilterChain (RfilterChain,nfilters)

    ! =================================
    ! Create output files
    ! =================================

    call output_line ("Writing postprocessing files...")

    ! Open / write / close; write the solution to a VTK file.
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(spostdir)//"/tutorial020b.vtk")
    else
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial020b.vtk")
    end if

    call ucd_addVectorByVertex (rexport, "solution", &
        UCD_VAR_STANDARD, rsolution%RvectorBlock(1))
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================

    ! Release the matrix/vectors
    call lsysbl_releaseMatrix (rmass)
    call lsysbl_releaseMatrix (rmatrix)

    call lsyssc_releaseMatrix (rmassLumped)
    call lsyssc_releaseMatrix (rtemplateFEM)

    call lsysbl_releaseVector (rtemp)
    call lsysbl_releaseVector (rrhs)
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

!##############################################################################
!# ****************************************************************************
!# <name> statoptc_method0 </name>
!# ****************************************************************************
!#
!# <purpose>
!# Solves a stationary optimal control problem
!#
!#    J(u) = 1/2 ||y-z||^2  +  alpha/2 ||u||^2
!#
!#        -laplace(y) + y^2 = u
!#
!# which, using the Lagrange multiplier method, gives the system
!#
!#        -laplace(y) + y^2   = u
!#        -laplace(l) + 2 y l = y-z
!#
!#    C(u) :=   alpha u + l = 0
!#
!# Using the integral method approach, a Newton method is set up
!# for the control equation
!#
!#    u_n+1  =  u_n  +  C'(u_n)^-1 ( - alpha u_n - l_n )
!#           =  u_n  +  u~
!#
!# with u~ being the solution of the linearised
!# system
!#
!#    C'(u_n) u~  =  f  :=  - alpha u_n - l_n
!#
!# which is
!#
!#    -laplace(y~) + 2 y~   = u~
!#    -laplace(l~) + 2 y l~ = y~ - 2 l y~
!#
!#            alpha u~ + l~ = f
!# 
!# </purpose>
!##############################################################################

module statoptc_method0

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use element
  use ucd
  use pprocerror
  use linearalgebra
  use stdoperators
  
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use collection
    
  use poisson2d_callback
  
  implicit none

  ! Assembly templates. Temmplate matrices and structures.
  type t_asmTemplates

    ! Boundary definition
    type(t_boundary) :: rboundary
    
    ! Mesh deifnition
    type(t_triangulation) :: rtriangulation
  
    ! Block discretisation of the space
    type(t_blockDiscretisation) :: rblockDiscr
    
    ! Cubature structure
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! Structure matrix
    type(t_matrixScalar) :: rmatrixFEM

    ! Laplace matrix
    type(t_matrixScalar) :: rmatrixLaplace

    ! Mass matrix
    type(t_matrixScalar) :: rmatrixMass
    
    ! Boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
    
    ! Temporary matrix with the structure of the FEM space.
    type(t_matrixBlock) :: rmatrixTemp
    
  end type  

contains

  ! ***************************************************************************

!<subroutine>

  subroutine getzerobd (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                        cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! Returns zero boundary conditions.
!</description>
  
!<input>
  integer, dimension(:), intent(in)                           :: Icomponents
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  integer, intent(in)                                         :: ielement
  integer, intent(in)                                         :: cinfoNeeded
  integer, intent(in)                                          :: iwhere
  real(DP), intent(in)                                        :: dwhere
  type(t_collection), intent(inout), optional                 :: rcollection
!</input>

!<output>
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! REAL(DP) :: dx,dy
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP
  
  end subroutine

  ! ***************************************************************************
  ! Generates the template structures
  
  subroutine init_templates (sboundary,smesh,nlmax,rtemplates)

!<output>  
  type(t_asmTemplates), intent(out) :: rtemplates
!</output>

!<input>
  character(len=*), intent(in) :: sboundary
  character(len=*), intent(in) :: smesh
  integer, intent(in) :: nlmax
!</input>
  
    ! local variables
    type(t_boundaryRegion) :: rboundaryRegion
  
    ! -------------------
    ! Domain
    ! -------------------
  
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rtemplates%rboundary, sboundary)
        
    ! -------------------
    ! Triangulation
    ! -------------------
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtemplates%rtriangulation, smesh, rtemplates%rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (nlmax-1,rtemplates%rtriangulation,rtemplates%rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtemplates%rtriangulation,rtemplates%rboundary)
  
    ! --------------
    ! Discretisation
    ! --------------
  
    ! Set up a Q1 space
    call spdiscr_initBlockDiscr (rtemplates%rblockDiscr,1,&
        rtemplates%rtriangulation, rtemplates%rboundary)
    call spdiscr_initDiscr_simple (rtemplates%rblockDiscr%RspatialDiscr(1), &
        EL_Q2,rtemplates%rtriangulation,rtemplates%rboundary)
  
    ! -------------------
    ! Boundary conditions
    ! -------------------
  
    ! Set up boundary conditions
    call bcasm_initDiscreteBC(rtemplates%rdiscreteBC)
    call boundary_createRegion(rtemplates%rboundary,1,0,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rtemplates%rblockDiscr,1,&
        rboundaryRegion,rtemplates%rdiscreteBC,getzerobd)

    ! -------------------
    ! Cubature data
    ! -------------------
    call spdiscr_createDefCubStructure(rtemplates%rblockDiscr%RspatialDiscr(1), &
        rtemplates%rcubatureInfo, CUB_GEN_AUTO_G4)

    ! -------------------
    ! Template matrices
    ! -------------------

    ! Create a template matrix
    call bilf_createMatrixStructure (rtemplates%rblockDiscr%RspatialDiscr(1),&
        LSYSSC_MATRIX9,rtemplates%rmatrixFEM)
    
    ! Support Neumann problems
    ! call mmod_expandToFullRow(rtemplates%rmatrixFEM,1)

    ! Create Laplace, Mass,... based on the FEM matrix
    call lsyssc_duplicateMatrix (rtemplates%rmatrixFEM,rtemplates%rmatrixMass,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call lsyssc_duplicateMatrix (rtemplates%rmatrixFEM,rtemplates%rmatrixLaplace,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
    call lsyssc_clearMatrix (rtemplates%rmatrixMass)
    call lsyssc_clearMatrix (rtemplates%rmatrixLaplace)
    
    call stdop_assembleLaplaceMatrix(rtemplates%rmatrixLaplace, &
        .false., 1.0_DP, rtemplates%rcubatureInfo)

    call stdop_assembleSimpleMatrix(rtemplates%rmatrixMass, DER_FUNC,&
        DER_FUNC,1.0_DP, .false., rtemplates%rcubatureInfo)

    ! Temporary block matrix
    call lsysbl_createMatBlockByDiscr (rtemplates%rblockDiscr,rtemplates%rmatrixTemp)
    call lsyssc_duplicateMatrix (rtemplates%rmatrixFEM,rtemplates%rmatrixTemp%RmatrixBlock(1,1),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

  end subroutine

  ! ***************************************************************************
  ! Clean up the template structures
  
  subroutine done_templates (rtemplates)
  
  type(t_asmTemplates), intent(inout) :: rtemplates

    call boundary_release(rtemplates%rboundary)
    call tria_done(rtemplates%rtriangulation)
    call spdiscr_releaseBlockDiscr(rtemplates%rblockDiscr)
    call spdiscr_releaseCubStructure(rtemplates%rcubatureInfo)
    call lsyssc_releaseMatrix(rtemplates%rmatrixFEM)
    call lsyssc_releaseMatrix(rtemplates%rmatrixLaplace)
    call lsyssc_releaseMatrix(rtemplates%rmatrixMass)
    call lsysbl_releaseMatrix(rtemplates%rmatrixTemp)
    call bcasm_releaseDiscreteBC(rtemplates%rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  ! Solves a linear system
  
  subroutine calc_discretebcdef (rtemplates,rdefect)

!<input>
  type(t_asmTemplates), intent(in) :: rtemplates
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>

    call vecfil_discreteBCdef(rdefect,rtemplates%rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  ! Solves a linear system
  
  subroutine calc_discretebcsol (rtemplates,rdefect)

!<input>
  type(t_asmTemplates), intent(in) :: rtemplates
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>

    call vecfil_discreteBCsol(rdefect,rtemplates%rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  ! Solves a linear system
  
  subroutine calc_discretebcrhs (rtemplates,rdefect)

!<input>
  type(t_asmTemplates), intent(in) :: rtemplates
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>

    call vecfil_discreteBCrhs(rdefect,rtemplates%rdiscreteBC)

  end subroutine

  ! ***************************************************************************
  ! Solves a linear system
  
  subroutine solve_linsys (rmatrix,rdefect)

!<input>
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    type(t_linsolNode), pointer :: p_rlinsol
    integer :: ierror
    
    ! Solve
    call linsol_initUmfpack4 (p_rlinsol)
    call linsol_setMatrix (p_rlinsol,rmatrix)
    call linsol_initStructure (p_rlinsol,ierror)
    if (ierror .ne. 0) stop
    call linsol_initData (p_rlinsol,ierror)
    if (ierror .ne. 0) stop
    call linsol_precondDefect (p_rlinsol,rdefect)
    call linsol_doneData (p_rlinsol)
    call linsol_doneStructure (p_rlinsol)
    call linsol_releaseSolver (p_rlinsol)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_forward(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator at position (x,y) in a block matrix.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = x-coordinate in the block matrix
    ! rcollection%IquickAccess(2) = y-coordinate in the block matrix
!</description>

!<inputoutput>
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    real(DP) :: dbasJ, dbasI
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:,:), pointer :: p_Dy1

    real(DP) :: dy1

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the matrix data
    p_rmatrixData => RmatrixData(1,1)
    p_DlocalMatrix => RmatrixData(1,1)%p_Dentry
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial

    ! Nonlinear data
    p_Dy1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        dy1 = p_Dy1(icubp,iel)

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( dy1 * dbasJ*dbasI )

          end do ! idofe

        end do ! jdofe

      end do ! icubp

    end do ! iel

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_forwardlin(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator at position (x,y) in a block matrix.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = x-coordinate in the block matrix
    ! rcollection%IquickAccess(2) = y-coordinate in the block matrix
!</description>

!<inputoutput>
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    real(DP) :: dbasJ, dbasI
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:,:), pointer :: p_Dy1

    real(DP) :: dy1

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the matrix data
    p_rmatrixData => RmatrixData(1,1)
    p_DlocalMatrix => RmatrixData(1,1)%p_Dentry
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial

    ! Nonlinear data
    p_Dy1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        dy1 = p_Dy1(icubp,iel)

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( 2.0_DP*dy1 *dbasJ*dbasI )

          end do ! idofe

        end do ! jdofe

      end do ! icubp

    end do ! iel

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_backward(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator at position (x,y) in a block matrix.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = x-coordinate in the block matrix
    ! rcollection%IquickAccess(2) = y-coordinate in the block matrix
!</description>

!<inputoutput>
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    real(DP) :: dbasJ, dbasI, dx, dy
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:,:), pointer :: p_Dy1
    real(DP), dimension(:,:,:), pointer :: p_Dpoints

    real(DP) :: dy1

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get the matrix data
    p_rmatrixData => RmatrixData(1,1)
    p_DlocalMatrix => RmatrixData(1,1)%p_Dentry
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial

    ! Nonlinear data
    p_Dy1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Values of the nonlinearity
        dy1 = p_Dy1(icubp,iel)
        !dy1 = 1.0_DP

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( 2.0_DP*dy1 *dbasJ*dbasI )

          end do ! idofe

        end do ! jdofe

      end do ! icubp

    end do ! iel

  end subroutine

  !****************************************************************************
  ! Matrix of the additional term in the linearised dual equation

!<subroutine>

  subroutine bma_fcalc_backwardlin(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<inputoutput>
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    real(DP) :: dbasJ, dbasI, dx, dy
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:,:), pointer :: p_Dy1,p_Dlambda1,p_Dylin1
    real(DP), dimension(:,:,:), pointer :: p_Dpoints

    real(DP) :: dy1,dlambda1,dylin1

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get the matrix data
    p_rmatrixData => RmatrixData(1,1)
    p_DlocalMatrix => RmatrixData(1,1)%p_Dentry
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial

    ! Nonlinear data
    p_Dlambda1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dy1 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)
    p_Dylin1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Values of the nonlinearity
        dy1 = p_Dy1(icubp,iel)
        dlambda1 = p_Dlambda1(icubp,iel)
        dylin1 = p_Dylin1(icubp,iel)

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( 2.0_DP*dy1 *dbasJ*dbasI )

          end do ! idofe

        end do ! jdofe

      end do ! icubp

    end do ! iel

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_backwardlinrhs(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<inputoutput>
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<subroutine>

    real(DP) :: dbasJ, dbasI, dx, dy
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:,:), pointer :: p_Dy1,p_Dlambda1,p_Dylin1
    real(DP), dimension(:,:,:), pointer :: p_Dpoints

    real(DP) :: dy1,dlambda1,dylin1

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get the matrix data
    p_rmatrixData => RmatrixData(1,1)
    p_DlocalMatrix => RmatrixData(1,1)%p_Dentry
    p_DbasTest => RmatrixData(1,1)%p_DbasTest
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial

    ! Nonlinear data
    p_Dlambda1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dy1 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)
    p_Dylin1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC)

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Values of the nonlinearity
        dy1 = p_Dy1(icubp,iel)
        dlambda1 = p_Dlambda1(icubp,iel)
        dylin1 = p_Dylin1(icubp,iel)

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( 2.0_DP*dlambda1 * dbasJ*dbasI )

          end do ! idofe

        end do ! jdofe

      end do ! icubp

    end do ! iel

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsTarget(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the negative target "-z".
!</description>

!<inputoutput>
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(1)
    p_DlocalVector => RvectorData(1)%p_Dentry
    p_DbasTest => RvectorData(1)%p_DbasTest
  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Calculate the values of the RHS using the coordinates
        ! of the cubature points.
        dval = -(dx+dy)
        
        ! DEBUG!!!
        !dval = 32.0_DP*dy*(1.0_DP-dy)+32.0_DP*dx*(1.0_DP-dx) &
        !      +256.0_DP*dx**2*(1.0_DP-dx)**2*dy**2*(1.0_DP-dy)**2
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
          
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval * dbasI
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsSemilin(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the semilinear RHS.
!</description>

!<inputoutput>
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    integer, intent(in) :: npointsPerElement
    integer, intent(in) :: nelements
    type(t_fev2Vectors), intent(in) :: revalVectors
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy, dy1, dlambda1, dylin1
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    real(DP), dimension(:,:), pointer :: p_Dylin
    real(DP), dimension(:,:), pointer :: p_Dlambda
    real(DP), dimension(:,:), pointer :: p_Dy
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Nonlinearity
    p_Dlambda => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dy => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)
    p_Dylin => revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC)
    
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(1)
    p_DlocalVector => RvectorData(1)%p_Dentry
    p_DbasTest => RvectorData(1)%p_DbasTest
  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
        
        dy1 = p_Dy(icubp,iel)
        dlambda1 = p_Dlambda(icubp,iel)
        dylin1 = p_Dylin(icubp,iel)

        ! Calculate the values of the RHS using the coordinates
        ! of the cubature points.
        dval = -2.0_DP*dlambda1*dylin1
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
          
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval * dbasI
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
      
  end subroutine

  ! ***************************************************************************
  ! Calculate the solution of the forward problem
  
  subroutine calc_forward (rtemplates,ry,ru,rdefect)

!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
  type(t_vectorBlock), intent(in) :: ru
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: ry
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    real(DP) :: dres
    type(t_fev2Vectors) :: revalVectors
    real(DP), dimension(:), pointer :: p_Dy,p_Du,p_Dlambda
    real(DP), dimension(:), pointer :: p_Ddef
    integer :: ite
    
    call lsysbl_getbase_double (ry,p_Dy)
    call lsysbl_getbase_double (ru,p_Du)

    call lsysbl_getbase_double (rdefect,p_Ddef)
    
    ! Nonlinear loop
    call lsysbl_clearVector (ry)
    call calc_discretebcsol (rtemplates,ry)
    
    ite = 0
    do
      ! Prepare the evaluation stucture
      call fev2_addVectorToEvalList(revalVectors,ry%RvectorBlock(1),0)

      ! Create a RHS based on ru.
      call lsyssc_scalarMatVec (rtemplates%rmatrixMass, &
          ru%RvectorBlock(1), rdefect%RvectorBlock(1), 1.0_DP, 0.0_DP)
      
      ! DEBUG!!!
      !call lsysbl_clearVector (rdefect)
      !call bma_buildVector (rdefect,BMA_CALC_STANDARD,&
      !    bma_fcalc_rhsTarget,rcubatureInfo=rtemplates%rcubatureInfo,&
      !    revalVectors=revalVectors)
      
      
      ! Set up the system matrix
      call lsysbl_clearMatrix (rtemplates%rmatrixTemp)
      call bma_buildMatrix (rtemplates%rmatrixTemp,BMA_CALC_STANDARD,&
          bma_fcalc_forward,rcubatureInfo=rtemplates%rcubatureInfo,&
          revalVectors=revalVectors)
      call lsyssc_matrixLinearComb (rtemplates%rmatrixLaplace,&
          rtemplates%rmatrixTemp%RmatrixBlock(1,1),&
          1.0_DP,1.0_DP,.false.,.false.,.true.,.true.)
      call matfil_discreteBC (rtemplates%rmatrixTemp,rtemplates%rdiscreteBC)
    
      ! Subtract the LHS
      call lsysbl_blockMatVec (rtemplates%rmatrixTemp, &
          ry, rdefect, -1.0_DP, 1.0_DP)
      
      ! Boundary conditions
      call calc_discretebcdef (rtemplates,rdefect)
      
      ! Check the defect, stopping criterion
      dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
      write (*,*) "  Forward, iteration ",ite,", res = ",dres
      if (dres .lt. 1E-15) exit

      ! Set up the system matrix of the linearised operator
      call lsysbl_clearMatrix (rtemplates%rmatrixTemp)
      call bma_buildMatrix (rtemplates%rmatrixTemp,BMA_CALC_STANDARD,&
          bma_fcalc_forwardlin,rcubatureInfo=rtemplates%rcubatureInfo,&
          revalVectors=revalVectors)
      call lsyssc_matrixLinearComb (rtemplates%rmatrixLaplace,&
          rtemplates%rmatrixTemp%RmatrixBlock(1,1),&
          1.0_DP,1.0_DP,.false.,.false.,.true.,.true.)
      call matfil_discreteBC (rtemplates%rmatrixTemp,rtemplates%rdiscreteBC)
    
      ! Solve the linear system.
      call solve_linsys (rtemplates%rmatrixTemp,rdefect)
      
      ! Solution update
      call lsysbl_vectorLinearComb (rdefect,ry,1.0_DP,1.0_DP)
      
      ite = ite + 1
      
      ! Cleanup
      call fev2_releaseVectorList(revalVectors)
    
    end do
  
  end subroutine

  ! ***************************************************************************
  ! Calculate the solution of the backward problem
  
  subroutine calc_backward (rtemplates,rlambda,ry,ru,rdefect)

!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
  type(t_vectorBlock), intent(in) :: ry
  type(t_vectorBlock), intent(in) :: ru
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rlambda
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    real(DP) :: dres
    type(t_fev2Vectors) :: revalVectors
    real(DP), dimension(:), pointer :: p_Dy,p_Du,p_Dlambda
    real(DP), dimension(:), pointer :: p_Ddef
    
    call lsysbl_getbase_double (ry,p_Dy)
    call lsysbl_getbase_double (ru,p_Du)
    call lsysbl_getbase_double (rlambda,p_Dlambda)

    call lsysbl_getbase_double (rdefect,p_Ddef)

    ! Prepare the linear system.
    call lsysbl_clearVector (rlambda)
    call calc_discretebcsol (rtemplates,rlambda)
    
    ! Prepare the evaluation stucture
    call fev2_addVectorToEvalList(revalVectors,ry%RvectorBlock(1),0)

    ! Create a RHS based on ry.
    call lsyssc_scalarMatVec (rtemplates%rmatrixMass, &
        ry%RvectorBlock(1), rdefect%RvectorBlock(1), 1.0_DP, 0.0_DP)
        
    ! Target modifier: -z.
    call bma_buildVector (rdefect,BMA_CALC_STANDARD,&
        bma_fcalc_rhsTarget,rcubatureInfo=rtemplates%rcubatureInfo,&
        revalVectors=revalVectors)
    
    ! Set up the system matrix
    call lsysbl_clearMatrix (rtemplates%rmatrixTemp)
    call lsyssc_matrixLinearComb (rtemplates%rmatrixLaplace,&
        rtemplates%rmatrixTemp%RmatrixBlock(1,1),&
        1.0_DP,1.0_DP,.false.,.false.,.true.,.true.)
    call bma_buildMatrix (rtemplates%rmatrixTemp,BMA_CALC_STANDARD,&
        bma_fcalc_backward,rcubatureInfo=rtemplates%rcubatureInfo,&
        revalVectors=revalVectors)
    call matfil_discreteBC (rtemplates%rmatrixTemp,rtemplates%rdiscreteBC)
  
    ! Subtract the LHS
    call lsysbl_blockMatVec (rtemplates%rmatrixTemp, &
        rlambda, rdefect, -1.0_DP, 1.0_DP)
    
    call calc_discretebcdef (rtemplates,rdefect)
    
    ! Check the defect, stopping criterion
    dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
    write (*,*) "  Backward: ",dres

    ! Solve the linear system.
    call solve_linsys (rtemplates%rmatrixTemp,rdefect)
    
    ! Solution update
    call lsysbl_vectorLinearComb (rdefect,rlambda,1.0_DP,1.0_DP)
    
    ! Cleanup
    call fev2_releaseVectorList(revalVectors)
    
  end subroutine

  ! ***************************************************************************
  ! Calculate the solution of the backward problem
  
  subroutine calc_forwardlin (rtemplates,rylin,rulin,rlambda,ry,ru,rdefect)

!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
  type(t_vectorBlock), intent(in) :: rulin
  type(t_vectorBlock), intent(in) :: rlambda
  type(t_vectorBlock), intent(in) :: ry
  type(t_vectorBlock), intent(in) :: ru
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rylin
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    real(DP) :: dres
    type(t_fev2Vectors) :: revalVectors
    real(DP), dimension(:), pointer :: p_Dy,p_Du
    real(DP), dimension(:), pointer :: p_Dylin,p_Dulin
    real(DP), dimension(:), pointer :: p_Ddef
    
    call lsysbl_getbase_double (ry,p_Dy)
    call lsysbl_getbase_double (ru,p_Du)
    call lsysbl_getbase_double (rylin,p_Dylin)
    call lsysbl_getbase_double (ru,p_Dulin)    

    call lsysbl_getbase_double (rdefect,p_Ddef)
    
    ! Prepare the linear system.
    call lsysbl_clearVector (rylin)
    call calc_discretebcdef (rtemplates,rylin)
    
    ! Prepare the evaluation stucture
    call fev2_addVectorToEvalList(revalVectors,ry%RvectorBlock(1),0)

    ! Create a RHS based on rulin.
    call lsyssc_scalarMatVec (rtemplates%rmatrixMass, &
        rulin%RvectorBlock(1), rdefect%RvectorBlock(1), 1.0_DP, 0.0_DP)
    
    ! Set up the system matrix
    call lsysbl_clearMatrix (rtemplates%rmatrixTemp)
    call lsyssc_matrixLinearComb (rtemplates%rmatrixLaplace,&
        rtemplates%rmatrixTemp%RmatrixBlock(1,1),&
        1.0_DP,1.0_DP,.false.,.false.,.true.,.true.)
    call bma_buildMatrix (rtemplates%rmatrixTemp,BMA_CALC_STANDARD,&
        bma_fcalc_forwardlin,rcubatureInfo=rtemplates%rcubatureInfo,&
        revalVectors=revalVectors)
    call matfil_discreteBC (rtemplates%rmatrixTemp,rtemplates%rdiscreteBC)
  
    ! Subtract the LHS
    call lsysbl_blockMatVec (rtemplates%rmatrixTemp, &
        rylin, rdefect, -1.0_DP, 1.0_DP)
    
    call calc_discretebcdef (rtemplates,rdefect)
    
    ! Check the defect, stopping criterion
    dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
    write (*,*) "    Forward linear: ",dres

    ! Solve the linear system.
    call solve_linsys (rtemplates%rmatrixTemp,rdefect)
    
    ! Solution update
    call lsysbl_vectorLinearComb (rdefect,rylin,1.0_DP,1.0_DP)
    
    ! Cleanup
    call fev2_releaseVectorList(revalVectors)
    
  end subroutine

  ! ***************************************************************************
  ! Calculate the solution of the backward problem
  
  subroutine calc_backwardlin (rtemplates,rlambdalin,rylin,rulin,rlambda,ry,ru,rdefect)

!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
  type(t_vectorBlock), intent(in) :: rylin
  type(t_vectorBlock), intent(in) :: rulin
  type(t_vectorBlock), intent(in) :: rlambda
  type(t_vectorBlock), intent(in) :: ry
  type(t_vectorBlock), intent(in) :: ru
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rlambdalin
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    real(DP) :: dres
    type(t_fev2Vectors) :: revalVectors
    real(DP), dimension(:), pointer :: p_Dy,p_Du,p_Dlambda
    real(DP), dimension(:), pointer :: p_Dylin,p_Dulin,p_Dlambdalin
    real(DP), dimension(:), pointer :: p_Ddef
    
    call lsysbl_getbase_double (ry,p_Dy)
    call lsysbl_getbase_double (ru,p_Du)
    call lsysbl_getbase_double (rlambda,p_Dlambda)
    call lsysbl_getbase_double (rylin,p_Dylin)
    call lsysbl_getbase_double (rulin,p_Dulin)
    call lsysbl_getbase_double (rlambdalin,p_Dlambdalin)
    
    call lsysbl_getbase_double (rdefect,p_Ddef)

    ! Prepare the linear system.
    call lsysbl_clearVector (rlambdalin)
    call calc_discretebcdef (rtemplates,rlambdalin)
    
    ! Prepare the evaluation stucture
    call fev2_addVectorToEvalList(revalVectors,rlambda%RvectorBlock(1),0)
    call fev2_addVectorToEvalList(revalVectors,ry%RvectorBlock(1),0)
    call fev2_addVectorToEvalList(revalVectors,rylin%RvectorBlock(1),0)

    ! Create a RHS based on rylin.
    call lsyssc_scalarMatVec (rtemplates%rmatrixMass, &
        rylin%RvectorBlock(1), rdefect%RvectorBlock(1), 1.0_DP, 0.0_DP)
    
!    call lsysbl_clearMatrix (rtemplates%rmatrixTemp)
!    call bma_buildMatrix (rtemplates%rmatrixTemp,BMA_CALC_STANDARD,&
!        bma_fcalc_backwardlinrhs,rcubatureInfo=rtemplates%rcubatureInfo,&
!        revalVectors=revalVectors)
!    call lsysbl_blockMatVec (rtemplates%rmatrixTemp, &
!        rylin, rdefect, -1.0_DP, 1.0_DP)

    call bma_buildVector (rdefect,BMA_CALC_STANDARD,&
        bma_fcalc_rhsSemilin,rcubatureInfo=rtemplates%rcubatureInfo,&
        revalVectors=revalVectors)
    
    ! Set up the system matrix
    call lsysbl_clearMatrix (rtemplates%rmatrixTemp)
    call lsyssc_matrixLinearComb (rtemplates%rmatrixLaplace,&
        rtemplates%rmatrixTemp%RmatrixBlock(1,1),&
        1.0_DP,1.0_DP,.false.,.false.,.true.,.true.)
    call bma_buildMatrix (rtemplates%rmatrixTemp,BMA_CALC_STANDARD,&
        bma_fcalc_backwardlin,rcubatureInfo=rtemplates%rcubatureInfo,&
        revalVectors=revalVectors)
    call matfil_discreteBC (rtemplates%rmatrixTemp,rtemplates%rdiscreteBC)
  
    ! Subtract the LHS
    call lsysbl_blockMatVec (rtemplates%rmatrixTemp, &
        rlambdalin, rdefect, -1.0_DP, 1.0_DP)
    
    call calc_discretebcdef (rtemplates,rdefect)
    
    ! Check the defect, stopping criterion
    dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
    write (*,*) "    Backward linear: ",dres

    ! Solve the linear system.
    call solve_linsys (rtemplates%rmatrixTemp,rdefect)
    
    ! Solution update
    call lsysbl_vectorLinearComb (rdefect,rlambdalin,1.0_DP,1.0_DP)
    
    ! Cleanup
    call fev2_releaseVectorList(revalVectors)
    
  end subroutine

  ! ***************************************************************************
  ! Calculate the residual in the control equation
  
  subroutine calc_resControlEquation (dalpha,rlambda,ru,rdefect)

!<input>
  real(DP), intent(in) :: dalpha
  type(t_vectorBlock), intent(in) :: rlambda
  type(t_vectorBlock), intent(in) :: ru
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    call lsysbl_vectorLinearComb (ru,rdefect,-dalpha,0.0_DP)
    call lsysbl_vectorLinearComb (rlambda,rdefect,-1.0_DP,1.0_DP)
    
  end subroutine

  ! ***************************************************************************
  ! Calculate the residual in the control equation
  
  subroutine calc_resControlEquationLin (dalpha,rrhs,rlambda,ru,rdefect)

!<input>
  real(DP), intent(in) :: dalpha
  type(t_vectorBlock), intent(in) :: rrhs
  type(t_vectorBlock), intent(in) :: rlambda
  type(t_vectorBlock), intent(in) :: ru
!</input>  
  
!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rdefect
!</inputoutput>
  
    ! local variables
    call lsysbl_vectorLinearComb (rrhs,rdefect,1.0_DP,0.0_DP)
    call lsysbl_vectorLinearComb (ru,rdefect,-dalpha,1.0_DP)
    call lsysbl_vectorLinearComb (rlambda,rdefect,-1.0_DP,1.0_DP)
    
  end subroutine

  ! ***************************************************************************
  ! Postprocessing
  
  subroutine postproc (rtemplates,ry,rlambda,ru)

!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
  type(t_vectorBlock), intent(in) :: ry
  type(t_vectorBlock), intent(in) :: rlambda
  type(t_vectorBlock), intent(in) :: ru
!</input>  

    ! local variables
    type(t_ucdExport) :: rexport
  
    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtemplates%rtriangulation,&
                       "gmv/sol.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "y", UCD_VAR_STANDARD, &
        ry%RvectorBlock(1))

    call ucd_addVectorByVertex (rexport, "lambda", UCD_VAR_STANDARD, &
        rlambda%RvectorBlock(1))

    call ucd_addVectorByVertex (rexport, "u", UCD_VAR_STANDARD, &
        ru%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
  end subroutine

  ! ***************************************************************************
  ! Determins the Newton direction associated to a rhs of the linearised
  ! integral equation

  subroutine solve_integraleqn_dirderiv (rtemplates,dalpha,ry,rlambda,ru,rrhs,rulin)
  
!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
  real(DP), intent(in) :: dalpha
  type(t_vectorBlock), intent(in) :: ry
  type(t_vectorBlock), intent(in) :: rlambda
  type(t_vectorBlock), intent(in) :: ru
  type(t_vectorBlock), intent(in) :: rrhs
!</input>

!<inputoutput>
  type(t_vectorBlock), intent(inout) :: rulin
!</inputoutput>

    ! local variables
    type(t_vectorBlock) :: rdefect
    type(t_vectorBlock) :: rylin
    type(t_vectorBlock) :: rlambdalin
    real(DP) :: dres, domega
    integer :: ite
    real(DP), dimension(:), pointer :: p_Dy,p_Du,p_Dlambda
    real(DP), dimension(:), pointer :: p_Dylin,p_Dulin,p_Dlambdalin
    real(DP), dimension(:), pointer :: p_Ddef
    
    domega = 0.3_DP/dalpha
    
    ! Create additional vectors.
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,rdefect,.true.)
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,rylin,.true.)
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,rlambdalin,.true.)
    
    call lsysbl_getbase_double (ry,p_Dy)
    call lsysbl_getbase_double (ru,p_Du)
    call lsysbl_getbase_double (rlambda,p_Dlambda)
    call lsysbl_getbase_double (rylin,p_Dylin)
    call lsysbl_getbase_double (rulin,p_Dulin)
    call lsysbl_getbase_double (rlambdalin,p_Dlambdalin)
    
    call lsysbl_getbase_double (rdefect,p_Ddef)

    ! Start the defect correction loop
    call lsysbl_clearVector (rulin)
    
    ite = 0
    write (*,*) "  Starting inner loop."
    
    do
    
      ! ---------------------------------------------------
      ! Step 1: Solve the promal eqn, dual eqn,
      ! determine the residual.
      ! ---------------------------------------------------
      call calc_forwardlin (rtemplates,rylin,rulin,rlambda,ry,ru,rdefect)
      call calc_backwardlin (rtemplates,rlambdalin,rylin,rulin,rlambda,ry,ru,rdefect)
      call calc_resControlEquationLin (dalpha,rrhs,rlambdalin,rulin,rdefect)

      ! ---------------------------------------------------
      ! Step 2: Stopping criterion
      ! ---------------------------------------------------
      dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
      write (*,*) "  Inner loop, iteration ",ite,", res = ",dres
      if (dres .lt. 1E-15) exit
      
      ! ---------------------------------------------------
      ! Step 3: Correction
      ! ---------------------------------------------------
      call lsysbl_vectorLinearComb (rdefect,rulin,domega,1.0_DP)
      
      ite = ite + 1
    
    end do
    
    call lsysbl_releaseVector (rdefect)
    call lsysbl_releaseVector (rylin)
    call lsysbl_releaseVector (rlambdalin)
    
  end subroutine

  ! ***************************************************************************
  ! Solve the problem.

  subroutine solve_integraleqn_newton (rtemplates,ry,rlambda,ru)
  
!<input>
  type(t_asmTemplates), intent(inout) :: rtemplates
!</input>

!<inputoutput>
  type(t_vectorBlock), intent(inout) :: ru
  type(t_vectorBlock), intent(inout) :: ry
  type(t_vectorBlock), intent(inout) :: rlambda
!</inputoutput>

    ! local variables
    type(t_vectorBlock) :: rdefect
    type(t_vectorBlock) :: rulin
    real(DP) :: dres,dalpha
    integer :: ite
    real(DP), dimension(:), pointer :: p_Dy,p_Du,p_Dlambda
    real(DP), dimension(:), pointer :: p_Ddef
    
    dalpha = 0.001_DP
    
    ! Create additional vectors.
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,rdefect,.true.)
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,rulin,.true.)
    
    call lsysbl_getbase_double (ry,p_Dy)
    call lsysbl_getbase_double (ru,p_Du)
    call lsysbl_getbase_double (rlambda,p_Dlambda)
    
    call lsysbl_getbase_double (rdefect,p_Ddef)
    
    ! Start the Newton iteration
    call lsysbl_clearVector (ru)
    ite = 0
    write (*,*) "Starting outer loop."
    do
    
      ! ---------------------------------------------------
      ! Step 1: Solve the promal eqn, dual eqn,
      ! determine nonlinear residual.
      ! ---------------------------------------------------
      call calc_forward (rtemplates,ry,ru,rdefect)
      call calc_backward (rtemplates,rlambda,ry,ru,rdefect)
      call calc_resControlEquation (dalpha,rlambda,ru,rdefect)
    
      ! ---------------------------------------------------
      ! Step 2: Stopping criterion
      ! ---------------------------------------------------
      dres = lsysbl_vectorNorm (rdefect,LINALG_NORML2)
      write (*,*) "Outer loop, iteration ",ite,", res = ",dres
      if (dres .lt. 1E-14) exit

      ! ---------------------------------------------------
      ! Step 3: Solve the linearised equation
      ! ---------------------------------------------------
      call solve_integraleqn_dirderiv (rtemplates,dalpha,ry,rlambda,ru,rdefect,rulin)

      ! ---------------------------------------------------
      ! Step 4: Correction
      ! ---------------------------------------------------
      call lsysbl_vectorLinearComb (rulin,ru,1.0_DP,1.0_DP)
      
      ite = ite + 1

    end do
  
    call lsysbl_releaseVector (rdefect)
    call lsysbl_releaseVector (rulin)
    
  end subroutine

  ! ***************************************************************************
  ! Solve the problem.

  subroutine statoptc_0
  
    ! local variables
    type(t_asmTemplates) :: rtemplates
    type(t_vectorBlock) :: ry
    type(t_vectorBlock) :: rlambda
    type(t_vectorBlock) :: ru
    integer :: NLMAX
  
    NLMAX = 4
  
    ! Initialisation
    call init_templates ("./pre/QUAD.prm","./pre/QUAD.tri",NLMAX,rtemplates)
    
    ! Create vectors
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,ru,.true.)
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,ry,.true.)
    call lsysbl_createVectorBlock (rtemplates%rblockDiscr,rlambda,.true.)
    
    ! Invoke the solver
    call solve_integraleqn_newton (rtemplates,ry,rlambda,ru)
    
    ! Postprocessing
    call postproc (rtemplates,ry,rlambda,ru)
    
    ! Cleanup
    call lsysbl_releaseVector (ry)
    call lsysbl_releaseVector (rlambda)
    call lsysbl_releaseVector (ru)
    call done_templates (rtemplates)
    
  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_smart </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson2d_method0_simple, however, it utilises
!# new kernel features to compress the code, especially regarding the assembly
!# of the operator, the right-hand-side and the boundary conditions.
!# </purpose>
!##############################################################################

module poisson2d_method0_smart

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use derivatives
  use element
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
  use stdoperators

  use poisson2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_0_smart
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to VTK file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let us see...
    !
    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir, sucddir

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtria

    ! An object specifying a Q1 discretisation.
    type(t_spatialDiscretisation) :: rdiscrQ1
    
    ! An object specifying the discretisation of the equation.
    type(t_blockDiscretisation) :: rdiscr

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubature
    
    ! A matrix, a RHS vector, a solution vector and a temporary vector.
    type(t_matrixBlock) :: rmatSystem
    type(t_vectorBlock) :: rvecSol,rvecRhs,rvecTmp

    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rregion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolver

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 7

    ! Fetch the system variables for the "pre" and "ucd" directories
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the boundary, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Read in the boundary
    call boundary_read_prm(rboundary, trim(spredir) // "/QUAD.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtria, trim(spredir) // "/QUAD.tri", rboundary)
     
    ! Refine it up to the desired level.
    call tria_quickRefine2LevelOrdering (NLMAX-1, rtria, rboundary)
    
    ! And initialise a standard mesh from it.
    call tria_initStandardMeshFromRaw (rtria, rboundary)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Initialise a simple Q1 space.
    call spdiscr_initDiscr_simple (rdiscrQ1, EL_Q1_2D, rtria, rboundary)
   
    ! Create a block discretisation for the equation - one Q1 component.
    call spdiscr_initBlockDiscr (rdiscr, rtria, rboundary)
    call spdiscr_appendBlockComponent (rdiscr, rdiscrQ1)
    call spdiscr_commitBlockDiscr (rdiscr)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create a cubature information structure: 3x3 Gauss rule
    call spdiscr_createDefCubStructure(rdiscr%RspatialDiscr(1), rcubature, CUB_GEN_AUTO_G3)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 1x1 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Create a block matrix based on the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscr, rmatSystem)
    
    ! Now assemble the matrix structure.
    call bilf_createMatrixStructure (rdiscr%RspatialDiscr(1), LSYSSC_MATRIX9, &
        rmatSystem%RmatrixBlock(1,1))
    
    ! And assemble a Laplace matrix
    call stdop_assembleLaplaceMatrix(rmatSystem%RmatrixBlock(1,1), .true., 1.0_DP, rcubature)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Initialise a discrete BC structure
    call bcasm_initDiscreteBC(rdiscreteBC)

    ! Fetch the boundary`s region. Passing 0 as the third parameter will give us a
    ! region representing the whole boundary component instead of a single segment.
    call boundary_createRegion(rboundary, 1, 0, rregion)

    ! Now assemble homogene Dirichlet BCs on that region
    call bcasm_newHomDirichletBConRealBd (rdiscr, 1, rregion, rdiscreteBC)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create 3 block vectors based on the matrix:
    ! the RHS vector, the solution vector and a temporary vector
    call lsysbl_createVectorBlock (rmatSystem, rvecRhs, .true.)
    call lsysbl_createVectorBlock (rmatSystem, rvecSol, .true.)
    call lsysbl_createVectorBlock (rmatSystem, rvecTmp, .true.)

    ! Now assemble the right-hand-side vector
    call linf_buildSimpleVector(rvecRhs%RvectorBlock(1), rcubature, coeff_RHS_2D)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Filter the system; this will impose the previously assembled boundary conditions
    call matfil_discreteBC (rmatSystem,rdiscreteBC)
    call vecfil_discreteBCrhs (rvecRhs,rdiscreteBC)
    call vecfil_discreteBCsol (rvecSol,rdiscreteBC)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create a Gauss elimination solver without any preconditioner.
    call linsol_initUmfpack4 (p_rsolver)
    
    ! Attach the system matrix to the solver.
    call linsol_setMatrix(p_rsolver, rmatSystem)

    ! Initialise structure of the solver and validate the output.
    call linsol_initStructure (p_rsolver, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if

    ! Initialise the data of the solver
    call linsol_initData (p_rsolver, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix singular!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solve the system
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Solve the system or die trying
    call linsol_solveAdaptively (p_rsolver, rvecSol, rvecRhs, rvecTmp)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Export the solution to VTK
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Start UCD export to VTK file:
    call ucd_startVTK(rexport, UCD_FLAG_STANDARD, rtria, trim(sucddir) // "/u2d_0_smart.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, rvecSol%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Perform error analysis
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR, derror, rvecSol%RvectorBlock(1), &
                       getReferenceFunction_2D, rcubatureInfo=rcubature)
    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR, derror, rvecSol%RvectorBlock(1), &
                       getReferenceFunction_2D, rcubatureInfo=rcubature)
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolver)
    call linsol_doneStructure (p_rsolver)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolver)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecSol)
    call lsysbl_releaseVector (rvecRhs)
    call lsysbl_releaseMatrix (rmatSystem)

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rcubature)

    ! Release the discretisation structures.
    call spdiscr_releaseBlockDiscr(rdiscr)
    call spdiscr_releaseDiscr(rdiscrQ1)
    
    ! Release the triangulation.
    call tria_done (rtria)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
  end subroutine

end module

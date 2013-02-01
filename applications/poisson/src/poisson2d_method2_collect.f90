!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method2_collect </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using a collection structure that saves problem-dependent data.
!# </purpose>
!##############################################################################

module poisson2d_method2_collect

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
  use ucd
  use pprocerror
  
  use collection
    
  use poisson2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_initParamTriang (ilv,rcollection)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<inputoutput>
  ! The level up to where we refine the coarse mesh.
  ! If the number is too large, ilv is reduced to the maximum allowed value.
  integer, intent(inout) :: ilv

  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to p_rboundary.
    allocate(p_rboundary)
    call boundary_read_prm(p_rboundary, trim(spredir)//"/QUAD.prm")

    ! Now read in the basic triangulation.
    allocate(p_rtriangulation)
    call tria_readTriFile2D (p_rtriangulation, trim(spredir)//"/QUAD.tri", p_rboundary)
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (ilv-1,p_rtriangulation,p_rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (p_rtriangulation,p_rboundary)
    
    ! The TRIAS(,)-array is now part pf the triangulation structure,
    ! we do not need it anymore.
    !
    ! *What* we need later is the definition of the boundary and the
    ! triangulation. Save these to the collection.
    call collct_setvalue_bdry(rcollection,"DOMAIN",p_rboundary,.true.)
    call collct_setvalue_tria(rcollection,"TRIA",p_rtriangulation,.true.)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_initDiscretisation (rcollection)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the collection.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! An object for the spatial discretisation
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Ask the collection to give us the boundary and triangulation.
    ! We need it for the discretisation.
    p_rboundary => collct_getvalue_bdry(rcollection,"DOMAIN")
    p_rtriangulation => collct_getvalue_tria(rcollection,"TRIA")
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    allocate(p_rdiscretisation)
    call spdiscr_initBlockDiscr (p_rdiscretisation,1,&
                                 p_rtriangulation, p_rboundary)
    
    ! p_rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple (p_rdiscretisation%RspatialDiscr(1), &
                                   EL_Q1,p_rtriangulation, p_rboundary)
                                   
    ! Add the discretisation structure to the collection so that
    ! we can use it later.
    call collct_setvalue_bldiscr(rcollection,"DISCR2D",p_rdiscretisation,.true.)
    
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    allocate(p_rcubatureInfo)    
    call spdiscr_createDefCubStructure(&  
        p_rdiscretisation%RspatialDiscr(1),p_rcubatureInfo,CUB_GEN_AUTO_G3)
    call collct_setvalue_cubinfo(rcollection,"CUBINFO2D",p_rcubatureInfo,.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_initMatVec (rcollection)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the collection.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form. The matrix will receive the discrete Laplace operator.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Ask the collection to give us the discretisation structure
    p_rdiscretisation => collct_getvalue_bldiscr(rcollection,"DISCR2D")

    ! Get the cubature information structure
    p_rcubatureInfo => collct_getvalue_cubinfo(rcollection,"CUBINFO2D")
    
    ! Create the matrix and the vectors on the heap
    allocate(p_rmatrix)
    allocate(p_rrhs)
    allocate(p_rvector)
    
    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)
    
    ! Save matrix and vectors to the collection.
    ! They maybe used later, expecially in nonlinear problems.
    call collct_setvalue_vec(rcollection,"RHS",p_rrhs,.true.)
    call collct_setvalue_vec(rcollection,"SOLUTION",p_rvector,.true.)
    call collct_setvalue_mat(rcollection,"LAPLACE",p_rmatrix,.true.)

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    call bilf_createMatrixStructure (&
              p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
              p_rmatrix%RmatrixBlock(1,1))
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical data.
    call bilf_buildMatrixScalar (rform,.true.,&
        p_rmatrix%RmatrixBlock(1,1),p_rcubatureInfo)
    
    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (p_rdiscretisation,p_rrhs,.true.)
    call lsysbl_createVectorBlock (p_rdiscretisation,p_rvector,.true.)

    ! The vector structure is done but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    call linf_buildVectorScalar (&
              rlinform,.true.,p_rrhs%RvectorBlock(1),p_rcubatureInfo,&
              coeff_RHS_2D)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_initDiscreteBC (rcollection)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    type(t_boundaryRegion) :: rboundaryRegion
    
    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC

    ! A pointer to the domain
    type(t_boundary), pointer :: p_rboundary
  
    ! Ask the collection to give us the discretisation structure and
    p_rdiscretisation => collct_getvalue_bldiscr(rcollection,"DISCR2D")
    
    ! Get the domain from the discretisation
    p_rboundary => p_rdiscretisation%p_rboundary

    ! Get our matrix and right hand side from the collection.
    p_rrhs    => collct_getvalue_vec(rcollection,"RHS")
    p_rvector => collct_getvalue_vec(rcollection,"SOLUTION")
    p_rmatrix => collct_getvalue_mat(rcollection,"LAPLACE")
    
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    allocate (p_rdiscreteBC)
    call bcasm_initDiscreteBC(p_rdiscreteBC)
    !
    ! We "know" already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    !
    ! We ask the boundary routines to create a "boundary region" - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent="1" to indicate that we set up the
    !   Dirichlet BC`s for the first (here: one and only) component in the
    !   solution vector.
    ! - Discretise the boundary condition so that the BC`s can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
        rboundaryRegion,p_rdiscreteBC,&
        getBoundaryValues_2D,rcollection)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
        rboundaryRegion,p_rdiscreteBC,&
        getBoundaryValues_2D,rcollection)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
        rboundaryRegion,p_rdiscreteBC,&
        getBoundaryValues_2D,rcollection)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
        rboundaryRegion,p_rdiscreteBC,&
        getBoundaryValues_2D,rcollection)
                             
    ! Save the structures to the collection, so we can access them
    ! without having the matrix or the vector.
    call collct_setvalue_discbc(rcollection,"DISCBC",p_rdiscreteBC,.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_implementBC (rcollection)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_discreteBC), pointer :: p_rdiscreteBC

    ! Get our matrix and right hand side from the collection.
    p_rrhs    => collct_getvalue_vec(rcollection,"RHS")
    p_rvector => collct_getvalue_vec(rcollection,"SOLUTION")
    p_rmatrix => collct_getvalue_mat(rcollection,"LAPLACE")
    p_rdiscreteBC => collct_getvalue_discbc(rcollection,"DISCBC")
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    call vecfil_discreteBCrhs (p_rrhs,p_rdiscreteBC)
    call vecfil_discreteBCsol (p_rvector,p_rdiscreteBC)
    call matfil_discreteBC (p_rmatrix,p_rdiscreteBC)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_solve (rcollection)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_vectorBlock), target :: rvecTmp

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner
    
    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => collct_getvalue_vec(rcollection,"RHS")
    p_rvector => collct_getvalue_vec(rcollection,"SOLUTION")
    p_rmatrix => collct_getvalue_mat(rcollection,"LAPLACE")
    
    ! Create a temporary vector for the solver - it needs that.
    call lsysbl_createVecBlockIndirect (p_rrhs, rvecTmp, .false.)
    
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.

    p_rdiscreteBC => collct_getvalue_discbc(rcollection,"DISCBC")

    call filter_clearFilterChain (RfilterChain,nfilters)
    call filter_newFilterDiscBCDef (RfilterChain,nfilters,p_rdiscreteBC)

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    nullify(p_rpreconditioner)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrix to the solver.
    call linsol_setMatrix(p_RsolverNode,p_rmatrix)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if

    call linsol_initData (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix singular!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rvecTmp)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the temporary vector
    call lsysbl_releaseVector (rvecTmp)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_postprocessing (rcollection)
  
!<description>
  ! Writes the solution into a VTK file.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir

    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! Error of FE function to reference function
    real(DP) :: derror

    ! Get the solution vector from the collection.
    p_rvector => collct_getvalue_vec(rcollection,"SOLUTION")
    
    ! Get the cubature information structure
    p_rcubatureInfo => collct_getvalue_cubinfo(rcollection,"CUBINFO2D")
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       trim(sucddir)//"/u2d_2_collect.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, &
        p_rvector%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR,derror,p_rvector%RvectorBlock(1),&
                       getReferenceFunction_2D, rcubatureInfo=p_rcubatureInfo)
    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR,derror,p_rvector%RvectorBlock(1),&
                       getReferenceFunction_2D, rcubatureInfo=p_rcubatureInfo)
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_doneMatVec (rcollection)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! Get our matrix and right hand side from the collection.
    p_rrhs    => collct_getvalue_vec(rcollection,"RHS")
    p_rvector => collct_getvalue_vec(rcollection,"SOLUTION")
    p_rmatrix => collct_getvalue_mat(rcollection,"LAPLACE")

    ! Release them from memory
    call lsysbl_releaseVector (p_rvector)
    call lsysbl_releaseVector (p_rrhs)
    call lsysbl_releaseMatrix (p_rmatrix)
    deallocate(p_rmatrix)
    deallocate(p_rrhs)
    deallocate(p_rvector)

    ! Delete the variables from the collection.
    call collct_deletevalue (rcollection,"RHS")
    call collct_deletevalue (rcollection,"SOLUTION")
    call collct_deletevalue (rcollection,"LAPLACE")

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine pm2_doneBC (rcollection)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables

    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC
    
    ! Get pointers to the discrete and analytic boundary conditions
    ! from the collection structure.
    p_rdiscreteBC => collct_getvalue_discbc(rcollection,"DISCBC")

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (p_rdiscreteBC)
    deallocate(p_rdiscreteBC)

    ! Delete the variables from the collection.
    call collct_deletevalue (rcollection,"BDCOND")
    call collct_deletevalue (rcollection,"DISCBC")

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine pm2_doneDiscretisation (rcollection)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
      ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    ! Ask the collection to give us the discretisation structure
    p_rdiscretisation => collct_getvalue_bldiscr(rcollection,"DISCR2D")
    
    ! Delete the discretisation...
    call spdiscr_releaseBlockDiscr(p_rdiscretisation)
    
    ! remove the allocated block discretisation structure
    deallocate(p_rdiscretisation)
    
    ! and remove it from the collection.
    call collct_deletevalue (rcollection,"DISCR2D")

    ! A pointer to the cubature information structure with the data.    
    p_rcubatureInfo => collct_getvalue_cubinfo(rcollection,"CUBINFO2D")
    
    ! Delete the structure...
    call spdiscr_releaseCubStructure(p_rcubatureInfo)
    
    ! remove the allocated block discretisation structure
    deallocate(p_rcubatureInfo)
    
    ! and remove it from the collection.
    call collct_deletevalue (rcollection,"CUBINFO2D")
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine pm2_doneParamTriang (rcollection)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A collection object for saving structural data and some problem-dependent
  ! information.
  type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!</subroutine>

    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation

    ! Ask the collection to give us the boundary and triangulation.
    ! We need it for the discretisation.
    p_rboundary => collct_getvalue_bdry(rcollection,"DOMAIN")
    p_rtriangulation => collct_getvalue_tria(rcollection,"TRIA")
    
    ! Release the triangulation.
    call tria_done (p_rtriangulation)
    deallocate(p_rtriangulation)
    
    ! Finally release the domain.
    call boundary_release (p_rboundary)
    deallocate(p_rboundary)
    
    ! Remove everything from the collection.
    call collct_deletevalue (rcollection,"TRIA")
    call collct_deletevalue (rcollection,"DOMAIN")
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_2_collect
  
!<description>
  ! This is a "separated" poisson solver for solving a Poisson
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a collection structure for the communication:
  ! All subroutines add their generated information to the collection, so that
  ! the other subroutines can work with them.
  !
  ! The following tasks are performed by the subroutines:
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

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! A collection structure for our problem
    type(t_collection) :: rcollection
    
    ! Ok, let us start.
    ! We want to solve our Poisson problem on level...

    NLMAX = 7
    
    ! Initialise the collection.
    call collct_init (rcollection)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call pm2_initParamTriang (NLMAX,rcollection)
    call pm2_initDiscretisation (rcollection)
    call pm2_initMatVec (rcollection)
    call pm2_initDiscreteBC (rcollection)
    
    ! Implementation of boundary conditions
    call pm2_implementBC (rcollection)
    
    ! Solve the problem
    call pm2_solve (rcollection)
    
    ! Postprocessing
    call pm2_postprocessing (rcollection)
    
    ! Cleanup
    call pm2_doneMatVec (rcollection)
    call pm2_doneBC (rcollection)
    call pm2_doneDiscretisation (rcollection)
    call pm2_doneParamTriang (rcollection)
    
    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ("Remaining collection statistics:")
    call output_line ("--------------------------------")
    call output_lbrk ()
    call collct_printStatistics (rcollection)
    
    ! Finally release the collection.
    call collct_done (rcollection)
       
  end subroutine

end module

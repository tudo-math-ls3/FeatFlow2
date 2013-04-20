!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# The problem is a pure Neumann problem, solved with a direct solver.
!#
!# The reference function is
!#      u = 256 * ( x*(1-x) * y*(1-y) )^2
!# and is realised by separate callback functions in this module.
!# </purpose>
!##############################################################################

module poisson2d_method0_neumann

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
  use matrixmodification

  use poisson2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_NEUMANN_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    !    u(x,y) = 256*(x*(1-x)*y*(1-y))^2
    ! => f(x,y) = ... (a bit complicated)
    Dcoefficients (1,:,:) = &
        - 512.0_DP*(1-Dpoints(1,:,:))**2*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))**2 &
        +2048.0_DP*Dpoints(1,:,:)*(1-Dpoints(1,:,:))*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))**2 &
        - 512.0_DP*Dpoints(1,:,:)**2*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))**2 &
        - 512.0_DP*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))**2*(1-Dpoints(2,:,:))**2 & 
        +2048.0_DP*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))**2*Dpoints(2,:,:)*(1-Dpoints(2,:,:)) &
        - 512.0_DP*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))**2*Dpoints(2,:,:)**2

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReference_NEUMANN_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)                      :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>

    select case (cderivative)
    case (DER_FUNC)
      ! u(x,y) = 256*(x*(1-x)*y*(1-y))^2 - c
      ! with c = int( 256*(x*(1-x)*y*(1-y))^2 ) / |Omega| = 64/225 
      ! due to the mean value constraint.
      Dvalues (:,:) = 256.0_DP * ( Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) * &
                                  Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) )**2 &
                      - 64.0_DP/225.0_DP
                      
    case (DER_DERIV_X)
      !    u(x,y)   = 256*(x*(1-x)*y*(1-y))^2 - c
      ! => u_x(x,y) = ... (a bit complicated)
      Dvalues (:,:) = &
          512.0_DP*Dpoints(1,:,:)*(1-Dpoints(1,:,:))**2*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))**2 &
          -512.0_DP*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))**2
          
    case (DER_DERIV_Y)
      !    u(x,y)   = 256*(x*(1-x)*y*(1-y))^2 - c
      ! => u_y(x,y) = ... (a bit complicated)
      Dvalues (:,:) = &
          512.0_DP*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))**2*Dpoints(2,:,:)*(1-Dpoints(2,:,:))**2 &
          -512.0_DP*Dpoints(1,:,:)**2*(1-Dpoints(1,:,:))**2*Dpoints(2,:,:)**2*(1-Dpoints(2,:,:))
          
    case default
      ! Unknown. Set the result to 0.0.
      Dvalues = 0.0_DP
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_0_neumann
  
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
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A matrix, a RHS vector, a solution vector and a temporary vector. 
    ! The RHS vector accepts the RHS of the problem, the solution vector
    ! accepts the solution. All are block vectors with only one block.
    type(t_matrixBlock) :: rmatSystem
    type(t_vectorBlock) :: rvecSol,rvecRhs,rvecTmp
    type(t_matrixScalar) :: rmassMatrixLumped

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 7
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, trim(spredir)//"/QUAD.tri", rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                 rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_Q1,rtriangulation, rboundary)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use.
    call spdiscr_createDefCubStructure(&  
        rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_GEN_AUTO_G4)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 1x1 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! At first, create a basic 1x1 block matrix based on the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatSystem)
    
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9,rmatSystem%RmatrixBlock(1,1))
    
    ! And now to the entries of the matrix. 
    ! We do a quick approach here applying the method from the stdoperators
    ! module.
    call stdop_assembleLaplaceMatrix(&
        rmatSystem%RmatrixBlock(1,1), .true., 1.0_DP, rcubatureInfo)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (rdiscretisation,rvecRhs,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecSol,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecTmp,.true.)

    ! Now assemble the right-hand-side vector
    call linf_buildSimpleVector(rvecRhs%RvectorBlock(1), rcubatureInfo, &
        coeff_RHS_NEUMANN_2D)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! There are no Dirichlet boundary conditions, this is a pure
    ! Neumann problem. To successfully process this type of problem,
    ! we modify the system as well as the RHS using a lumped mass matrix
    ! approach. Note that this does only works for Lagrangian based
    ! finite elements (i.e., DOFs corresponds to values in vertices)
    
    ! Step 1: Generate a diagonal lumped mass matrix based on the structure of the 
    !         Laplace matrix
    call lsyssc_duplicateMatrix (rmatSystem%RmatrixBlock(1,1),rmassMatrixLumped,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    
    call lsyssc_clearMatrix (rmassMatrixLumped)
    
    call stdop_assembleSimpleMatrix (rmassMatrixLumped,DER_FUNC,DER_FUNC,1.0_DP)
    call lsyssc_lumpMatrix (rmassMatrixLumped,LSYSSC_LUMP_DIAG,.true.)
    
    ! Step 2: Filter the RHS to fulfil "int(f)=0".
    ! This means, subtract the integral mean value: f:=f-int(f)/|Omega|
    call vecfil_rhsL1To0ByLmass (rvecRhs%RvectorBlock(1),rmassMatrixLumped)
    
    ! Step 3: Modify the matrix/rhs and impose the mean value condition.
    !
    ! We impose the diagonal of a lumped mass matrix as first row into the
    ! matrix. Furthermore, we modify the first entry of the RHS to be =0.
    ! This imposes the condition "int(u)=0" into the system.
    !
    ! a) Modification of the RHS, impose a zero at f_1.
    call vecfil_oneEntryZero (rvecRhs,1,1)
    
    ! b) Impose the diagonal of the lumped mass matrix as first row
    !    of the system matrix. Note that this changes the structure
    !    of the matrix! The first row is full afterwards.
    call mmod_replaceLineByLumpedMass (rmatSystem%RmatrixBlock(1,1),1,rmassMatrixLumped)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create an UMFPACK direct solver.
    call linsol_initUmfpack4 (p_rsolverNode)
    
    ! Attach the system matrix to the solver.
    call linsol_setMatrix(p_rsolverNode,rmatSystem)
    
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
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solve the system
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rvecSol,rvecRhs,rvecTmp)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! That is it, rvecSol now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//"/u2d_0_neumann.vtk")
    
    ! Add the solution to the UCD exporter
    call ucd_addVectorByVertex (rexport, "sol", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(1))
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! Calculate the error to the reference function.
    call pperr_scalar (PPERR_L2ERROR,derror,rvecSol%RvectorBlock(1),&
                       getReference_NEUMANN_2D, rcubatureInfo=rcubatureInfo)
    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    call pperr_scalar (PPERR_H1ERROR,derror,rvecSol%RvectorBlock(1),&
                       getReference_NEUMANN_2D, rcubatureInfo=rcubatureInfo)
    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Clean up
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecSol)
    call lsysbl_releaseVector (rvecRhs)
    call lsysbl_releaseMatrix (rmatSystem)
    call lsyssc_releaseMatrix (rmassMatrixLumped)

    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# optimal control problem with constant coefficients on a simple domain.
!# </purpose>
!#
!# Problem and system definition
!# -----------------------------
!# The system to solve is the following:
!#
!#  -Laplace(y) = f + u
!#  -Laplace(p) = u - z
!#            u = P(-1/alpha p)
!#
!# with:
!#   u = primal solution
!#   p = dual solution
!#   z = target solution
!#
!# The boundary conditions for the dual solution are =0 on the
!# Dirichlet boundary.
!# The projection operator
!#
!#   P(g) = P(g)g = P[a,b](g)g = max(a,min(b,g))
!#
!# with
!#
!#   P(h)g = a     , for all x with h(x) < a
!#         = g(x)  , for all x with a <= h(x) <= b
!#         = b     , for all x with h(x) > b
!#
!# for a function g=g(x) defines the projection to the interval [a,b].
!# In the moment where this operator get's active, the system
!# is going to be nonlinear!
!#
!# In the simple case b=-a=infinity, the projection operator is the
!# identity, so we get the system
!#
!#  -Laplace(y) + 1/alpha p = f
!#  -Laplace(p) - u = -z
!#
!# which corresponds to the system
!#
!#  (  A    1/alpha M ) ( y )= (  f )
!#  ( -M    A         ) ( p )  ( -z )
!#
!# A simple UMFPACK Gauss elimination can be used to solve this
!# system, as it is completely linear.
!#
!# In the moment where the projection is active, the whole
!# thing gets a little bit harder. In that moment, we have
!# to solve
!#
!#  (  A    -P(u) ) ( y ) = (  f )
!#  ( -M    A     ) ( p )   ( -z )
!#
!# To solve this system, we can use a Newton iteration.
!# This is written down as follows:
!#
!#  ( y_n+1 ) = ( y_n ) + (  A   -P'(u_n)(-1/alpha .) )^-1 [ (  f ) - (  A    -P(u_n)(-1/alpha .) ) ( y_n ) ]
!#  ( p_n+1 )   ( p_n )   ( -M    A                   )    [ ( -z )   ( -M    A                   ) ( p_n ) ]
!#
!# <=>
!#
!#  ( y_n+1 ) = ( y_n ) + (  A    P'(u_n)(1/alpha .) )^-1 [ (  f ) - (  A    P(u_n)(1/alpha .) ) ( y_n ) ]
!#  ( p_n+1 )   ( p_n )   ( -M    A                  )    [ ( -z )   ( -M    A                 ) ( p_n ) ]
!#
!# where u_n = -1/alpha p_n and P'(g) (and its corresponding matrix) is
!# the semismooth Newton operator the projection P(g). This operator is defined by
!#
!#   P'(h)g (x) = g(x)  , if a <= h(x) <= b
!#              = 0     , elsewhere
!#
!# or in another way of writing:
!#
!#   P'(h) = Id  , for all x with a <= h(x) <= b
!#         = 0   , elsewhere.
!#
!# For the corresponding matrix, this can be interpreted as
!#
!#   P'(h)g = M  , for all nodes (vertices) x where a <= h(x) <= b
!#          = 0  , elsewhere (-> insert a zero row to the mass matrix!)
!#
!# Iterating this Newton iteration to convergence gives the optimal control u.
!#
!# -----
!#
!# Current problem:
!#
!# 1.) bpointBC=TRUE:
!#
!#  min J(y,u) = |y-z|^2 + alpha/2 |u|
!#
!# with alpha=0.001 and
!#
!#    -Laplace(y) = u
!#              y = 0    in (0,0)
!#          dy/dn = 0    on boundary of [0,1]^2 \ (0,0)
!#       -1 <= u <= 1
!#              z = x1 + x2
!#
!# 2.) bpointBC=FALSE:
!#
!#  min J(y,u) = |y-z|^2 + alpha/2 |u|
!#
!# with alpha=0.001 and
!#
!#    -Laplace(y) = u
!#              y = 0    on boundary of [0,1]^2
!#        0 <= u <= 15
!#              z = x1 + x2
!#
!##############################################################################

module poisson2d_method0_simple

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use linearalgebra
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use pprocerror
  use genoutput
  use stdoperators
  use spdiscprojection
  use convection
  use analyticprojection
  use optcontrolconvection
  use matrixio
  use collection
  use filtersupport
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use linearsolver
  use matrixmodification
    
  use poisson2d_callback
  
  implicit none

contains

  ! -----------------------------------------------------

  subroutine massmatfilter (rmatrix, rvector, dalphaC, dmin, dmax)
      
  ! Filters a mass matrix. The lines in the matrix rmatrix corresponding
  ! to all entries in the (control-)vector violating the constraints
  ! of the problem.
  
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! Vector containing a dial solution lambda. Whereever -1/alpha*lambda
  ! violates the control constraints given by rmatrixComponents, the corresponding
  ! lines are set to 0.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! ALPHA regularisation parameter from the space-time matrix
  real(dp), intent(in) :: dalphaC
  
  ! minimum bound for the control
  real(dp), intent(in) :: dmin

  ! maximum bound for the control
  real(dp), intent(in) :: dmax
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata
    integer, dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du
    
    ! Get the vector data
    call lsyssc_getbase_double (rvector,p_Ddata)
    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rvector%NEQ))
    
    nviolate = 0
    do i=1,rvector%NEQ
      du = -p_Ddata(i)/dalphaC
      if ((du .le. dmin) .or. (du .ge. dmax)) then
        nviolate = nviolate + 1
        p_Idofs(nviolate) = i
      end if
    end do
    
    if (nviolate .gt. 0) then
      ! Filter the matrix
      call mmod_replaceLinesByZero (rmatrix,p_Idofs(1:nviolate))
    end if
    
    deallocate(p_Idofs)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine projectControlTimestep (rdualSolution,dumin,dumax)

!<description>
  ! Projects a dual solution vector u in such a way, that
  ! dumin <= u <= dumax holds.
!</description>

!<input>
  ! Minimum value for u
  real(DP), intent(in) :: dumin

  ! Maximum value for u
  real(DP), intent(in) :: dumax
!</input>

!<inputoutput>
  ! Vector to be restricted
  type(t_vectorScalar), intent(inout) :: rdualSolution
!</inputoutput>

!</subroutine>
 
    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i
    
    ! Get the vector array
    call lsyssc_getbase_double (rdualSolution,p_Ddata)
    
    ! Restrict the vector
    do i=1,rdualSolution%NEQ
      p_Ddata(i) = min(max(p_Ddata(i),dumin),dumax)
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_0_simple
  
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
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    type(t_spatialDiscretisation) :: rdiscrOutput
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixLaplace, rmatrixMass

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempRhsBlock
    type(t_vectorScalar) :: rvectorTmp,rvectorOutput

    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Data for the Newton iteration
    integer :: ite, nmaxiterations
    real(dp) :: dinitRes, dcurrentRes
    
    ! Relaxation parameters
    real(DP) :: dalpha
    logical :: bboundsActive,bpointBC
    real(dp) :: dmax,dmin
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata

    ! Ok, let's start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 6
    
    ! Newton iteration counter
    nmaxiterations = 10
    
    ! Relaxation parameter
    dalpha = 0.001_DP
    
    ! TRUE=BC only on (0,0). FALSE=BC on whole boundary.
    bpointBC = .true.
    
    if (bpointBC) then
      ! Bounds on the control
      bboundsActive = .true.
      dmin = -1.0_DP
      dmax = 1.0_DP
    else
      ! Bounds on the control
      bboundsActive = .true.
      dmin = 0.0_DP
      dmax = 15.0_DP
    end if
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! Set up a block discretisation structure for two components:
    ! Primal and dual solution vector.
    call spdiscr_initBlockDiscr (rdiscretisation,2,&
                                 rtriangulation, rboundary)
    
    ! Set up the blocks. Both are discretised with the same finite element.
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
        EL_Q1,CUB_G3X3,rtriangulation, rboundary)
    call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
        rdiscretisation%RspatialDiscr(2), .true.)
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixLaplace)
                                     
    ! Assemble the Laplace matrix.
    call stdop_assembleLaplaceMatrix (rmatrixLaplace)
    
    ! And a mass matrix; we need it for the coupling.
    ! Share the structure with the Laplace matrix.
    call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixMass,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call stdop_assembleSimpleMatrix (rmatrixMass,DER_FUNC,DER_FUNC)

    ! -----
    ! Boundary conditions
    !
    ! In our example, we have pure Dirichlet-0-BC on all boundary components.
    
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Edge 1 of boundary component 1 the domain.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    if (bpointBC) rboundaryRegion%dmaxParam = 0.0_DP
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValuesPrimal_2D)
           
    if (.not. bpointBC) then
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesPrimal_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesPrimal_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesPrimal_2D)
    end if

    ! ---
    ! The same for the dual variable.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    if (bpointBC) rboundaryRegion%dmaxParam = 0.0_DP
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValuesDual_2D)
    if (.not. bpointBC) then
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesDual_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesDual_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesDual_2D)
    end if

    ! -----
    ! Linear system: Basic structure.
    !
    ! Create a 2x2 matrix and 2-block vectors from the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrixBlock)
    call lsysbl_createVecBlockByDiscr (rdiscretisation,rrhsBlock)
    call lsysbl_createVecBlockByDiscr (rdiscretisation,rvectorBlock)
    call lsysbl_createVecBlockByDiscr (rdiscretisation,rtempRhsBlock)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rvectorTmp)
    
    ! Create the RHS vector f of the primal equation
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(1),&
                                  coeff_RHS_primal_2D)
                                  
    ! Create the RHS vector -z of the dual equation
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(2),&
                                  coeff_RHS_dual_2D)
  
    call vecfil_discreteBCrhs (rrhsBlock,rdiscreteBC)
    
    ! -----
    ! Newton iteration
    !
    ! Initialise the solution vector.
    call lsysbl_clearVector (rvectorBlock)
    
    call vecfil_discreteBCsol (rvectorBlock,rdiscreteBC)
    
    do ite = 1,nmaxIterations
    
      ! Create the nonliear defect
      !  (d1) =  (  f ) - (  A    -P(u)(-1/alpha .) ) ( y )
      !  (d2)    ( -z )   ( -M    A                 ) ( p )
      ! We do this manually...
      
      call lsysbl_copyVector (rrhsBlock,rtempRhsBlock)
      
      call lsyssc_scalarMatVec (rmatrixLaplace,&
          rvectorBlock%RvectorBlock(1),rtempRhsBlock%RvectorBlock(1),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rmatrixLaplace,&
          rvectorBlock%RvectorBlock(2),rtempRhsBlock%RvectorBlock(2),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rmatrixMass,&
          rvectorBlock%RvectorBlock(1),rtempRhsBlock%RvectorBlock(2),1.0_DP,1.0_DP)
          
      ! No constraints: -P(u_n) = 1/alpha p_n -> multiply p_n by weighted Mass matrix.
      ! Constraints active: Create the projection and multiply by the weighted mass matrix.
      call lsyssc_copyVector (rvectorBlock%RvectorBlock(2),rvectorTmp)
      call lsyssc_scaleVector (rvectorTmp,-1.0_DP/dalpha)
      if (bboundsActive) then
        call projectControlTimestep (rvectorTmp,dmin,dmax)
      end if
      call lsyssc_scalarMatVec (rmatrixMass,&
          rvectorTmp,rtempRhsBlock%RvectorBlock(1),1.0_DP,1.0_DP)
      
      ! Prepare the preconditioner matrix (Newton).
      call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(2,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      rmatrixBlock%RmatrixBlock(2,1)%dscaleFactor = -1.0_DP
      
      ! Copy also the mass matrix for the 4th block. In the case of no
      ! control constraints, that's it!
      call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(1,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      rmatrixBlock%RmatrixBlock(1,2)%dscaleFactor = -(-1.0_DP/dalpha)

      if (bboundsActive) then
        ! Filter the mass matrix. Set those rows to zero where the DOF's are out
        ! of bounds. The result is the derivative of the projection operator...
        call massmatfilter (rmatrixBlock%RmatrixBlock(1,2), &
            rvectorBlock%RvectorBlock(2), dalpha, dmin, dmax)
      end if
    
      ! Include boundary conditions
      call vecfil_discreteBCdef (rtempRhsBlock,rdiscreteBC)
      
      ! Check for convergence
      if (ite .eq. 1) then
        dinitRes = lsysbl_vectorNorm(rtempRhsBlock,LINALG_NORML2)
        call output_line ('Iteration: '//TRIM(sys_siL(ite-1,10))//&
            ', Initial defect: '//sys_sdEL(dinitRes,10))
      else
        dcurrentRes = lsysbl_vectorNorm(rtempRhsBlock,LINALG_NORML2)
        call output_line ('Iteration: '//TRIM(sys_siL(ite-1,10))//&
            ', Current defect: '//sys_sdEL(dcurrentRes,10))
        if (dcurrentRes .lt. dinitRes * 1E-10_DP) exit
      end if
      
      call matfil_discreteBC (rmatrixBlock,rdiscreteBC)
      
      ! Prepare an UMFPACK solver for the system
      call linsol_initUMFPACK4(p_rsolverNode)
    
      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2
      
      ! Attach the system matrix to the solver.
      ! First create an array with the matrix data (on all levels, but we
      ! only have one level here), then call the initialisation
      ! routine to attach all these matrices.
      ! Remark: Don't make a call like
      !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
      ! This doesn't work on all compilers, since the compiler would have
      ! to create a temp array on the stack - which does not always work!
      Rmatrices = (/rmatrixBlock/)
      call linsol_setMatrices(p_RsolverNode,Rmatrices)
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Preconditioning of the defect by the inverse of the Newton matrix
      call linsol_precondDefect (p_rsolverNode,rtempRhsBlock)
      
      ! Release solver data and structure
      call linsol_doneData (p_rsolverNode)
      call linsol_doneStructure (p_rsolverNode)
      
      ! Release the solver node and all subnodes attached to it (if at all):
      call linsol_releaseSolver (p_rsolverNode)
      
      ! Sum up the correction to the current solution.
      call lsysbl_vectorLinearComb (rtempRhsBlock,rvectorBlock,1.0_DP,1.0_DP)
      
      ! Plug in the BC's to the solution
      call vecfil_discreteBCsol (rvectorBlock,rdiscreteBC)
      
    end do
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/ups2d_0_simple.gmv')
                       
    call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(1),&
      rvectorOutput,rdiscrOutput)
    
    call lsyssc_getbase_double (rvectorOutput,p_Ddata)
    
    call ucd_addVariableVertexBased (rexport,'solprimal',UCD_VAR_STANDARD, p_Ddata)
    call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(2),&
        rvectorOutput,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'soldual',UCD_VAR_STANDARD, p_Ddata)
    
    call lsyssc_copyVector (rvectorBlock%RvectorBlock(2),rvectorTmp)
    call lsyssc_scaleVector (rvectorTmp,-1.0_DP/dalpha)
    if (bboundsActive) then
      call projectControlTimestep (rvectorTmp,dmin,dmax)
    end if
    call spdp_stdProjectionToP1Q1Scalar (rvectorTmp,&
        rvectorOutput,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'control',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )

    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    
    ! Release the block matrix/vectors
    call lsyssc_releaseVector (rvectorOutput)
    call lsyssc_releaseVector (rvectorTmp)
    call lsysbl_releaseVector (rtempRhsBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseMatrix (rmatrixBlock)

    call lsyssc_releaseMatrix (rmatrixMass)
    call lsyssc_releaseMatrix (rmatrixLaplace)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseDiscr(rdiscrOutput)
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

end module

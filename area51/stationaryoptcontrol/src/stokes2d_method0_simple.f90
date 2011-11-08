!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Stokes
!# optimal control problem with constant coefficients on a simple domain.
!# </purpose>
!#
!# Problem and system definition
!# -----------------------------
!# The system to solve is the following:
!#
!#        -Laplace(y) + grad(p) = f + u
!#                       -div y = 0
!#  -Laplace(lambda) + grad(xi) = u - z
!#                  -div lambda = 0
!#                            u = P(-1/alpha l)
!#
!# with:
!#   u/p = primal solution
!#   l/xi = dual solution
!#   z = target velocity
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
!#  -Laplace(y) + 1/alpha l = f
!#  -Laplace(p) - u = -z
!#
!# which corresponds to the system
!#
!#  (  A    1/alpha M ) ( y )= (  f )
!#  ( -M    A         ) ( l )  ( -z )
!#
!# A simple UMFPACK Gauss elimination can be used to solve this
!# system, as it is completely linear.
!#
!# In the moment where the projection is active, the whole
!# thing gets a little bit harder. In that moment, we have
!# to solve
!#
!#  (  A    -P(u) ) ( y ) = (  f )
!#  ( -M    A     ) ( l )   ( -z )
!#
!# To solve this system, we can use a Newton iteration.
!# This is written down as follows:
!#
!#  ( y_n+1 ) = ( y_n ) + (  A   -P'(u_n)(-1/alpha .) )^-1 [ (  f ) - (  A    -P(u_n)(-1/alpha .) ) ( y_n ) ]
!#  ( l_n+1 )   ( l_n )   ( -M    A                   )    [ ( -z )   ( -M    A                   ) ( l_n ) ]
!#
!# <=>
!#
!#  ( y_n+1 ) = ( l_n ) + (  A    P'(u_n)(1/alpha .) )^-1 [ (  f ) - (  A    P(u_n)(1/alpha .) ) ( y_n ) ]
!#  ( l_n+1 )   ( p_n )   ( -M    A                  )    [ ( -z )   ( -M    A                 ) ( l_n ) ]
!#
!# where u_n = -1/alpha l_n and P'(g) (and its corresponding matrix) is
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
!##############################################################################

module stokes2d_method0_simple

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
  use matrixio
    
  use stokes2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine stokes2d_0_simple
  
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
    type(t_bilinearForm) :: rbilinearForm
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixLaplace, rmatrixMass
    type(t_matrixScalar) :: rmatrixB1, rmatrixB2

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock), target :: rvectorBlock
    type(t_vectorBlock) :: rrhsBlock,rtempRhsBlock
    type(t_vectorScalar) :: rvectorTmp,rvectorOutputX,rvectorOutputY

    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rprecond
    TYPE(t_filterChain), DIMENSION(2), TARGET :: RfilterChain

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Data for the Newton iteration
    integer :: ite, nmaxiterations, idiscretisation
    real(dp) :: dinitRes, dcurrentRes
    
    ! Relaxation parameters
    real(DP) :: dalpha
    logical :: bboundsActive, bcontrolactive, bdualcoupledtoprimal
    logical :: bemulateTimestep, bnewton,bexactderiv
    real(dp) :: dmin1,dmax1,dmin2,dmax2, dnu,dt
    
    type(t_collection) :: rcollection
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataX,p_DdataY,p_DdataRhs,p_DdataRhsTemp

    ! Ok, let's start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 2
    
    ! Newton iteration counter
    nmaxiterations = 1000
    
    ! Relaxation parameter
    dalpha = 0.01_DP
    
    ! Control active or not.
    bcontrolactive = .false.
    bdualcoupledtoprimal = .false.
    
    ! Bounds on the control
    bboundsActive = .false.
    dmin1 = -1E99_DP
    dmax1 = 0.8_DP !8.5_DP
    dmin2 = -1E99_DP
    dmax2 = 1E99_DP
    
    ! Viscosity constant
    dnu = 1.0_DP
    
    ! Discretisation. 1=EM30, 2=Q2
    idiscretisation = 1
    
    ! TRUE emulates a timestep with implicit Euler, timestep length dt,
    ! starting from a zero solution
    bemulateTimestep = .false.
    
    ! Timestep length
    dt = 1.0_DP
    
    ! TRUE: Use Newton iteration
    bnewton = .true.
    
    ! TRUE: Use exact derivative of the semismooth operator
    ! (mass matrix set up with an appropriate coefficient) instead
    ! of zero rows in the mass matrix.
    bexactderiv = .true.
    
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
    call spdiscr_initBlockDiscr (rdiscretisation,6,&
                                 rtriangulation, rboundary)
    
    ! Set up the blocks. Both are discretised with the same finite element.
    if (idiscretisation .eq. 1) then
    
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          EL_EM30,CUB_G3X3,rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(4), .true.)
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(5), .true.)
          
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
          EL_Q0,CUB_G2X2,rdiscretisation%RspatialDiscr(3))
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(3), &
          rdiscretisation%RspatialDiscr(6), .true.)
          
    else if (idiscretisation .eq. 2) then
    
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          EL_Q2,CUB_G4X4,rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(4), .true.)
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(5), .true.)
          
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
          EL_QP1,CUB_G2X2,rdiscretisation%RspatialDiscr(3))
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(3), &
          rdiscretisation%RspatialDiscr(6), .true.)
    end if
                 
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

    ! Create B-matrices
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(3),&
        LSYSSC_MATRIX9,rmatrixB1,rdiscretisation%RspatialDiscr(1))
    call lsyssc_duplicateMatrix (rmatrixB1,rmatrixB2,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
        
    call stdop_assembleSimpleMatrix (rmatrixB1,DER_FUNC,DER_DERIV_X,-1.0_DP,.true.)
    call stdop_assembleSimpleMatrix (rmatrixB2,DER_FUNC,DER_DERIV_Y,-1.0_DP,.true.)

    ! -----
    ! Boundary conditions
    !
    ! In our example, we have pure Dirichlet-0-BC on all boundary components.
    
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Edge 1 of boundary component 1 the domain.
    rcollection%IquickAccess(1) = 1
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
                             
    ! Edge 3 of boundary component 1.
    rcollection%IquickAccess(1) = 2
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
    

    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    ! ---
    ! The same for the dual variable.
    rcollection%IquickAccess(1) = 4
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
                                       
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
                             
    ! Edge 3 of boundary component 1.
    rcollection%IquickAccess(1) = 5
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D,rcollection)

    ! -----
    ! Linear system: Basic structure.
    !
    ! Create a 2x2 matrix and 2-block vectors from the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrixBlock)
    call lsysbl_createVecBlockByDiscr (rdiscretisation,rrhsBlock)
    call lsysbl_createVecBlockByDiscr (rdiscretisation,rvectorBlock)
    call lsysbl_createVecBlockByDiscr (rdiscretisation,rtempRhsBlock)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rvectorTmp)
    
    call lsysbl_getbase_double (rtempRhsBlock,p_DdataRhsTemp)
    call lsysbl_getbase_double (rRhsBlock,p_DdataRhs)
    
    ! Create the RHS vector f of the primal equation
    rcollection%DquickAccess(1) = dalpha
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    rcollection%IquickAccess(1) = 1
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(1),&
                                  coeff_RHS_primal_2D,rcollection)
    
    rcollection%IquickAccess(1) = 2
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(2),&
                                  coeff_RHS_primal_2D,rcollection)
    call lsyssc_clearVector (rrhsBlock%RvectorBlock(3))
                                  
    ! Create the RHS vector -z of the dual equation
    rcollection%IquickAccess(1) = 1
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(4),&
                                  coeff_RHS_dual_2D,rcollection)
    rcollection%IquickAccess(1) = 2
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(5),&
                                  coeff_RHS_dual_2D,rcollection)
    call lsyssc_clearVector (rrhsBlock%RvectorBlock(6))
  
    ! Implement BC's
    call vecfil_discreteBCrhs (rrhsBlock,rdiscreteBC)
    
    ! -----
    ! Newton iteration
    !
    ! Initialise the solution vector.
    call lsysbl_clearVector (rvectorBlock)
    
    call vecfil_discreteBCsol (rvectorBlock,rdiscreteBC)
    
    do ite = 1,nmaxIterations+1
    
!      if (ite .eq. 2) then
!        bcontrolactive = .false.
!        rmatrixBlock%RmatrixBlock(1,4)%dscaleFactor = 0.0_DP
!        rmatrixBlock%RmatrixBlock(2,5)%dscaleFactor = 0.0_DP
!      end if
    
      ! Create the nonliear defect
      !  (d1) =  (  f ) - (  A    -P(u)(-1/alpha .) ) ( y )
      !  (d2)    ( -z )   ( -M    A                 ) ( p )
      ! We do this manually...
      
      call lsysbl_copyVector (rrhsBlock,rtempRhsBlock)
      
      call lsyssc_scalarMatVec (rmatrixLaplace,&
          rvectorBlock%RvectorBlock(1),rtempRhsBlock%RvectorBlock(1),-dnu,1.0_DP)
      call lsyssc_scalarMatVec (rmatrixLaplace,&
          rvectorBlock%RvectorBlock(2),rtempRhsBlock%RvectorBlock(2),-dnu,1.0_DP)
    
      if (bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorBlock%RvectorBlock(1),rtempRhsBlock%RvectorBlock(1),-1.0_DP/dt,1.0_DP)
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorBlock%RvectorBlock(2),rtempRhsBlock%RvectorBlock(2),-1.0_DP/dt,1.0_DP)
      end if

      call lsyssc_scalarMatVec (rmatrixB1,&
          rvectorBlock%RvectorBlock(3),rtempRhsBlock%RvectorBlock(1),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rmatrixB2,&
          rvectorBlock%RvectorBlock(3),rtempRhsBlock%RvectorBlock(2),-1.0_DP,1.0_DP)

      call lsyssc_scalarMatVec (rmatrixB1,&
          rvectorBlock%RvectorBlock(1),rtempRhsBlock%RvectorBlock(3),-1.0_DP,1.0_DP,&
          .true.)
      call lsyssc_scalarMatVec (rmatrixB2,&
          rvectorBlock%RvectorBlock(2),rtempRhsBlock%RvectorBlock(3),-1.0_DP,1.0_DP,&
          .true.)
      
      call lsyssc_scalarMatVec (rmatrixLaplace,&
          rvectorBlock%RvectorBlock(4),rtempRhsBlock%RvectorBlock(4),-dnu,1.0_DP)
      call lsyssc_scalarMatVec (rmatrixLaplace,&
          rvectorBlock%RvectorBlock(5),rtempRhsBlock%RvectorBlock(5),-dnu,1.0_DP)

      if (bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorBlock%RvectorBlock(4),rtempRhsBlock%RvectorBlock(4),-1.0_DP/dt,1.0_DP)
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorBlock%RvectorBlock(5),rtempRhsBlock%RvectorBlock(5),-1.0_DP/dt,1.0_DP)
      end if

      call lsyssc_scalarMatVec (rmatrixB1,&
          rvectorBlock%RvectorBlock(6),rtempRhsBlock%RvectorBlock(4),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rmatrixB2,&
          rvectorBlock%RvectorBlock(6),rtempRhsBlock%RvectorBlock(5),-1.0_DP,1.0_DP)

      call lsyssc_scalarMatVec (rmatrixB1,&
          rvectorBlock%RvectorBlock(4),rtempRhsBlock%RvectorBlock(6),-1.0_DP,1.0_DP,&
          .true.)
      call lsyssc_scalarMatVec (rmatrixB2,&
          rvectorBlock%RvectorBlock(5),rtempRhsBlock%RvectorBlock(6),-1.0_DP,1.0_DP,&
          .true.)
      
      if (bdualcoupledtoprimal) then
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorBlock%RvectorBlock(1),rtempRhsBlock%RvectorBlock(4),1.0_DP,1.0_DP)
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorBlock%RvectorBlock(2),rtempRhsBlock%RvectorBlock(5),1.0_DP,1.0_DP)
      end if
      
      if (bcontrolactive) then
        ! No constraints: -P(u_n) = 1/alpha p_n -> multiply p_n by weighted Mass matrix.
        ! Constraints active: Create the projection and multiply by the weighted mass matrix.
        call lsyssc_copyVector (rvectorBlock%RvectorBlock(4),rvectorTmp)
        call lsyssc_scaleVector (rvectorTmp,-1.0_DP/dalpha)
        if (bboundsActive) then
          call projectControlTimestep (rvectorTmp,dmin1,dmax1)
        end if
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorTmp,rtempRhsBlock%RvectorBlock(1),1.0_DP,1.0_DP)

        call lsyssc_copyVector (rvectorBlock%RvectorBlock(5),rvectorTmp)
        call lsyssc_scaleVector (rvectorTmp,-1.0_DP/dalpha)
        if (bboundsActive) then
          call projectControlTimestep (rvectorTmp,dmin2,dmax2)
        end if
        call lsyssc_scalarMatVec (rmatrixMass,&
            rvectorTmp,rtempRhsBlock%RvectorBlock(2),1.0_DP,1.0_DP)
            
      end if

      ! ---
      ! Prepare the preconditioner matrix (Newton).
      call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(1,1),dnu)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(2,2),dnu)

      if (bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_matrixLinearComb (rmatrixMass,1.0_DP/dt,&
            rmatrixBlock%RmatrixBlock(1,1),1.0_DP,rmatrixBlock%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)
        call lsyssc_matrixLinearComb (rmatrixMass,1.0_DP/dt,&
            rmatrixBlock%RmatrixBlock(2,2),1.0_DP,rmatrixBlock%RmatrixBlock(2,2),&
            .false.,.false.,.true.,.true.)
      end if

      call lsyssc_duplicateMatrix (rmatrixB1,rmatrixBlock%RmatrixBlock(1,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rmatrixB2,rmatrixBlock%RmatrixBlock(2,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_transposeMatrix (rmatrixB1,rmatrixBlock%RmatrixBlock(3,1),LSYSSC_TR_ALL)
      call lsyssc_transposeMatrix (rmatrixB2,rmatrixBlock%RmatrixBlock(3,2),LSYSSC_TR_ALL)

      call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,4),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,5),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(4,4),dnu)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(5,5),dnu)

      if (bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_matrixLinearComb (rmatrixMass,1.0_DP/dt,&
            rmatrixBlock%RmatrixBlock(4,4),1.0_DP,rmatrixBlock%RmatrixBlock(4,4),&
            .false.,.false.,.true.,.true.)
        call lsyssc_matrixLinearComb (rmatrixMass,1.0_DP/dt,&
            rmatrixBlock%RmatrixBlock(5,5),1.0_DP,rmatrixBlock%RmatrixBlock(5,5),&
            .false.,.false.,.true.,.true.)
      end if

      call lsyssc_duplicateMatrix (rmatrixB1,rmatrixBlock%RmatrixBlock(4,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rmatrixB2,rmatrixBlock%RmatrixBlock(5,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_transposeMatrix (rmatrixB1,rmatrixBlock%RmatrixBlock(6,4),LSYSSC_TR_ALL)
      call lsyssc_transposeMatrix (rmatrixB2,rmatrixBlock%RmatrixBlock(6,5),LSYSSC_TR_ALL)

      if (bdualcoupledtoprimal) then
        call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(4,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(5,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        rmatrixBlock%RmatrixBlock(4,1)%dscaleFactor = -1.0_DP
        rmatrixBlock%RmatrixBlock(5,2)%dscaleFactor = -1.0_DP
      end if

      if (bcontrolactive) then
        if (.not. bexactderiv) then
          ! Copy also the mass matrix for the 4th block. In the case of no
          ! control constraints, that's it!
          call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(1,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(2,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          rmatrixBlock%RmatrixBlock(1,4)%dscaleFactor = -(-1.0_DP/dalpha)
          rmatrixBlock%RmatrixBlock(2,5)%dscaleFactor = -(-1.0_DP/dalpha)

          if (bboundsActive .and. bnewton) then
            ! Filter the mass matrix. Set those rows to zero where the DOF's are out
            ! of bounds. The result is the derivative of the projection operator...
            call massmatfilter (rmatrixBlock%RmatrixBlock(1,4), &
                rvectorBlock%RvectorBlock(4), dalpha, dmin1, dmax1)
            call massmatfilter (rmatrixBlock%RmatrixBlock(2,5), &
                rvectorBlock%RvectorBlock(5), dalpha, dmin2, dmax2)
          end if
        
        else

          call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(1,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
          call lsyssc_duplicateMatrix (rmatrixMass,rmatrixBlock%RmatrixBlock(2,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
          rbilinearForm%itermCount = 1
          rbilinearForm%Idescriptors(1,1) = DER_FUNC
          rbilinearForm%Idescriptors(2,1) = DER_FUNC
          rbilinearForm%ballCoeffConstant = .false.
          rbilinearForm%BconstantCoeff(1) = .false.
          rcollection%p_rvectorQuickAccess1 => rvectorBlock
          rcollection%DquickAccess(3) = dalpha
          rcollection%IquickAccess(1) = 1
          if (bboundsActive .and. bnewton) then
            rcollection%DquickAccess(1) = dmin1
            rcollection%DquickAccess(2) = dmax1
          else
            rcollection%DquickAccess(1) = -SYS_MAXREAL_DP
            rcollection%DquickAccess(2) = SYS_MAXREAL_DP
          end if
          call bilf_buildMatrixScalar (rbilinearForm,.true.,rmatrixBlock%RmatrixBlock(1,4),&
                                      coeff_ProjMass,rcollection)

          rcollection%IquickAccess(1) = 2
          if (bboundsActive .and. bnewton) then
            rcollection%DquickAccess(1) = dmin2
            rcollection%DquickAccess(2) = dmax2
          else
            rcollection%DquickAccess(1) = -SYS_MAXREAL_DP
            rcollection%DquickAccess(2) = SYS_MAXREAL_DP
          end if
          call bilf_buildMatrixScalar (rbilinearForm,.true.,rmatrixBlock%RmatrixBlock(2,5),&
                                      coeff_ProjMass,rcollection)
                                      
        end if
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
        if (dcurrentRes .lt. 1E-14_DP) exit
        if (dcurrentRes .lt. dinitRes * 1E-13_DP) exit
        if (ite .gt. nmaxIterations) exit
      end if
      
      call matfil_discreteBC (rmatrixBlock,rdiscreteBC)
      
      RfilterChain(1)%ifilterType = FILTER_SMALLL1TO0
      RfilterChain(1)%ismallL1to0component = 3
      RfilterChain(2)%ifilterType = FILTER_SMALLL1TO0
      RfilterChain(2)%ismallL1to0component = 6
    
      ! Prepare an UMFPACK solver for the system
      !
      !call linsol_initUMFPACK4(p_rsolverNode)
      call linsol_initVANKA(p_rprecond,0.7_DP,LINSOL_VANKA_GENERAL)
      CALL linsol_initDefCorr (p_rsolverNode,p_rprecond,RfilterChain)
      p_rsolverNode%ioutputLevel = 2
      p_rsolverNode%nmaxIterations = 1000
      p_rsolverNode%depsRel = 1E-13_DP
      
      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2
      
      !call matio_writeBlockMatrixHR (rmatrixBlock, 'matrix',&
      !    .true., 0, "matrix.txt","(E20.10)")
      
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
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/ust2d_0_simple.vtk')
                       
    call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(1),&
      rvectorOutputX,rdiscrOutput)
    call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(2),&
      rvectorOutputY,rdiscrOutput)
    
    call lsyssc_getbase_double (rvectorOutputX,p_DdataX)
    call lsyssc_getbase_double (rvectorOutputY,p_DdataY)
    
    call ucd_addVarVertBasedVec (rexport, 'primalvel', p_DdataX, p_DdataY)

    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'primalP',UCD_VAR_STANDARD, p_Ddata)
    
    call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(4),&
        rvectorOutputX,rdiscrOutput)
    call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(5),&
        rvectorOutputY,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'dualvelX',UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableVertexBased (rexport,'dualvelY',UCD_VAR_STANDARD, p_DdataY)

    call ucd_addVarVertBasedVec (rexport, 'dualvel', p_DdataX, p_DdataY)

    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(6),p_Ddata)
    call ucd_addVariableElementBased (rexport,'dualP',UCD_VAR_STANDARD, p_Ddata)
    
    call lsyssc_copyVector (rvectorBlock%RvectorBlock(4),rvectorTmp)
    call lsyssc_scaleVector (rvectorTmp,-1.0_DP/dalpha)
    if (bboundsActive) then
      call projectControlTimestep (rvectorTmp,dmin1,dmax1)
    end if
    call spdp_stdProjectionToP1Q1Scalar (rvectorTmp,rvectorOutputX,rdiscrOutput)
    
    call lsyssc_copyVector (rvectorBlock%RvectorBlock(5),rvectorTmp)
    call lsyssc_scaleVector (rvectorTmp,-1.0_DP/dalpha)
    if (bboundsActive) then
      call projectControlTimestep (rvectorTmp,dmin2,dmax2)
    end if
    call spdp_stdProjectionToP1Q1Scalar (rvectorTmp,rvectorOutputY,rdiscrOutput)
    
    call ucd_addVariableVertexBased (rexport,'controlX',UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableVertexBased (rexport,'controlY',UCD_VAR_STANDARD, p_DdataY)
    
    ! residual
    call spdp_stdProjectionToP1Q1Scalar (rtempRhsBlock%RvectorBlock(1),&
        rvectorOutputX,rdiscrOutput)
    call spdp_stdProjectionToP1Q1Scalar (rtempRhsBlock%RvectorBlock(2),&
        rvectorOutputY,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'respx',UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableVertexBased (rexport,'respy',UCD_VAR_STANDARD, p_DdataY)

    call spdp_stdProjectionToP1Q1Scalar (rtempRhsBlock%RvectorBlock(4),&
        rvectorOutputX,rdiscrOutput)
    call spdp_stdProjectionToP1Q1Scalar (rtempRhsBlock%RvectorBlock(5),&
        rvectorOutputY,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'resdx',UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableVertexBased (rexport,'resdy',UCD_VAR_STANDARD, p_DdataY)
    
    ! rhs
    call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(1),&
        rvectorOutputX,rdiscrOutput)
    call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(1),&
        rvectorOutputY,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'rhspx',UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableVertexBased (rexport,'rhspy',UCD_VAR_STANDARD, p_DdataY)

    call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(4),&
        rvectorOutputX,rdiscrOutput)
    call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(5),&
        rvectorOutputY,rdiscrOutput)
    call ucd_addVariableVertexBased (rexport,'rhsdx',UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableVertexBased (rexport,'rhsdy',UCD_VAR_STANDARD, p_DdataY)

    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    
    ! Release the block matrix/vectors
    call lsyssc_releaseVector (rvectorOutputY)
    call lsyssc_releaseVector (rvectorOutputX)
    call lsyssc_releaseVector (rvectorTmp)
    call lsysbl_releaseVector (rtempRhsBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseMatrix (rmatrixBlock)

    call lsyssc_releaseMatrix (rmatrixB2)
    call lsyssc_releaseMatrix (rmatrixB1)
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
  
end module

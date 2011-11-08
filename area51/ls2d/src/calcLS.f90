
module CalcLS

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use pprocerror
  use genoutput

  use ccbasic
  use timestepping
  use trilinearformevaluation
  use convection
  use cccallback
  use levelset
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC

    type(t_blockDiscretisation) :: rdiscretisationStabil
  
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine calcls_mg(rproblem,rvector,&
                 rvectorVelocity,rtimestepping)


    type(t_problem), intent(in) :: rproblem
    type(t_vectorBlock),intent(in) :: rvectorVelocity
    TYPE(t_explicitTimeStepping),intent(in)  :: rtimestepping

    type(t_vectorBlock),intent(inout) :: rvector
    

    type(t_blockDiscretisation),pointer :: p_rdiscretisation

    

    ! An object for saving the domain:
    type(t_level),dimension(:),allocatable,target :: rlevels
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_trilinearForm) :: rtriform
    
    
    type(t_vectorBlock) :: rrhsBlock,rtempBlock

    ! A variable that is used to specify a region on the boundary.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), allocatable :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfo
    
    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    
    real(DP), dimension(:), pointer :: p_Ddata

    type(t_jumpStabilisation) :: rjumpstabil
    type(t_vectorblock) :: rvectortemp
    
    ! Some temporary variables
    integer :: i

    ! Ok, let's start.
    !
    ! We want to solve our Poisson problem on level...
    NLMIN = rproblem%NLMIN
    NLMAX = rproblem%NLMAX
    
    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))
    
    
   call lsysbl_copyvector (rvectorvelocity,rvectortemp)
  
    do i = NLMAX, NLMIN,-1


      
      p_rdiscretisation => rproblem%rlevelinfo(i)%rdiscretisationLS

      call spdiscr_deriveBlockDiscr (p_rdiscretisation,&
               rlevels(i)%rdiscretisationStabil)

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (&
          p_rdiscretisation,Rlevels(i)%rmatrix)

      CALL lsyssc_duplicateMatrix (rproblem%rlevelinfo(i)%rmatrixLS,&
              Rlevels(i)%rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_copy)

     
      call lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)

      rtriform%itermcount=1
      rtriform%BconstantCoeff= .true.
      rtriform%ballcoeffconstant = .true.
      rtriform%Dcoefficients(1)= rtimestepping%dtstep
      rtriform%Idescriptors(1,1)=DER_FUNC
      rtriform%Idescriptors(2,1)=DER_DERIV_X
      rtriform%Idescriptors(3,1)=DER_FUNC
         
     CALL trilf_buildMatrixScalar (rtriform,.FALSE.,Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                                   rvectortemp%rvectorBlock(1))
       
     rtriform%Idescriptors(2,1)=DER_DERIV_Y
        
     CALL trilf_buildMatrixScalar (rtriform,.FALSE.,Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                              rvectortemp%rvectorBlock(2))


     rjumpstabil%dnu = 0.0_DP
     rjumpstabil%dgamma = 0.01_DP
     rjumpstabil%dgammastar = 0.01_DP
     rjumpstabil%dtheta = rtimestepping%dtstep
    
     !call conv_JumpStabilisation2d ( &
     !  rjumpstabil, conv_modmatrix, rlevels(i)%rmatrix%rmatrixBlock(1,1),&
     !    rdiscretisation=rlevels(i)%rdiscretisationStabil%RspatialDiscr(1))

    IF(i.GT.rproblem%NLMIN) &
        CALL prjF2C(rvectorTemp,rproblem%rlevelinfo(i-1)%rdiscretisationLS)
    
    end do
      
    
    !call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhsBlock, .true.)
    call lsysbl_createVecBlockIndirect (rvector, rrhsBlock, .true.)

    CALL lsyssc_scalarMatVec (rproblem%rlevelinfo(nlmax)%rmatrixLS, &
                                   rvector%rvectorBlock(1),rrhsBlock%rvectorBlock(1),1.0_DP, rtimestepping%dtstep)
    
    do i = NLMIN, NLMAX
    
      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

     
     CALL boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
     CALL bcasm_newDirichletBConRealBD (rproblem%rlevelinfo(i)%rdiscretisationLS,1,&
                                       rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                      getBoundaryValues_trace)
    
      
      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
  
      ! Also implement the boundary conditions into the matrix.
      call matfil_discreteBC (Rlevels(i)%rmatrix)
      
    end do

    ! Our right-hand-side also needs to know the boundary conditions.
    rrhsBlock%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
    
   
    !call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
    
 
    call vecfil_discreteBCrhs (rrhsBlock)
    call vecfil_discreteBCsol (rvector)
    
 
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    call mlprj_initProjectionMat (rprojection,Rlevels(NLMAX)%rmatrix)

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    call linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Set up a coarse grid solver.
    call linsol_initUMFPACK4 (p_rcoarsegridSolver)
    
    ! Add the coarse grid level.
    call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                  null(), null(), p_rcoarseGridSolver)

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      ! Create a Jacobi smoother
      !CALL linsol_initJacobi(p_rsmoother)
      
      ! Create an ILU(0) smoother
      call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                    p_rsmoother, p_rsmoother, null())
      
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(NLMIN:NLMAX))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release Rlevels(i)%rmatrix!
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,rvector,rrhsBlock,rtempBlock)
    
 
    
    ! Calculate the error to the reference function.
    !call pperr_scalar (rvector%RvectorBlock(1),PPERR_L2ERROR,derror,&
    !                   getReferenceFunction_2D)
    !call output_line ('L2-error: ' // sys_sdEL(derror,10) )

    !call pperr_scalar (rvector%RvectorBlock(1),PPERR_H1ERROR,derror,&
    !                   getReferenceFunction_2D)
    !call output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
   
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the interlevel projection structure
    call mlprj_doneProjection (rprojection)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rrhsBlock)

    do i = NLMAX, NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
      call spdiscr_releaseBlockDiscr (Rlevels(i)%rdiscretisationstabil)
    end do

    ! Release our discrete version of the boundary conditions
    do i = NLMAX, NLMIN, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    !do i = NLMAX, NLMIN, -1
    !  call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    !end do
    
    ! Release the triangulation.
    !do i = NLMAX, NLMIN, -1
    !  call tria_done (Rlevels(i)%rtriangulation)
    !end do
    
    deallocate(Rlevels)
    
    call lsysbl_releasevector(rvectortemp)

  end subroutine
  
! *********************************************************************
!  <subroutine>
  subroutine CalcLS_simple(rproblem,rvector,&
             rvectorvelocity_new,rvectorvelocity_old,rtimestepping)
  
  type(t_problem), intent(in),target :: rproblem
  type(t_vectorBlock), intent(in), target :: rvectorvelocity_new,rvectorvelocity_old
  type(t_explicitTimeStepping),intent(in)  :: rtimestepping
  type(t_vectorBlock), intent(inout) :: rvector
     
  type(t_vectorBlock),target :: rrhsblock,rtempBlock
  type(t_matrixBlock) :: rmatrix,rmatrix2
! type(t_convStreamlineDiffusion) :: rstreamline
! type(t_convUpwind) :: rupwind
  type(t_jumpStabilisation) :: rjumpstabil
  type(t_matrixBlock), DIMENSION(1) :: Rmatrices
  type(t_matrixScalar), POINTER :: rmatrixMass
  type(t_Blockdiscretisation), POINTER :: p_rdiscretisation
  type(t_blockDiscretisation) :: rdiscretisationStabil
  type(t_discreteBC),TARGET :: rdiscreteBC
  type(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner
  type(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
  type(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  type(t_boundaryRegion) :: rboundaryRegion
  integer :: i
  integer :: ierror
  type(t_trilinearForm) :: rtriform
 
     
     
 
    
     i = rproblem%nlmax
     p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisationLS
     rmatrixMass => rproblem%RlevelInfo(i)%rmatrixLS
     
                               
     call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
     call lsyssc_duplicateMatrix (rmatrixMass,&
              rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_copy)
     call lsysbl_updateMatStrucInfo (rmatrix)
     
     rtriform%itermcount=1
     rtriform%BconstantCoeff= .true.
     rtriform%ballcoeffconstant = .true.
     rtriform%Dcoefficients(1)= rtimestepping%dweightMatrixLHS
     rtriform%Idescriptors(1,1)=DER_FUNC
     rtriform%Idescriptors(2,1)=DER_DERIV_X
     rtriform%Idescriptors(3,1)=DER_FUNC
    
     call trilf_buildMatrixScalar (rtriform,.FALSE.,&
         rmatrix%RmatrixBlock(1,1), rvectorvelocity_new%rvectorBlock(1))
     
     rtriform%Idescriptors(2,1)=DER_DERIV_Y
     call trilf_buildMatrixScalar (rtriform,.FALSE.,&
         rmatrix%RmatrixBlock(1,1), rvectorvelocity_new%rvectorBlock(2))
             
    
     call spdiscr_deriveBlockDiscr (p_rdiscretisation,rdiscretisationStabil)
     rjumpstabil%dnu = 0.0_DP
     rjumpstabil%dgamma = 0.01_DP
     rjumpstabil%dgammastar = 0.01_DP
     rjumpstabil%dtheta = rtimestepping%dweightMatrixLHS
    
     !call conv_JumpStabilisation2d ( &
     !     rjumpstabil, conv_modmatrix, rmatrix%rmatrixBlock(1,1),&
     !         rdiscretisation=rdiscretisationStabil%RspatialDiscr(1))
     
     call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix2)
     call lsyssc_duplicateMatrix (rmatrixMass,&
              rmatrix2%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_copy)
     rtriform%Dcoefficients(1)= rtimestepping%dweightMatrixRHS
     call trilf_buildMatrixScalar (rtriform,.false.,&
         rmatrix2%RmatrixBlock(1,1), rvectorvelocity_old%rvectorBlock(2))
     rtriform%Idescriptors(2,1)=DER_DERIV_X
     call trilf_buildMatrixScalar (rtriform,.FALSE.,&
         rmatrix2%RmatrixBlock(1,1), rvectorvelocity_old%rvectorBlock(1))
         
     call lsysbl_createVecBlockIndirect (rvector, rrhsBlock, .true.)
     !call lsyssc_scalarMatVec (rmatrixMass,rvector%rvectorBlock(1),&
     !           rrhsBlock%rvectorBlock(1),1.0_DP, rtimestepping%dtstep)
      
     call lsyssc_scalarMatVec (rmatrix2%rmatrixblock(1,1),rvector%rvectorBlock(1),&
                rrhsBlock%rvectorBlock(1),1.0_DP, rtimestepping%dtstep)
     
                  
 ! BC
     CALL bcasm_initDiscreteBC(rdiscreteBC)
     
     CALL boundary_createRegion(rproblem%rboundary,1,4,rboundaryRegion)
     CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_trace)
                                    
                                       
     rmatrix%p_rdiscreteBC =>  rdiscreteBC
     rrhsblock%p_rdiscreteBC => rdiscreteBC
     rvector%p_rdiscreteBC =>  rdiscreteBC
                                       
     call vecfil_discreteBCrhs (rrhsblock)
     call vecfil_discreteBCsol (rvector)
     
     call matfil_discreteBC (rmatrix)
  
     RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
     p_RfilterChain => RfilterChain
     NULLIFY(p_rpreconditioner)
     call linsol_initUMFPACK4 (p_rsolverNode)
     p_rsolverNode%ioutputLevel = 2
     Rmatrices = (/rmatrix/)
     call linsol_setMatrices(p_rsolverNode, Rmatrices)
     call linsol_initStructure(p_rsolverNode, ierror)
     call linsol_initData (p_rsolverNode,ierror)
      
     if (ierror .NE. LINSOL_ERR_NOERROR) stop
        call linsol_initData (p_rsolverNode, ierror)
     if (ierror .NE. LINSOL_ERR_NOERROR) stop
     
     call lsysbl_duplicateVector (rvector,rtempblock,&
                            LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
    
     call linsol_solveAdaptively (p_rsolverNode,rvector,rrhsblock,rtempBlock)
      
     call linsol_doneData (p_rsolverNode)
     call linsol_doneStructure (p_rsolverNode)
     call linsol_releaseSolver (p_rsolverNode)
      
     call lsysbl_releaseVector (rrhsblock)
     call lsysbl_releaseVector (rtempBlock)
      
     call lsysbl_releaseMatrix (rmatrix)
     call lsysbl_releaseMatrix (rmatrix2)
     call bcasm_releaseDiscreteBC (rdiscreteBC)
      

  
  END SUBROUTINE
                                 
  ! *********************************************************************

end module

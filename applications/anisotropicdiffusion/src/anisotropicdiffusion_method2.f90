!##############################################################################
!# ****************************************************************************
!# <name> anisotropicdiffusion_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This program is a simple test program for discretising the equation
!#
!#              - div (A grad(u) = f
!#
!# on a 2D domain for a scalar function u and a matrix
!#
!#   A =  ( cos(t) -sin(t) )  ( a11 a12 )  (  cos(t) sin(t) )
!#        ( sin(t)  cos(t) )  ( a21 a22 )  ( -sin(t) cos(t) )
!#
!# All parameters (e.g. a11,...,a22) are read from a .DAT file.
!#
!# The example  discretises and solves this equation in a direct way,
!# just listing all commands necessary  for initialisation, discretisation,
!# solving (with Gauss elimination = UMFPACK) and cleanup.
!#
!# The module uses h-adaptivity to adapt the mesh. The parameters in the
!# section [adaptivity] in the DAT file control the h-adaptivity!
!# </purpose>
!##############################################################################

module anisotropicdiffusion_method2

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
  use sortstrategy
  use coarsegridcorrection
  use scalarpde
  use ucd
  use pprocerror
  use paramlist
  use meshmodification
  use spatialdiscretisation
  use spdiscprojection
  
  use collection

  use hadaptivity
  use jumpstabilisation
    
  use anisotropicdiffusion_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine anisotropicdiffusion2
  
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
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Rotation angle, X-/Y-diffusion
    real(DP) :: dtheta,ddiffusionA11,ddiffusionA12
    real(DP) :: ddiffusionA21,ddiffusionA22
    
    ! Type of solution
    integer :: isolution
    
    ! Stabilisation
    real(DP) :: dgamma
    integer :: istabilisation
    
    ! Mesh distortion
    real(DP) :: dmeshDistortion
    
    ! Element type of the discretisation
    integer :: ieltype

    ! Type of error estimator
    integer :: ierrorestimator
    
    ! Output error as GMV file
    integer :: ioutputerror

    ! Whether to convert to triangle mesh
    integer :: iconvertToTriangleMesh
    
    ! Final diffusion matrix after rotation
    real(DP), dimension(2,2) :: DdiffusionMatrix
    
    ! PRM/TRI file
    character(LEN=SYS_STRLEN) :: sstring,sfilePRM,sfileTRI
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: dmin, dmax
    integer :: i
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir,smaster
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataQ1
    type(t_discreteBC), target :: rdiscreteBCPostProc
    type(t_vectorScalar) :: rvectorPostProc
    type(t_vectorBlock) :: rvectorPostProcBlock
    type(t_blockDiscretisation) :: rdiscretisationPostProc
    
    ! Structure for saving parameters from the DAT file
    type (t_parlist) :: rparams
    
    ! Collection structure for setting up the RHS
    type(t_collection) :: rcollection

    ! An object for saving the adaptivity structure
    type(t_hadapt) :: rhadapt
    
    ! A scalar vector for storing the indicator for h-adaptivity
    type(t_vectorScalar) :: rindicator

    ! Maximum number of h-adaptivity cycles
    integer :: MAXhRefinementSteps

    ! Set halt mode for debugging
    sys_haltmode = SYS_HALT_THROWFPE

    ! Ok, let us start.
    !
    ! +------------------------------------------------------------------------
    ! | READ DAT FILE PARAMETERS
    ! +------------------------------------------------------------------------
    !
    ! Initialise the parameter structure and read the DAT file.
    call parlst_init(rparams)

    ! Get the data file.
    call sys_getcommandLineArg(1,smaster,sdefault='./dat/anisotropicdiffusion.dat')
    call parlst_readfromfile (rparams, smaster)
    
    ! Get the parameters...
    !
    ! PRM file
    call parlst_getvalue_string (rparams, '', &
                                 'sfilePRM', sstring)
    read(sstring,*) sfilePRM
                                 
    ! TRI file
    call parlst_getvalue_string (rparams, '', &
                                 'sfileTRI', sstring)
    read(sstring,*) sfileTRI

    ! Whether to convert to a triangle mesh
    call parlst_getvalue_int (rparams, '', &
      'iconvertToTriangleMesh', iconvertToTriangleMesh, 0)

    ! Mesh distortion
    call parlst_getvalue_double (rparams, '', &
      'dmeshDistortion', dmeshDistortion, 0.0_DP)

    ! Element type
    call parlst_getvalue_int (rparams, '', &
                              'ielementType', ieltype, 7)

    ! Type of error estimator
    call parlst_getvalue_int (rparams, '', &
                              'ierrorestimator', ierrorestimator, 1)
    
    ! Output error as GMV file
    call parlst_getvalue_int (rparams, '', &
                              'ioutputerror', ioutputerror, 0)

    ! Get the path where to write gmv`s to.
    call parlst_getvalue_string (rparams, '', &
                                 'sucddir', sstring)
    read(sstring,*) sucddir

    ! Type of stabilisation
    call parlst_getvalue_int (rparams, '', &
      'istabilisation', istabilisation, 0)
    call parlst_getvalue_double (rparams, '', &
      'dgamma', dgamma, 0.01_DP)

    ! Maximum level
    call parlst_getvalue_int (rparams, '', &
                                     'NLMAX', NLMAX, 7)

    ! Maximum number of refinement steps
    call parlst_getvalue_int (rparams, '', &
                                     'MAXhRefinementSteps', MAXhRefinementSteps, 0)

    ! Rotation angle
    call parlst_getvalue_double (rparams, '', 'dtheta', dtheta, 0.0_DP)
    
    ! in rad...
    dtheta = 2.0_DP*SYS_PI*dtheta/360.0_DP
    
    ! diffusion matrix
    call parlst_getvalue_double (rparams, '', &
        'ddiffusionA11', ddiffusionA11, 1.0_DP)
    call parlst_getvalue_double (rparams, '', &
        'ddiffusionA12', ddiffusionA12, 0.0_DP)
    call parlst_getvalue_double (rparams, '', &
        'ddiffusionA21', ddiffusionA21, 0.0_DP)
    call parlst_getvalue_double (rparams, '', &
        'ddiffusionA22', ddiffusionA22, 1.0_DP)
        
    ! solution
    call parlst_getvalue_int (rparams, '', &
                              'isolution', isolution, 0)
    
    ! Create the actual diffusion matrix:
    !
    ! A = (  cos t  -sin t )  ( d11 d12 ) (  cos t   sin t )
    !     (  sin t   cos t )  ( d21 d22 ) ( -sin t   cos t )
    !
    DdiffusionMatrix(1,1) = &
      (cos(dtheta)*ddiffusionA11-sin(dtheta)*ddiffusionA21)*cos(dtheta)-&
      (cos(dtheta)*ddiffusionA12-sin(dtheta)*ddiffusionA22)*sin(dtheta)

    DdiffusionMatrix(1,2) = &
        (cos(dtheta)*ddiffusionA11-sin(dtheta)*ddiffusionA21)*sin(dtheta)+&
        (cos(dtheta)*ddiffusionA12-sin(dtheta)*ddiffusionA22)*cos(dtheta)

    DdiffusionMatrix(2,1) = &
       (sin(dtheta)*ddiffusionA11+cos(dtheta)*ddiffusionA21)*cos(dtheta)-&
       (sin(dtheta)*ddiffusionA12+cos(dtheta)*ddiffusionA22)*sin(dtheta)

    DdiffusionMatrix(2,2) = &
        (sin(dtheta)*ddiffusionA11+cos(dtheta)*ddiffusionA21)*sin(dtheta)+&
        (sin(dtheta)*ddiffusionA12+cos(dtheta)*ddiffusionA22)*cos(dtheta)
    
    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------
    !
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, sfilePRM)
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, sfileTRI, rboundary)
    
    ! Convert to triangle mesh?
    if (iconvertToTriangleMesh .eq. 1) then
      call tria_rawGridToTri(rtriangulation)
    end if
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)

    ! Mesh distortion?
    if (dmeshDistortion .ne. 0.0_DP) then
      call meshmod_disturbMesh (rtriangulation,dmeshDistortion)
    end if
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! +------------------------------------------------------------------------
    ! | SETUP OF THE H-ADAPTIVITY
    ! +------------------------------------------------------------------------
    
    call hadapt_initFromParameterlist(rhadapt,rparams,"adaptivity")
    call hadapt_initFromTriangulation(rhadapt,rtriangulation)

    ! +------------------------------------------------------------------------
    ! | SETUP OF THE DISCRETISATION
    ! +------------------------------------------------------------------------
    !
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
        rtriangulation, rboundary)
    
    ! Repeat the procedure until the maximum number of refinement
    ! steps has been reached. This will be checked below.
    do

      ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
      ! structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      select case (ieltype)
      case (1)
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
            EL_E001,SPDISC_CUB_AUTOMATIC,rtriangulation, rboundary)
      case (2)
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
            EL_E002,SPDISC_CUB_AUTOMATIC,rtriangulation, rboundary)
      case (11)
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
            EL_E011,SPDISC_CUB_AUTOMATIC,rtriangulation, rboundary)
      case (13)
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
            EL_E013,SPDISC_CUB_AUTOMATIC,rtriangulation, rboundary)
      case (-30)
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
            EL_EM30,SPDISC_CUB_AUTOMATIC,rtriangulation, rboundary)
      case (-1)
        call spdiscr_initDiscr_triquad (rdiscretisation%RspatialDiscr(1), &
            EL_E001,EL_E011,SPDISC_CUB_AUTOMATIC,SPDISC_CUB_AUTOMATIC,&
            rtriangulation,rboundary)
      case (-2)
        call spdiscr_initDiscr_triquad (rdiscretisation%RspatialDiscr(1), &
            EL_E002,EL_E013,SPDISC_CUB_AUTOMATIC,SPDISC_CUB_AUTOMATIC,&
            rtriangulation,rboundary)
      case DEFAULT
        call output_line('Unsupproted element type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'anisotropicdiffusion_method2')
        call sys_halt()
      end select

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create a scalar matrix, based on the discretisation structure
      ! for our one and only solution component.
      select case (istabilisation)
      case (1)
        ! Jump stabilisation. Needs an extended matrix stencil.
        call bilf_createMatrixStructure (&
            rdiscretisation%RspatialDiscr(1),&
            LSYSSC_MATRIX9,rmatrix,cconstrType=BILF_MATC_EDGEBASED)
      case DEFAULT
        ! No stabilisation
        call bilf_createMatrixStructure (&
            rdiscretisation%RspatialDiscr(1),&
            LSYSSC_MATRIX9,rmatrix)
      end select
      
      ! +------------------------------------------------------------------------
      ! | DISCRETISATION OF MATRICES/VECTORS
      ! +------------------------------------------------------------------------
      !
      ! And now to the entries of the matrix. For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      
      rform%itermCount = 4
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      
      rform%Idescriptors(1,2) = DER_DERIV_X
      rform%Idescriptors(2,2) = DER_DERIV_Y
      
      rform%Idescriptors(1,3) = DER_DERIV_Y
      rform%Idescriptors(2,3) = DER_DERIV_X
      
      rform%Idescriptors(1,4) = DER_DERIV_Y
      rform%Idescriptors(2,4) = DER_DERIV_Y
      
      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = DdiffusionMatrix(1,1)
      rform%Dcoefficients(2)  = DdiffusionMatrix(1,2)
      rform%Dcoefficients(3)  = DdiffusionMatrix(2,1)
      rform%Dcoefficients(4)  = DdiffusionMatrix(2,2)
      
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical
      ! data.
      select case (istabilisation)
      case (1)
        ! Jump stabilisation. CReate the Laplace matrix and add the
        ! stabilisation
        call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace)
        call jstab_calcUEOJumpStabilisation (&
            rmatrix,dgamma,dgamma,2.0_DP,1.0_DP,CUB_G3_1D,1.0_DP)
      case DEFAULT
        ! No stabilisation. Create the Laplace matrix directly.
        call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace)
      end select
      
      ! Create a collection; used for passing parameters to the RHS.
      call collct_init (rcollection)
      
      ! Put the parameters to the quick-access array
      rcollection%IquickAccess(1) = isolution
      rcollection%DquickAccess(1) = DdiffusionMatrix(1,1)
      rcollection%DquickAccess(2) = DdiffusionMatrix(1,2)
      rcollection%DquickAccess(3) = DdiffusionMatrix(2,1)
      rcollection%DquickAccess(4) = DdiffusionMatrix(2,2)
      
      ! The same has to be done for the right hand side of the problem.
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! ... and then discretise the RHS to get a discrete version of it.
      ! Again we simply create a scalar vector based on the one and only
      ! discretisation structure.
      ! This scalar vector will later be used as the one and only first
      ! component in a block vector.
      call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
          rlinform,.true.,rrhs,coeff_RHS,rcollection)
      
      ! The linear solver only works for block matrices/vectors - but above,
      ! we created scalar ones. So the next step is to make a 1x1 block
      ! system from the matrices/vectors above which the linear solver
      ! understands.
      call lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
      call lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
      
      ! +------------------------------------------------------------------------
      ! | DISCRETISATION AND IMPLEMENTATION OF BOUNDARY CONDITIONS
      ! +------------------------------------------------------------------------
      !
      ! Now we have the raw problem. What is missing is the definition of the boundary
      ! conditions.
      ! For implementing boundary conditions, we use a `filter technique with
      ! discretised boundary conditions`. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.

      call bcasm_initDiscreteBC(rdiscreteBC)

      ! Now discretise the boundary conditions. The result is saved to rdiscreteBC.
      call ad2_discretiseBC (isolution,rdiscretisation,rcollection,rdiscreteBC)

      ! Hang the pointer into the vector and matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      rmatrixBlock%p_rdiscreteBC => rdiscreteBC
      rrhsBlock%p_rdiscreteBC => rdiscreteBC
      
      ! Now we have block vectors for the RHS and the matrix. What we
      ! need additionally is a block vector for the solution and
      ! temporary data. Create them using the RHS as template.
      ! Fill the solution vector with 0:
      call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
      call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
      
      ! Next step is to implement boundary conditions into the RHS,
      ! solution and matrix. This is done using a vector/matrix filter
      ! for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! vectors/matrix. Call the appropriate vector/matrix filter that
      ! modifies the vectors/matrix according to the boundary conditions.
      call vecfil_discreteBCrhs (rrhsBlock)
      call vecfil_discreteBCsol (rvectorBlock)
      call matfil_discreteBC (rmatrixBlock)
      
      ! +------------------------------------------------------------------------
      ! | INVOKE THE SOLVER
      ! +------------------------------------------------------------------------
      !
      ! Initialise the solver
      call linsol_initUMFPACK4 (p_rsolverNode)
      
      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 2
      
      ! Attach the system matrix to the solver.
      ! First create an array with the matrix data (on all levels, but we
      ! only have one level here), then call the initialisation
      ! routine to attach all these matrices.
      ! Remark: Do not make a call like
      !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
      ! This does not work on all compilers, since the compiler would have
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
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,&
          rrhsBlock,rtempBlock)

      ! Do we have to perform one step of h-adaptivity?
      if (rhadapt%nRefinementSteps .ge. MAXhRefinementSteps) exit
      
      ! +----------------------------------------------------------------------
      ! | COMPUTE INDICATOR FOR H-ADAPTIVITY
      ! +----------------------------------------------------------------------

      ! Perform a posteriori error estimation
      call lsyssc_createVector(rindicator,rtriangulation%NEL,.true.)
      call getMonitorFunction(rtriangulation,rvectorBlock%RvectorBlock(1),&
          ieltype,ierrorestimator,rindicator)

      ! Should the error be written to GMV file
      if (ioutputerror .gt. 0) then
        call lsyssc_getbase_double(rindicator,p_Ddata)
        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,trim(sucddir)//'/u2.'//&
            trim(sys_siL(rhadapt%nRefinementSteps,3))//'.gmv')
        call ucd_addVariableElementBased (rexport,'error',UCD_VAR_STANDARD, p_Ddata)
        call ucd_write (rexport)
        call ucd_release (rexport)
      end if

      ! Perform one step h-adaptivity
      call hadapt_refreshAdaptation(rhadapt,rtriangulation)
      call hadapt_performAdaptation(rhadapt,rindicator)
      
      ! Release the indicator vector again
      call lsyssc_releaseVector(rindicator)
      
      ! Generate raw mesh from adaptivity structure
      call hadapt_generateRawMesh(rhadapt,rtriangulation)
      
      ! Create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
      
      ! +----------------------------------------------------------------------
      ! | TEMPORAL CLEANUP
      ! +----------------------------------------------------------------------

      ! Release solver data and structure
      call linsol_doneData (p_rsolverNode)
      call linsol_doneStructure (p_rsolverNode)
      
      ! Release the solver node and all subnodes attached to it (if at all):
      call linsol_releaseSolver (p_rsolverNode)
      
      ! Release the block matrix/vectors
      call lsysbl_releaseVector (rtempBlock)
      call lsysbl_releaseVector (rvectorBlock)
      call lsysbl_releaseVector (rrhsBlock)
      call lsysbl_releaseMatrix (rmatrixBlock)

      ! Release the scalar matrix/rhs vector which were used to create
      ! the block matrices/vectors. These must exist as long as the
      ! block matrices/vectors exist, as the block matrices/vectors are
      ! only 'copies' of the scalar ones, sharing the same handles!
      call lsyssc_releaseVector (rrhs)
      call lsyssc_releaseMatrix (rmatrix)

      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      
      ! Release the collection
      call collct_done (rcollection)

    end do
      
    ! +------------------------------------------------------------------------
    ! | POSTPROCESSING
    ! +------------------------------------------------------------------------
    !
    ! That is it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing.

    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    dmin = p_Ddata(1)
    dmax = dmin

    ! Start UCD export to GMV file:
    select case (ieltype)
    case (-1,1,11)
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,trim(sucddir)//'/u2.gmv')
      call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    case (2)
      call ucd_startGMV (rexport,UCD_FLAG_ONCEREFINED,rtriangulation,trim(sucddir)//'/u2.gmv')
      call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
          p_Ddata(1:rtriangulation%NVT),&
          p_Ddata(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT))
    case (13)
      call ucd_startGMV (rexport,UCD_FLAG_ONCEREFINED,rtriangulation,trim(sucddir)//'/u2.gmv')
      call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
          p_Ddata(1:rtriangulation%NVT),&
          p_Ddata(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT),&
          p_Ddata(rtriangulation%NVT+rtriangulation%NMT+1:&
            rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NEL))
    case (-2)
      call ucd_startGMV (rexport,UCD_FLAG_ONCEREFINED,rtriangulation,trim(sucddir)//'/u2.gmv')
      call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
          p_Ddata(1:rtriangulation%NVT),&
          p_Ddata(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT),&
          p_Ddata(rtriangulation%NVT+rtriangulation%NMT+1:))
    case (-30)
      !CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,trim(sucddir)//'/u2.gmv')
      call ucd_startGMV (rexport,UCD_FLAG_BULBQUADRATIC,rtriangulation,trim(sucddir)//'/u2.gmv')
      
      ! Project the solution to P1/Q1 -> rvectorPostProc
      call spdiscr_initBlockDiscr (rdiscretisationPostProc,1,&
                                   rtriangulation, rboundary)
      call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(1),&
          rvectorPostProc,rdiscretisationPostProc%RspatialDiscr(1))
          
      ! Discretise the boundary conditions and include them into the Q1
      ! solution vector rvectorPostProc using the corresponding vector filter.
      call bcasm_initDiscreteBC (rdiscreteBCPostProc)
      call ad2_discretiseBC(isolution,rdiscretisationPostProc,rcollection,&
          rdiscreteBCPostProc)

      ! Implement the discrete BC into the projected solution vector.
      call lsysbl_createVecFromScalar(rvectorPostProc,rvectorPostProcBlock)
      call vecfil_discreteBCsol (rvectorPostProcBlock,rdiscreteBCPostProc)
      call bcasm_releaseDiscreteBC (rdiscreteBCPostProc)
          
      ! Put the vector to the postprocessing file
      call lsyssc_getbase_double (rvectorPostProc,p_DdataQ1)
      call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, &
         p_DdataQ1,p_Ddata)

      do i=1,size(p_DdataQ1)
        dmin = min(dmin,p_DdataQ1(i))
        dmax = max(dmax,p_DdataQ1(i))
      end do
         
      ! Release the allocated information
      call lsysbl_releaseVector (rvectorPostProcBlock)
      call lsyssc_releaseVector (rvectorPostProc)
      call spdiscr_releaseBlockDiscr (rdiscretisationPostProc)

    case DEFAULT
      call output_line('Unsupproted element type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'anisotropicdiffusion_method2')
      call sys_halt()
    end select
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    do i=1,size(p_Ddata)
      dmin = min(dmin,p_Ddata(i))
      dmax = max(dmax,p_Ddata(i))
    end do
    
    call output_line ("Min: "//sys_sdEL(dmin,10))
    call output_line ("Max: "//sys_sdEL(dmax,10))
    
    ! +------------------------------------------------------------------------
    ! | CLEANUP
    ! +------------------------------------------------------------------------
    !
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release adaptivity structure
    call hadapt_releaseAdaptation(rhadapt)

    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseMatrix (rmatrixBlock)

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    call lsyssc_releaseVector (rrhs)
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

    ! Release the collection
    call collct_done (rcollection)

    ! Release the parameters from the DAT file.
    call parlst_done(rparams)
    
  end subroutine

  ! ---------------------------------------------------------------------------
  
!<subroutine>
  
  subroutine ad2_discretiseBC (isolution,rdiscretisation,rcollection,rdiscreteBC)
  
!<description>
  ! Discretises the current boundary conditions.
!</description>
  
!<input>
  ! Type of solution.
  integer, intent(in) :: isolution
  
  ! Current discretisation
  type(t_blockDiscretisation), intent(in), target :: rdiscretisation
!</input>

!<inputoutput>
  ! A collection structure with problem specific parameters. This structure
  ! is passed to the callback function that defines the values on the boundary.
  type(t_collection), intent(inout) :: rcollection

  ! Structure that receives the discrete boundary conditions.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_boundary), pointer :: p_rboundary
    
    ! Get a pointer to the domain
    p_rboundary => rdiscretisation%p_rboundary
  
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    
    select case (isolution)
    case (0,1,3)
      
      ! We ask the boundary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC`s for the first (here: one and only) component in the
      !   solution vector.
      ! - Discretise the boundary condition so that the BC`s can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
    case (2,4)

      ! Edge 1 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
                               
      ! Edge 2 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
      
      ! Edge 4 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)

      ! Edge 1 of boundary component 2.
      call boundary_createRegion(p_rboundary,2,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
                               
      ! Edge 2 of boundary component 2.
      call boundary_createRegion(p_rboundary,2,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
                               
      ! Edge 3 of boundary component 2.
      call boundary_createRegion(p_rboundary,2,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)
      
      ! Edge 4 of boundary component 2.
      call boundary_createRegion(p_rboundary,2,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues,rcollection)

    end select

  end subroutine
  
end module

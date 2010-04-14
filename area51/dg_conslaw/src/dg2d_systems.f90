!##############################################################################
!# ****************************************************************************
!# <name> dg2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

module dg2d_systems

  use fsystem
  use stdoperators
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
  use genoutput
  use dg2d_routines
  use collection
  use linearalgebra
  use paramlist
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine dg2d_sys
  
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
    
    ! Number of Variables (h, hu, hv)
    ! (h=Waterheights, u/v=speed in x/y-direction)
    integer, parameter :: nvar2d = 3
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinformconv, rlinformedge, rlinformIC
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixMC
    type(t_vectorScalar) :: rrhs,rsol,redge,rconv,rsoltemp,rrhstemp,rsolUp,rsolold

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock), target :: rvectorBlock,rrhsBlock,rtempBlock,rsolBlock,redgeBlock,rconvBlock,rsolTempBlock,rsolUpBlock,rsolOldBlock,rsolLimiterBlock,rsolSaveBlock
    type(t_vectorBlock), target :: rk1, rk2, rdefBlock, rimf1, rimf2


    ! A set of variables describing the discrete boundary conditions.    
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    integer :: ilimiting
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Command line and name of the paramater file
    character(LEN=SYS_STRLEN) :: cbuffer, sparameterfileName
    
    ! Strings describing initial condition and inlet condition for the function parser
    character(len=SYS_STRLEN) :: sstring, sic, sinlet
    
    real(DP) :: ttime, dt, ttfinal
    
    real(DP), dimension(2) :: vel
    
    real(dp) :: dL2updnorm
    
    integer :: ielementType
    
    type(t_collection) :: rcollection
    
    ! Parameter list
    type(t_parlist) :: rparlist
    
    character(LEN=*), dimension(2), parameter ::&
         cvariables = (/ (/'x'/), (/'y'/) /)
    
    ! How many extra points for output
    integer :: iextraPoints
    
    integer :: ivar
    
    ! For time measurement
    real(dp) :: dtime1, dtime2
    
    type(t_additionalTriaData) :: raddTriaData
    
    integer :: ilimiter
    
    real(dp) :: ddefnorm
    
    integer :: idef
    
    ! The profiler to measure the time
    type(t_profiler) :: rprofiler
    
    integer :: itimestepping
    
    integer :: iwithoutlimiting
    
    real(dp) :: dCFL, dgravconst
    
    ! Name of output file
    character (LEN=SYS_STRLEN) :: sofile
    
    integer :: imakeVideo, ifilenumber, ioutputtype
    
    real(dp) :: dvideotimestep, dvideotime
    
    integer, dimension(6) :: IdofGlob
    
    integer :: iel
     
    ! Start time measurement
    call cpu_time(dtime1)
    
        
    ! Get command line arguments and extract name of parameter file
    if (command_argument_count() .eq. 0) then
      call output_lbrk()
      call output_line('Using standart parameterfile: ./dat/1.dat')
      call output_lbrk()
      sparameterfileName = './dat/1.dat'
    else
      call get_command_argument(command_argument_count(), cbuffer)
      sparameterfileName = adjustl(cbuffer)
    end if
    
    ! Read parameter file
    call parlst_init(rparlist)
    call parlst_readFromFile(rparlist,sparameterfileName)
    
    
    ! We want to solve our problem on level... Default=1
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'NLMAX', nlmax, 1)    
    
    ! And with timestepsize
    call parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'dt', dt)

    ! To the final time
    call parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'ttfinal', ttfinal)
    
    ! To the final time
    call parlst_getvalue_double(rparlist, 'PROBLEM', 'gravconst', dgravconst)
    
    ! Type of finite element to use
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'FEkind', ielementType, 2)
    
    ! How many extra points for output
    call parlst_getvalue_int(rparlist, 'OUTPUT', 'extrapoints', iextraPoints, 3)
    
    ! Get string describing the initial condition for the function parser
    call parlst_getvalue_string (rparlist, 'PROBLEM', 'ic', sic)
    
    ! Get string describing the initial condition for the function parser
    call parlst_getvalue_string (rparlist, 'PROBLEM', 'inlet', sinlet)
    
    ! Which kind of time stepping
    call parlst_getvalue_int(rparlist, 'METHOD', 'timestepping', itimestepping, 1)
    
    ! What type of limiter to use
    call parlst_getvalue_int(rparlist, 'METHOD', 'limiter', ilimiter, 0)
    
    ! The output file
    call parlst_getvalue_string (rparlist, 'OUTPUT', 'ofile', sofile, 'gmv/u2d')
    
    ! Make gmv snapshots for video? Default=No.
    call parlst_getvalue_int(rparlist, 'OUTPUT', 'makevideo', imakevideo, 0)
    
    ! Type of output files (1 = gmv, 2 = vtk)
    call parlst_getvalue_int(rparlist, 'OUTPUT', 'outtype', ioutputtype, 1)
    
    ! Make gmv snapshot every ... seconds (should be n*dt)
    call parlst_getvalue_double(rparlist, 'OUTPUT', 'videotimestep', dvideotimestep)
    
    select case (ielementType)
    case (0)
      ielementType = EL_DG_T0_2D
      ilimiting = 0
      dCFL = 1.0_dp
    case (1)
      ielementType = EL_DG_T1_2D
      ilimiting = 1
      dCFL = 0.3_dp
    case (2)
      ielementType = EL_DG_T2_2D
      ilimiting = 2
      dCFL = 0.18_dp
    end select
    
    ! For video files
    ifilenumber = -1
    dvideotime = 0.0_dp
    
    
    
!    dt = 0.01_DP
!    ttfinal = 0.0_dp
!    !ttfinal = 5*SYS_PI
!    
!    ielementType = EL_DG_T1_2D
!    
!    
!    ilimiting = 1
!    
!        
!    vel(1)=1.0_DP
!    vel(2)=1.0_DP
!
!    ! Ok, let us start. 
!    !
!    ! We want to solve our problem on level...
!    NLMAX = 4

    ! Read in parametrisation of the boundary
    call parlst_getvalue_string (rparlist, 'TRIANGULATION', &
         'prmname', sstring)
    call boundary_read_prm(rboundary, sstring)
    
    ! Now read in the basic triangulation.
    call parlst_getvalue_string (rparlist, 'TRIANGULATION', &
         'triname', sstring)
    call tria_readTriFile2D (rtriangulation, sstring, rboundary, .true.)    
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! Calculate additional triangulation data, as the normal vectors and local edge numbers
    call addTriaData(rtriangulation,raddTriaData)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. For the shallow water problem we need three (nvar2d) blocks
    call spdiscr_initBlockDiscr (rdiscretisation, nvar2d,&
                                 rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   ielementType,CUB_G6_2D,rtriangulation, rboundary)
                                   
    ! Now copy this initialised block into the other ones
    ! But only create a derived copy, which shares the handles of the original one
    do ivar = 2, nvar2d
       call spdiscr_duplicateDiscrSc (rdiscretisation%Rspatialdiscr(1), &
            rdiscretisation%Rspatialdiscr(ivar), &
            .true.)
    end do
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixMC)!,&
                                     !rdiscretisation%RspatialDiscr(1),&
                                     !BILF_MATC_EDGEBASED)
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.

    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    
    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = 1.0 
    rform%Dcoefficients(2)  = 1.0 

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    call bilf_buildMatrixScalar (rform,.true.,rmatrixMC,coeff_Laplace_2D)
    
    ! Now create the preconditioner block matrix P
    ! First create an empty block matrix structur with nvar2d x nvar2d blocks
    call lsysbl_createEmptyMatrix (rmatrixBlock, nvar2d)

    ! Next create the diagonal blocks of P as empty matrices, using the 
    ! matrix structur of the matrix MC
    ! We will only need the diagonal blocks, as we employ a block jacobi
    ! method here
    do ivar = 1, nvar2d
      call lsyssc_duplicateMatrix (rmatrixMC, &
           rmatrixBlock%Rmatrixblock(ivar,ivar), &
           LSYSSC_DUP_SHARE, &
           !LSYSSC_DUP_EMPTY) !!!!!!!!!!
           LSYSSC_DUP_SHARE)
    end do
    
    

    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinformedge%itermCount = 1
    rlinformedge%Idescriptors(1) = DER_FUNC2D
    
    
    
    ! First create the block vectors
    ! rtempBlock, rrhsBlock, rsolBlock, rdefBlock and rstempBlock, rsolDotBlock
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rrhsBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,redgeBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rconvBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolTempBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolLimiterBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolOldBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolUpBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rTempBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rk1,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rk2,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rdefBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rimf1,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rimf2,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolSaveBlock,.true.,&
         ST_DOUBLE)
         
    
!    ! Now we have the raw problem. What is missing is the definition of the boundary
!    ! conditions.
!    ! For implementing boundary conditions, we use a `filter technique with
!    ! discretised boundary conditions`. This means, we first have to calculate
!    ! a discrete version of the analytic BC, which we can implement into the
!    ! solution/RHS vectors using the corresponding filter.
!    !
!    ! Create a t_discreteBC structure where we store all discretised boundary
!    ! conditions.
!    call bcasm_initDiscreteBC(rdiscreteBC)
!    !
!    ! We 'know' already (from the problem definition) that we have four boundary
!    ! segments in the domain. Each of these, we want to use for enforcing
!    ! some kind of boundary condition.
!    !
!    ! We ask the bondary routines to create a 'boundary region' - which is
!    ! simply a part of the boundary corresponding to a boundary segment.
!    ! A boundary region roughly contains the type, the min/max parameter value
!    ! and whether the endpoints are inside the region or not.
!    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
!    
!    ! We use this boundary region and specify that we want to have Dirichlet
!    ! boundary there. The following call does the following:
!    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
!    !   We specify icomponent='1' to indicate that we set up the
!    !   Dirichlet BC`s for the first (here: one and only) component in the 
!    !   solution vector.
!    ! - Discretise the boundary condition so that the BC`s can be applied
!    !   to matrices and vectors
!    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!                             
!    ! Now to the edge 2 of boundary component 1 the domain.
!    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!                             
!    ! Edge 3 of boundary component 1.
!    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!    
!    ! Edge 4 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!                             
!    ! Hang the pointer into the vector and matrix. That way, these
!    ! boundary conditions are always connected to that matrix and that
!    ! vector.
!    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
!    rrhsBlock%p_rdiscreteBC => rdiscreteBC
!                             
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    !call lsysbl_createVecBlockIndirect (rrhsBlock, rsolBlock, .true.)
    !call lsysbl_createVecBlockIndirect (rrhsBlock, rsolTempBlock, .true.)
    !call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
!    
!    ! Next step is to implement boundary conditions into the RHS,
!    ! solution and matrix. This is done using a vector/matrix filter
!    ! for discrete boundary conditions.
!    ! The discrete boundary conditions are already attached to the
!    ! vectors/matrix. Call the appropriate vector/matrix filter that
!    ! modifies the vectors/matrix according to the boundary conditions.
!    call vecfil_discreteBCrhs (rrhsBlock)
!    call vecfil_discreteBCsol (rvectorBlock)
!    call matfil_discreteBC (rmatrixBlock)
!
!    ! During the linear solver, the boundary conditions are also
!    ! frequently imposed to the vectors. But as the linear solver
!    ! does not work with the actual solution vectors but with
!    ! defect vectors instead.
!    ! So, set up a filter chain that filters the defect vector
!    ! during the solution process to implement discrete boundary conditions.
!    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
!
    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    nullify(p_rpreconditioner)
    !call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    call linsol_initUMFPACK4 (p_rsolverNode)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 0
    
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
    
    
    ! Initialise the collection
    call collct_init(Rcollection)
    
    
    
do iwithoutlimiting = 1,1
if (iwithoutlimiting==2) ilimiter = 0
    
    ! Now set the initial conditions via L2 projection
    rlinformIC%itermCount = 1
    rlinformIC%Idescriptors(1) = DER_FUNC2D
    !rcollection%SquickAccess(2) = cvariables
    rcollection%SquickAccess(1) = sic

    
    
    do ivar = 1, nvar2d
          
      rcollection%IquickAccess(1) = ivar
         
      call linf_buildVectorScalar2 (rlinformIC, .true., rrhsBlock%RvectorBlock(ivar),&
                                  coeff_RHS_IC, rcollection)

    end do
    
    call linsol_solveAdaptively (p_rsolverNode,rsolBlock,rrhsBlock,rtempBlock)    
    
    
!    ! Kill all quadratic parts (needed, if unsteady initial condition is applied with dg_t2)
!    call lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
!    do iel = 1, size(p_Ddata,1)/6
!    ! Get global DOFs of the element
!    call dof_locGlobMapping(rdiscretisation%RspatialDiscr(1), iel, IdofGlob(:))
!        p_Ddata(IdofGlob(4:6)) = 0.0_dp
!    end do
    
    
    
    rlinformconv%itermCount = 2
    rlinformconv%Idescriptors(1) = DER_DERIV_X
    rlinformconv%Idescriptors(2) = DER_DERIV_Y
    
    if (ilimiter .eq. 4) call dg_linearLimiterBlockIndicatorVar (rsolBlock, 1)
    if (ilimiter .eq. 6) call dg_quadraticLimiterBlockIndicatorVar (rsolBlock, 1)
    if (ilimiter .eq. 7) call dg_quadraticLimiterBlockIndicatorVar_2 (rsolBlock, 1)
    if (ilimiter .eq. 5) call dg_linearLimiterBlockCharVar (rsolBlock)
    if (ilimiter .eq. 8) call dg_quadraticLimiterBlockCharVar (rsolBlock, raddTriaData)
    if (ilimiter .eq. 9) call dg_linearLimiterBlockCharVar_mixedJacobian (rsolBlock)
    if (ilimiter .eq. 10) call dg_quadraticLimiterBlockCharVar_mixedJacobian (rsolBlock, raddTriaData)
    if (ilimiter .eq. 11) call dg_kuzminLimiterBlockCharVar_mixedJacobian (rsolBlock, raddTriaData)
        
    ! Write first video file (the initial conditions)
    ! If we want to make a video
    if (imakevideo == 1) then

      write(*,*) ''
      write(*,*) 'Writing videofile'
      
      ifilenumber = ifilenumber + 1
      
      select case (ioutputtype)
        case (1)
          ! Output solution to gmv file
          call dg2gmv(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
        case (2)
          ! Output solution to vtk file
          call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
      end select
      
      dvideotime = dvideotimestep
    end if
    
    
    
   select case (itimestepping)
   case (1)
    ttime = 0.0_DP
    
    ! Initialise the profiler with 5 timers
    call profiler_init(rprofiler, 5)
    
    if (ttfinal > 0.0_dp)then
    timestepping: do
    
       call getDtByCfl (rsolBlock, raddTriaData, dCFL, dt, dgravconst)
       dt= 0.5_dp*dt
       
       
       if (dt>ttfinal-ttime) dt = ttfinal-ttime

       ! Compute solution from time step t^n to time step t^{n+1}
       write(*,*)
       write(*,*)
       write(*,*) 'TIME STEP: t =', ttime,'dt =', dt
       write(*,*)
       
       call profiler_measure(rprofiler,1)
       
!       call lsyssc_getbase_double (rsol,p_Ddata)
!       p_Ddata(1)=0.0_DP       
       
       call lsysbl_copyVector(rsolBlock,rsoltempBlock)
       call lsysbl_copyVector(rsolBlock,rsolOldBlock)
       
       ! Step 1/3
       
       ! Create RHS-Vector
             
       ! First use the dg-function for the edge terms
       !rcollection%p_rvectorQuickAccess1 => rsolTempBlock
       !rcollection%SquickAccess(1) = sinlet
       call profiler_measure(rprofiler,2)
       call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhsBlock,rsolTempBlock,&
                                              raddTriaData,&
                                              flux_dg_buildVectorBlEdge2D_sim,&
                                              rcollection)
       
       call profiler_measure(rprofiler,1)
       call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
       ! Then add the convection terms
       call profiler_measure(rprofiler,3)
       if(ielementType .ne. EL_DG_T0_2D) then
       
         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
         call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                                       flux_sys_block,rcollection)
                                              
!         do ivar = 1, nvar2d
!         
!           rcollection%p_rvectorQuickAccess1 => rsolTempBlock
!           rcollection%IquickAccess(1) = ivar
!         
!
!           call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(ivar),&
!                                         flux_sys,rcollection)
!           !call lsyssc_vectorLinearComb (rrhstemp,rrhs,1.0_DP,1.0_DP)
!         end do
                                            
       end if
       
       call profiler_measure(rprofiler,4)
       ! Solve for solution update
       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
       
       call profiler_measure(rprofiler,1)
       ! Get new temp solution
       call lsysbl_vectorLinearComb (rsolBlock,rsolUpBlock,1.0_DP,dt,rsoltempBlock)
       
       call profiler_measure(rprofiler,5)
              ! Limit the solution vector
       !if (ilimiting.eq.1) call dg_linearLimiter (rsoltemp)
       !if (ilimiting.eq.2) call dg_quadraticLimiter (rsoltemp)
       if (ilimiter .eq. 4) call dg_linearLimiterBlockIndicatorVar (rsoltempBlock, 1)
       if (ilimiter .eq. 6) call dg_quadraticLimiterBlockIndicatorVar (rsoltempBlock, 1)
       if (ilimiter .eq. 7) call dg_quadraticLimiterBlockIndicatorVar_2 (rsoltempBlock, 1)
       if (ilimiter .eq. 5) call dg_linearLimiterBlockCharVar (rsoltempBlock)
       if (ilimiter .eq. 8) call dg_quadraticLimiterBlockCharVar (rsolTempBlock, raddTriaData)
       if (ilimiter .eq. 9) call dg_linearLimiterBlockCharVar_mixedJacobian (rsolTempBlock)
       if (ilimiter .eq. 10) call dg_quadraticLimiterBlockCharVar_mixedJacobian (rsolTempBlock, raddTriaData)
       if (ilimiter .eq. 11) call dg_kuzminLimiterBlockCharVar_mixedJacobian (rsolTempBlock, raddTriaData)
    
       
!       ! If we just wanted to use explicit euler, we would use this line instead of step 2 and 3
!       call lsysbl_copyVector (rsoltempBlock,rsolBlock)
       
       
       ! Step 2/3
       
       ! Create RHS-Vector
             
       call profiler_measure(rprofiler,2)
       ! First use the dg-function for the edge terms
       rcollection%p_rvectorQuickAccess1 => rsolTempBlock
       rcollection%SquickAccess(1) = sinlet
       call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhsBlock,rsolTempBlock,&
                                              raddTriaData,&
                                              flux_dg_buildVectorBlEdge2D_sim,&
                                              rcollection)
                                              
       call profiler_measure(rprofiler,1)
       call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
       
       ! Then add the convection terms
       call profiler_measure(rprofiler,3)
      if(ielementType .ne. EL_DG_T0_2D) then
       
         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
         call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                                       flux_sys_block,rcollection)
!         do ivar = 1, nvar2d
!         
!           rcollection%p_rvectorQuickAccess1 => rsolTempBlock
!           rcollection%IquickAccess(1) = ivar
!         
!
!           call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(ivar),&
!                                         flux_sys,rcollection)
!           !call lsyssc_vectorLinearComb (rrhstemp,rrhs,1.0_DP,1.0_DP)
!         end do
                                            
       end if
              
       ! Solve for solution update
       call profiler_measure(rprofiler,4)
       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
       
       call profiler_measure(rprofiler,1)
       ! Get new temp solution
       call lsysbl_vectorLinearComb (rsoltempBlock,rsolUpBlock,1.0_DP,dt)
       call lsysbl_vectorLinearComb (rsolBlock,rsolUpBlock,0.75_DP,0.25_DP,rsoltempBlock)
       
       call profiler_measure(rprofiler,5)
!              ! Limit the solution vector
!       if (ilimiting.eq.1) call dg_linearLimiter (rsoltemp)
!       if (ilimiting.eq.2) call dg_quadraticLimiter (rsoltemp)
        if (ilimiter .eq. 4) call dg_linearLimiterBlockIndicatorVar (rsoltempBlock, 1)
        if (ilimiter .eq. 6) call dg_quadraticLimiterBlockIndicatorVar (rsoltempBlock, 1)
        if (ilimiter .eq. 7) call dg_quadraticLimiterBlockIndicatorVar_2 (rsoltempBlock, 1)
        if (ilimiter .eq. 5) call dg_linearLimiterBlockCharVar (rsoltempBlock)
        if (ilimiter .eq. 8) call dg_quadraticLimiterBlockCharVar (rsolTempBlock, raddTriaData)
        if (ilimiter .eq. 9) call dg_linearLimiterBlockCharVar_mixedJacobian (rsolTempBlock)
        if (ilimiter .eq. 10) call dg_quadraticLimiterBlockCharVar_mixedJacobian (rsolTempBlock, raddTriaData)
        if (ilimiter .eq. 11) call dg_kuzminLimiterBlockCharVar_mixedJacobian (rsolTempBlock, raddTriaData)


       ! Step 3/3
       
       ! Create RHS-Vector
             
       ! First use the dg-function for the edge terms
       rcollection%p_rvectorQuickAccess1 => rsolTempBlock
       rcollection%SquickAccess(1) = sinlet
       call profiler_measure(rprofiler,2)
       call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhsBlock,rsolTempBlock,&
                                              raddTriaData,&
                                              flux_dg_buildVectorBlEdge2D_sim,&
                                              rcollection)
       call profiler_measure(rprofiler,1)
       call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
       
       ! Then add the convection terms
       call profiler_measure(rprofiler,3)
      if(ielementType .ne. EL_DG_T0_2D) then
         
         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
         call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                                       flux_sys_block,rcollection)
         
!         do ivar = 1, nvar2d
!         
!           rcollection%p_rvectorQuickAccess1 => rsolTempBlock
!           rcollection%IquickAccess(1) = ivar
!         
!
!           call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(ivar),&
!                                         flux_sys,rcollection)
!           !call lsyssc_vectorLinearComb (rrhstemp,rrhs,1.0_DP,1.0_DP)
!         end do
                                            
       end if
              
       ! Solve for solution update
       call profiler_measure(rprofiler,4)
       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
       
       ! Get new temp solution
       call profiler_measure(rprofiler,1)
       call lsysbl_vectorLinearComb (rsoltempBlock,rsolUpBlock,1.0_DP,dt)
       call lsysbl_vectorLinearComb (rsolUpBlock,rsolBlock,2.0_DP/3.0_DP,1.0_DP/3.0_DP)       
       
!       ! Limit the solution vector
        call profiler_measure(rprofiler,5)
!       if (ilimiting.eq.1) call dg_linearLimiter (rsol)
!       if (ilimiting.eq.2) call dg_quadraticLimiter (rsol)
       if (ilimiter .eq. 4) call dg_linearLimiterBlockIndicatorVar (rsolBlock, 1)
       if (ilimiter .eq. 6) call dg_quadraticLimiterBlockIndicatorVar (rsolBlock, 1)
       if (ilimiter .eq. 7) call dg_quadraticLimiterBlockIndicatorVar_2 (rsolBlock, 1)
       if (ilimiter .eq. 5) call dg_linearLimiterBlockCharVar (rsolBlock)
       if (ilimiter .eq. 8) call dg_quadraticLimiterBlockCharVar (rsolBlock, raddTriaData)
       if (ilimiter .eq. 9) call dg_linearLimiterBlockCharVar_mixedJacobian (rsolBlock)
       if (ilimiter .eq. 10) call dg_quadraticLimiterBlockCharVar_mixedJacobian (rsolBlock, raddTriaData)
       if (ilimiter .eq. 11) call dg_kuzminLimiterBlockCharVar_mixedJacobian (rsolBlock, raddTriaData)

!       
!       ! Test, if the solution has converged
!       call lsyssc_vectorLinearComb (rsol,rsolOld,-1.0_DP,1.0_dp)
!       dL2updnorm = lsyssc_vectorNorm (rsolOld,LINALG_NORML2) /dt/lsyssc_vectorNorm (rsol,LINALG_NORML2)
!       write(*,*) dL2updnorm
       
    
    
    call profiler_measure(rprofiler,1)
       ! Go on to the next time step
       ttime = ttime + dt
       ! If we would go beyond the final time in our next time step,
       ! then reduce the timestep
       !if (ttfinal-ttime<dt) dt = ttfinal-ttime
       
       ! If we want to make a video
        if ((imakevideo == 1).and.(ttime>dvideotime-0.001_DP*dt)) then

          write(*,*) ''
          write(*,*) 'Writing videofile'
          
          ifilenumber = ifilenumber + 1
          
          select case (ioutputtype)
            case (1)
              ! Output solution to gmv file
              call dg2gmv(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
            case (2)
              ! Output solution to vtk file
              call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
          end select
          
          dvideotime = dvideotime + dvideotimestep
        end if

       ! Leave the time stepping loop if final time is reached
       if (ttime .ge. ttfinal-0.001_DP*dt) exit timestepping
       
       !if (dL2updnorm.le.1.0e-6) exit timestepping

    end do timestepping
   end if
   
   ! Release the profiler and print statistics
   call profiler_release(rprofiler)
   
   
   
 case(2)  
   
  
  ! Initialise the profiler with 5 timers
  call profiler_init(rprofiler, 5)
   
  ttime = 0.0_DP
    
  if (ttfinal > 0.0_dp)then
    timestepping2: do

      ! Compute solution from time step t^n to time step t^{n+1}
      write(*,*)
      write(*,*)
      write(*,*) 'TIME STEP:', ttime
      write(*,*)
      
      idef = 0
      
      ! Calculate rimf2
      call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                          rrhsBlock,rsolBlock,&
                                           raddTriaData,&
                                           flux_dg_buildVectorBlEdge2D_sim,&
                                           rcollection)
       
      call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
       
      if(ielementType .ne. EL_DG_T0_2D) then
       
        rcollection%p_rvectorQuickAccess1 => rsolBlock
        call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                                     flux_sys_block,rcollection)
                                            
      end if
       
      call linsol_solveAdaptively (p_rsolverNode,rimf2,rrhsBlock,rtempBlock)
      
      
      
      ! Calculate K1 and K2 
      call lsysbl_vectorLinearComb (rsolBlock,rimf2,0.5_dp,dt/8.0_dp,rk1)
      call lsysbl_vectorLinearComb (rsolBlock,rimf2,1.0_dp,dt/6.0_dp,rk2)
      
      defcorr: do
      
      ! Calculate U1
      call lsysbl_vectorLinearComb (rsolBlock,rk1,0.5_dp,1.0_dp,rsolTempBlock)
      call lsysbl_vectorLinearComb (rimf2,rsolTempBlock,-dt/8.0_dp,1.0_dp)
      
      ! calculate rimf1 
      call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhsBlock,rsolTempBlock,&
                                              raddTriaData,&
                                              flux_dg_buildVectorBlEdge2D_sim,&
                                              rcollection)
       
       call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
       
       if(ielementType .ne. EL_DG_T0_2D) then
       
         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
         call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                                       flux_sys_block,rcollection)
                                            
       end if
       
       call linsol_solveAdaptively (p_rsolverNode,rimf1,rrhsBlock,rtempBlock)
      
      
      ! Calculate U2
      call lsysbl_vectorLinearComb (rk2,rimf2,1.0_dp,dt/6.0_dp,rsolBlock)
      call lsysbl_vectorLinearComb (rimf1,rsolBlock,2.0_dp/3.0_dp*dt,1.0_dp)
      
      ! Calculate rimf2
      call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhsBlock,rsolBlock,&
                                              raddTriaData,&
                                              flux_dg_buildVectorBlEdge2D_sim,&
                                              rcollection)
       
       call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
       
       if(ielementType .ne. EL_DG_T0_2D) then
       
         rcollection%p_rvectorQuickAccess1 => rsolBlock
         call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                                       flux_sys_block,rcollection)
                                            
       end if
       
       call linsol_solveAdaptively (p_rsolverNode,rimf2,rrhsBlock,rtempBlock)
      
      ! Calculate defect
      call lsysbl_vectorLinearComb (rk2,rimf2,1.0_dp,dt/6.0_dp,rdefBlock)
      call lsysbl_vectorLinearComb (rimf1,rdefBlock,2.0_dp/3.0_dp*dt,1.0_dp)
      call lsysbl_vectorLinearComb (rsolBlock,rdefBlock,-1.0_dp,1.0_dp)
      
      ddefNorm = lsysbl_vectorNorm (rdefBlock,LINALG_NORML2)
      write(*,*) '   Defect:', ddefNorm
      
      if (ddefnorm < 0.00000000001) exit defcorr
      
      idef = idef +1
      
      if (idef > 100) exit defcorr
      
      end do defcorr
      
            ! Limit the solution vector
      !if (ilimiting.eq.1) call dg_linearLimiter (rsoltemp)
      !if (ilimiting.eq.2) call dg_quadraticLimiter (rsoltemp)
      if (ilimiter .eq. 4) call dg_linearLimiterBlockIndicatorVar (rsolBlock, 1)
      if (ilimiter .eq. 6) call dg_quadraticLimiterBlockIndicatorVar (rsolBlock, 1)
      if (ilimiter .eq. 7) call dg_quadraticLimiterBlockIndicatorVar_2 (rsolBlock, 1)
      if (ilimiter .eq. 5) call dg_linearLimiterBlockCharVar (rsolBlock)
      if (ilimiter .eq. 8) call dg_quadraticLimiterBlockCharVar (rsolBlock, raddTriaData)
      if (ilimiter .eq. 9) call dg_linearLimiterBlockCharVar_mixedJacobian (rsolBlock)
      if (ilimiter .eq. 10) call dg_quadraticLimiterBlockCharVar_mixedJacobian (rsolBlock, raddTriaData)
      if (ilimiter .eq. 11) call dg_kuzminLimiterBlockCharVar_mixedJacobian (rsolBlock, raddTriaData)



     ! Go on to the next time step
     ttime = ttime + dt
     ! If we would go beyond the final time in our next time step,
     ! then reduce the timestep
     !if (ttfinal-ttime<dt) dt = ttfinal-ttime

     ! Leave the time stepping loop if final time is reached
     if (ttime .ge. ttfinal-0.001_DP*dt) exit timestepping2
       

    end do timestepping2
   end if
   
   
      ! Release the profiler and print statistics
   call profiler_release(rprofiler)
   end select
   
   
   
   
   
   if (iwithoutlimiting==1) call lsysbl_copyVector(rsolBlock,rsolSaveBlock)
   
   
 end do ! iwithoutlimiting
   
   
   
   
   
   
   
   
   
    
    
!    
!    
!    
!    
!    call lsyssc_copyVector(rsol,rsoltemp)
!    call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G3_1D, .true.,&
!                                              rrhs,rsolTemp,&
!                                              flux_dg_buildVectorScEdge2D_sim)
!       

!       ! Solve for solution update
!       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
!       
!       ! Get new temp solution
!       call lsyssc_vectorLinearComb (rsol,rsolUp,1.0_DP,dt,rsoltemp)
!    call lsyssc_copyVector(rsoltemp,rsol)
!    
!    
    
    
    write(*,*) ''
    write(*,*) 'Writing solution to file'
      
    ifilenumber = -1
     
    !select case (ioutputtype)
    !  case (1)
        ! Output solution to gmv file
        call dg2gmv(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
    !  case (2)
        ! Output solution to vtk file
        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
    !end select
      
    write(*,*) 'Writing steady solution to file'
    ! And output the steady projection
     call dg_proj2steady(rsolBlock,rtriangulation, rboundary)
     
    ! Saving the solution DOFs to file
    write(*,*) 'Writing DOFs to file'
    call saveSolutionData(rsolBlock%Rvectorblock(1),sofile,ifilenumber)

    
!    sofile = 'l8'
!    call loadSolutionData(rsolBlock%Rvectorblock(1),sofile)
   
 
    ! Calculate the error to the reference function.
    call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('L1-error: ' // sys_sdEL(derror,10) )

    ! Calculate the error to the reference function.
    call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )    


!    ! Calculate the difference of limited and non-limited solution
!    rcollection%p_rvectorQuickAccess1 => rsolSaveBlock
!    call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror,&
!                       getCompareFunction_2D,rcollection)
!    call output_line ('L1-error: ' // sys_sdEL(derror,10) )
!
!    ! Calculate the error to the reference function.
!    call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L2ERROR,derror,&
!                       getCompareFunction_2D,rcollection)
!    call output_line ('L2-error: ' // sys_sdEL(derror,10) )    

    
!
!    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
!                       getReferenceFunction_2D)
!    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
!    
!    ! We are finished - but not completely!
!    ! Now, clean up so that all the memory is available again.
!    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    

    ! Initialise the collection
    call collct_done(rcollection)

!    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseVector (rsolBlock)
    call lsysbl_releaseVector (rsolLimiterBlock)
    call lsysbl_releaseVector (rsolTempBlock)
    call lsysbl_releaseVector (rsolUpBlock)

    call lsysbl_releaseMatrix (rmatrixBlock)

!    call lsyssc_releaseVector (rrhs)
!    call lsyssc_releaseVector (redge)
!    call lsyssc_releaseVector (rconv)
!    call lsyssc_releaseVector (rsol)
!    call lsyssc_releaseVector (rsoltemp)
!    call lsyssc_releaseVector (rrhstemp)
!    call lsyssc_releaseVector (rsolUp)
!    call lsyssc_releaseVector (rsolOld)
    
    call lsyssc_releaseMatrix (rmatrixMC)
    
    
!
!    ! Release the scalar matrix/rhs vector which were used to create
!    ! the block matrices/vectors. These must exist as long as the
!    ! block matrices/vectors exist, as the block matrices/vectors are
!    ! only 'copies' of the scalar ones, sharing the same handles!
!    call lsyssc_releaseVector (rrhs)
!    call lsyssc_releaseMatrix (rmatrix)
!    
!    ! Release our discrete version of the boundary conditions
!    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    call tria_done (rtriangulation)
    
    ! Release additional triangulation data, as the normal vectors and local edge numbers
    call releaseTriaData(raddTriaData)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
    ! End of time measurement
    call cpu_time(dtime2)
    write (*,*) 'Calculation took', dtime2-dtime1
    
  end subroutine

end module

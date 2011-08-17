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
  use matrixio

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

    ! Number of Variables 
    ! Shallow water : 3 (h, hu, hv)
    ! (h=Waterheights, u/v=speed in x/y-direction)
    ! Euler: 4 (rho, rho u, rho v, rho E)
    integer, parameter :: nvar2d = 1

    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinformconv, rlinformedge, rlinformIC, rlinformSource, rlinform

    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixMC, rmatrixiMC, rmatrixA
    type(t_vectorScalar) :: rrhs,rsol,redge,rconv,rsoltemp,rrhstemp,rsolUp,rsolold,rdef

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock, rmatrixiBlock, rmatrixABlock
    type(t_vectorBlock), target :: rvectorBlock,rrhsBlock,rtempBlock,rsolBlock,redgeBlock,rconvBlock,rsolTempBlock,rsolUpBlock,rsolOldBlock,rsolLimiterBlock,rsolSaveBlock,rsourceTermBlock
    type(t_vectorBlock), target :: rk0, rk1, rk2, rk3, rdefBlock, rimf1, rimf2


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

    integer :: iel, NEL

    integer :: irhstype

    real(dp) :: dtheta

    integer :: iinsertSourceTerm = 0

    integer :: ilimitEveryStep = 1

    real(dp) , dimension(:), pointer :: p_DiMCdata, p_DMCdata

    integer :: i, j
    
    integer :: inumSolverCalls = 0, iiterations =0

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

    ! Theta for theta-scheme
    call parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'theta', dtheta)

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

    ! What type of building the rhs
    call parlst_getvalue_int(rparlist, 'METHOD', 'rhstype', irhstype)

    ! The output file
    call parlst_getvalue_string (rparlist, 'OUTPUT', 'ofile', sofile, 'gmv/u2d')

    ! Make snapshots for video? Default=No.
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
    case (3)
       ielementtype = EL_DG_Q1_2D
    case (4)
       ielementtype = EL_DG_Q2_2D
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
    !CUB_G3x3

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
            !LSYSSC_DUP_EMPTY)
       LSYSSC_DUP_SHARE)
    end do

    ! Now we calculate the inverse matrix

    ! Calculate the inverse of the mass matrix
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixiMC,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
    call lsyssc_copyMatrix (rmatrixMC, rmatrixiMC)
    if ((ielementType == EL_DG_T2_2D).or.(ielementType == EL_DG_T1_2D).or.(ielementType == EL_DG_T0_2D))&
    call dg_invertMassMatrix(rmatrixiMC)  

    !    ! Output MC and iMC
    !    call matio_writeMatrixHR (rmatrixiMC, 'iMC',&
    !                                  .true., 0, './gmv/iMC.txt', '(E20.10)')
    !    call matio_writeMatrixHR (rmatrixMC, 'MC',&
    !                                  .true., 0, './gmv/MC.txt', '(E20.10)')


    !    ! Calculate the condition number of MC    
    !    call lsyssc_getbase_double (rmatrixMC, p_DMCdata)
    !    call lsyssc_getbase_double (rmatrixiMC, p_DiMCdata)
    !    write(*,*) p_DMCdata
    !    pause
    !    write(*,*) p_DiMCdata
    !    pause
    !    write(*,*) 'Konditionszahl: ', maxval(p_DiMCdata)*maxval(p_DMCdata)
    !    pause

    ! First create an empty block matrix structur with nvar2d x nvar2d blocks
    call lsysbl_createEmptyMatrix (rmatrixiBlock, nvar2d)

    ! Next create the diagonal blocks of P as empty matrices, using the 
    ! matrix structur of the matrix MC
    ! We will only need the diagonal blocks, as we employ a block jacobi
    ! method here
    do ivar = 1, nvar2d
       call lsyssc_duplicateMatrix (rmatrixiMC, &
            rmatrixiBlock%Rmatrixblock(ivar,ivar), &
            LSYSSC_DUP_SHARE, &
            !LSYSSC_DUP_EMPTY)
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
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rk0,.true.,&
         ST_DOUBLE)     
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rk1,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rk2,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rk3,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rdefBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rimf1,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rimf2,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolSaveBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsourceTermBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rdefBlock,.true.,&
         ST_DOUBLE)

    ! Get number of elements in the triangulation
    NEL = rsolBlock%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL



    !
    ! Create a UMFPACK-solver.
    nullify(p_rpreconditioner)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    !call linsol_initUMFPACK4 (p_rsolverNode)

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! The linear solver stops, when this relative or absolut norm of
    ! the residual is reached.
    p_rsolverNode%depsRel = 1.0e-14
    p_rsolverNode%depsAbs = 1.0e-14

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




    !    ! Now calculate the source term
    !    if (iinsertSourceTerm == 1) then
    !      rsolBlock%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(1)%ccubTypeEval=CUB_G6_2d
    !      rlinformSource%itermCount = 1
    !      rlinformSource%Idescriptors(1) = DER_FUNC2D
    !      call linf_buildVectorBlock2 (rlinformSource, .true., rsourceTermBlock,&
    !                                         Callback_Sourceterm,rcollection)
    !      rsolBlock%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(1)%ccubTypeEval=CUB_G3x3
    !    end if


    do iwithoutlimiting = 1,1
       if (iwithoutlimiting==2) ilimiter = 0

       ! Now set the initial conditions via L2 projection
       rlinformIC%itermCount = 1
       rlinformIC%Idescriptors(1) = DER_FUNC2D
       !rcollection%SquickAccess(2) = cvariables
       rcollection%SquickAccess(1) = sic



       do ivar = 1, nvar2d

          rcollection%IquickAccess(1) = ivar

          !rrhsBlock%p_rblockDiscr%RspatialDiscr(ivar)%RelementDistr(1)%ccubTypeLinForm=CUB_G6_2D
          call linf_buildVectorScalar2 (rlinformIC, .true., rrhsBlock%RvectorBlock(ivar),&
               Euler_coeff_RHS_IC, rcollection)
          !rrhsBlock%p_rblockDiscr%RspatialDiscr(ivar)%RelementDistr(1)%ccubTypeLinForm=CUB_G3x3

       end do

       call linsol_solveAdaptively (p_rsolverNode,rsolBlock,rrhsBlock,rtempBlock)    


       !    ! Kill all quadratic parts (can be needed, if unsteady initial condition is applied with dg_t2)
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
       if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)

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

                !call getDtByCfl (rsolBlock, raddTriaData, dCFL, dt, dgravconst)
                !dt = 0.5_dp*dt


                if (dt>ttfinal-ttime) dt = ttfinal-ttime

                ! Compute solution from time step t^n to time step t^{n+1}
                write(*,*)
                write(*,*)
                write(*,*) 'TIME STEP: t =', ttime,'dt =', dt
                write(*,*)

                ! Set time level in collection
                rcollection%Dquickaccess(1) = ttime

                call profiler_measure(rprofiler,1)

                !       call lsyssc_getbase_double (rsol,p_Ddata)
                !       p_Ddata(1)=0.0_DP       

                call lsysbl_copyVector(rsolBlock,rsoltempBlock)
                call lsysbl_copyVector(rsolBlock,rsolOldBlock)

                ! Step 1/3

                ! Output total mass
                call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
                call output_line ('Total mass: ' // sys_sdEL(derror,10) )


                ! Create RHS-Vector

                ! First use the dg-function for the edge terms
                !rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                !rcollection%SquickAccess(1) = sinlet
                call profiler_measure(rprofiler,2)
                call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G5_1D, .true.,&
                     rrhsBlock,rsolTempBlock,&
                     raddTriaData,&
                     Euler_flux_dg_buildVectorBlEdge2D_sim,&
                     rcollection)


                write(*,*) lsyssc_vectornorm(rrhsBlock%Rvectorblock(1),1)                                     


                !       !!!!!!!!!!!!!!!!!!!!!!                                       
                !       call lsysbl_scaleVector (rrhsBlock,-1.0_DP)                                       
                !       rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                !       rcollection%SquickAccess(1) = sinlet
                !       call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G3_1D, .false.,&
                !                                              rrhsBlock%rvectorBlock(1),rsolTempBlock%rvectorBlock(1),&
                !                                              raddTriaData,&
                !                                              flux_dg_buildVectorScEdge2D_sim,&
                !                                              rcollection)
                !       call pperr_scalar (rrhsBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
                !       write(*,*) 'Difference1:', derror  
                !       
                !       ! If we want to make a video
                !        if ((ttime>dvideotime-0.001_DP*dt)) then
                !
                !          write(*,*) ''
                !          write(*,*) 'Writing videofile'
                !          
                !          ifilenumber = ifilenumber + 1
                !          
                !          select case (ioutputtype)
                !            case (1)
                !              ! Output solution to gmv file
                !              call dg2gmv(rrhsBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !            case (2)
                !              ! Output solution to vtk file
                !              call dg2vtk(rrhsBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !          end select
                !          
                !          dvideotime = dvideotime + dvideotimestep
                !        end if
                !
                !                                              
                !                                              
                !         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                !         call linf_buildVectorBlock2 (rlinformconv, .true., rrhsBlock,&
                !                                       flux_sys_block,rcollection)
                !         
                !         call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
                !         
                !         rcollection%IquickAccess(1) = 1
                !         call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(1),&
                !                                       flux_sys,rcollection)
                !                call pperr_scalar (rrhsBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
                !       write(*,*) 'Difference2:', derror
                !                                               
                !
                !       call profiler_measure(rprofiler,2)
                !       call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                !                                              rrhsBlock,rsolTempBlock,&
                !                                              raddTriaData,&
                !                                              flux_dg_buildVectorBlEdge2D_sim,&
                !                                              rcollection)
                !
                !      !!!!!!!!!!!!!!!!!!!!!!!                                        



                call profiler_measure(rprofiler,1)
                call lsysbl_scaleVector (rrhsBlock,-1.0_DP)
                ! Then add the convection terms
                call profiler_measure(rprofiler,3)
                if(ielementType .ne. EL_DG_T0_2D) then
                   select case (irhstype)
                   case (1)
                      rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                      call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                           Euler_flux_sys_block,rcollection)
                      write(*,*) lsyssc_vectornorm(rrhsBlock%Rvectorblock(1),1)
                   case(2)                             
                      do ivar = 1, nvar2d

                         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                         rcollection%IquickAccess(1) = ivar

                         call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(ivar),&
                              Euler_flux_sys_sc,rcollection)
                         !call lsyssc_vectorLinearComb (rrhstemp,rrhs,1.0_DP,1.0_DP)
                      end do
                   end select
                end if

                ! Implement source term
                if (iinsertSourceTerm==1) call lsysbl_vectorLinearComb (rsourceTermBlock,rrhsBlock,1.0_dp,1.0_dp)

                call profiler_measure(rprofiler,4)
                ! Solve for solution update
                !call lsysbl_blockMatVec (rmatrixiBlock, rrhsBlock, rsolUpBlock, 1.0_dp, 0.0_dp)
                if ((ielementType == EL_DG_T2_2D).or.(ielementType == EL_DG_T1_2D).or.(ielementType == EL_DG_T0_2D)) then
                do ivar = 1, nvar2d
                   call lsyssc_scalarMatVec (rmatrixiMC, rrhsBlock%RvectorBlock(ivar), rsolUpBlock%RvectorBlock(ivar), 1.0_dp, 0.0_dp)
                end do
                else
                   call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
                end if
                

                call profiler_measure(rprofiler,1)
                ! Get new temp solution
                call lsysbl_vectorLinearComb (rsolBlock,rsolUpBlock,1.0_DP,dt,rsoltempBlock)

                call profiler_measure(rprofiler,5)
                ! Limit the solution vector

                if (ilimitEveryStep.eq.1) then
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
                   if (ilimiter .eq. 12) call dg_realKuzmin (rsolTempBlock, raddTriaData)
                end if


                !       ! If we just wanted to use explicit euler, we would use this line instead of step 2 and 3
                !       call lsysbl_copyVector (rsoltempBlock,rsolBlock)
                  
                  
                ! Step 2/3

                ! Create RHS-Vector

                call profiler_measure(rprofiler,2)
                ! First use the dg-function for the edge terms
                rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                rcollection%SquickAccess(1) = sinlet
                call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G5_1D, .true.,&
                     rrhsBlock,rsolTempBlock,&
                     raddTriaData,&
                     Euler_flux_dg_buildVectorBlEdge2D_sim,&
                     rcollection)

                call profiler_measure(rprofiler,1)
                call lsysbl_scaleVector (rrhsBlock,-1.0_DP)

                ! Then add the convection terms
                call profiler_measure(rprofiler,3)
                if(ielementType .ne. EL_DG_T0_2D) then
                   select case (irhstype)
                   case (1)
                      rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                      call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                           Euler_flux_sys_block,rcollection)
                   case(2)
                      do ivar = 1, nvar2d

                         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                         rcollection%IquickAccess(1) = ivar


                         call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(ivar),&
                              Euler_flux_sys_sc,rcollection)
                         !call lsyssc_vectorLinearComb (rrhstemp,rrhs,1.0_DP,1.0_DP)
                      end do
                   end select

                end if

                ! Implement source term
                if (iinsertSourceTerm==1) call lsysbl_vectorLinearComb (rsourceTermBlock,rrhsBlock,1.0_dp,1.0_dp)

                ! Solve for solution update
                call profiler_measure(rprofiler,4)
                !call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
                !call lsysbl_blockMatVec (rmatrixiBlock, rrhsBlock, rsolUpBlock, 1.0_dp, 0.0_dp)
                if ((ielementType == EL_DG_T2_2D).or.(ielementType == EL_DG_T1_2D).or.(ielementType == EL_DG_T0_2D)) then
                do ivar = 1, nvar2d
                   call lsyssc_scalarMatVec (rmatrixiMC, rrhsBlock%RvectorBlock(ivar), rsolUpBlock%RvectorBlock(ivar), 1.0_dp, 0.0_dp)
                end do
                else
                   call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
                end if

                call profiler_measure(rprofiler,1)
                ! Get new temp solution
                call lsysbl_vectorLinearComb (rsoltempBlock,rsolUpBlock,1.0_DP,dt)
                call lsysbl_vectorLinearComb (rsolBlock,rsolUpBlock,0.75_DP,0.25_DP,rsoltempBlock)

                call profiler_measure(rprofiler,5)
                if (ilimitEveryStep.eq.1) then
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
                   if (ilimiter .eq. 12) call dg_realKuzmin (rsolTempBlock, raddTriaData)
                end if

                ! Step 3/3

                ! Create RHS-Vector

                ! First use the dg-function for the edge terms
                rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                rcollection%SquickAccess(1) = sinlet
                call profiler_measure(rprofiler,2)
                call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G5_1D, .true.,&
                     rrhsBlock,rsolTempBlock,&
                     raddTriaData,&
                     Euler_flux_dg_buildVectorBlEdge2D_sim,&
                     rcollection)
                call profiler_measure(rprofiler,1)
                call lsysbl_scaleVector (rrhsBlock,-1.0_DP)

                ! Then add the convection terms
                call profiler_measure(rprofiler,3)
                if(ielementType .ne. EL_DG_T0_2D) then
                   select case (irhstype)
                   case (1)
                      rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                      call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                           Euler_flux_sys_block,rcollection)
                   case(2)

                      do ivar = 1, nvar2d

                         rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                         rcollection%IquickAccess(1) = ivar


                         call linf_buildVectorScalar2 (rlinformconv, .false., rrhsBlock%RvectorBlock(ivar),&
                              Euler_flux_sys_sc,rcollection)
                         !call lsyssc_vectorLinearComb (rrhstemp,rrhs,1.0_DP,1.0_DP)
                      end do
                   end select

                end if

                ! Implement source term
                if (iinsertSourceTerm==1) call lsysbl_vectorLinearComb (rsourceTermBlock,rrhsBlock,1.0_dp,1.0_dp)

                ! Solve for solution update
                call profiler_measure(rprofiler,4)
                !call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
                !call lsysbl_blockMatVec (rmatrixiBlock, rrhsBlock, rsolUpBlock, 1.0_dp, 0.0_dp)      
                if ((ielementType == EL_DG_T2_2D).or.(ielementType == EL_DG_T1_2D).or.(ielementType == EL_DG_T0_2D)) then
                do ivar = 1, nvar2d
                   call lsyssc_scalarMatVec (rmatrixiMC, rrhsBlock%RvectorBlock(ivar), rsolUpBlock%RvectorBlock(ivar), 1.0_dp, 0.0_dp)
                end do
                else
                   call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
                end if

                !       derror = lsyssc_vectorNorm (rsolUpBlock%Rvectorblock(1),LINALG_NORML1) 
                !       call pperr_scalar (rsolUpBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
                !       write(*,*) 'SolUp:', derror
                !       
                !       
                call lsysbl_vectorLinearComb (rsolBlock,rsolOldBlock,1.0_DP,-1.0_dp)
                call lsysbl_scaleVector (rsolOldBlock,1.0_DP/dt)
                call lsyssc_scalarMatVec (rmatrixMC, rsolOldBlock%Rvectorblock(1), rrhsBlock%Rvectorblock(1), -1.0_dp, 1.0_dp)
                !       call pperr_scalar (rrhsBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
                !       write(*,*) 'Res  :', derror 


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
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)


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
                !if ((ttime .ge. ttfinal-0.001_DP*dt).or.(abs(derror)<1e-12)) exit timestepping
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
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)



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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case(3)  


          ! Initialise the profiler with 5 timers
          call profiler_init(rprofiler, 5)

          ttime = 0.0_DP

          if (ttfinal > 0.0_dp)then
             timestepping3: do

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

                ! Calculate K1
                call lsysbl_vectorLinearComb (rsolBlock,rimf2,1.0_dp,dt*(1-dtheta),rk1)

                ! Initalise temp solution with current solution
                call lsysbl_copyVector(rsolBlock,rsoltempBlock)

                cn_iteration: do

                   call lsysbl_copyVector(rsolTempBlock,rimf1)

                   ! Calculate RHS
                   call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                        rrhsBlock, rsolTempBlock,&
                        raddTriaData,&
                        flux_dg_buildVectorBlEdge2D_sim,&
                        rcollection)

                   call lsysbl_scaleVector (rrhsBlock,-1.0_DP)

                   if(ielementType .ne. EL_DG_T0_2D) then

                      rcollection%p_rvectorQuickAccess1 => rsolTempBlock
                      call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                           flux_sys_block,rcollection)

                   end if

                   call linsol_solveAdaptively (p_rsolverNode,rk2,rrhsBlock,rtempBlock)

                   ! Calculate new temp solution
                   call lsysbl_vectorLinearComb (rk1,rk2,1.0_dp,dt*dtheta,rsolTempBlock)

                   ! Calculate solution update
                   call lsysbl_vectorLinearComb (rsolTempBlock,rimf1,-1.0_dp,1.0_dp)

                   ddefNorm = lsysbl_vectorNorm (rimf1,LINALG_NORML2)/lsysbl_vectorNorm (rsolTempBlock,LINALG_NORML2)
                   write(*,*) '   Defect:', ddefNorm

                   if (ddefnorm < 1e-12) exit cn_iteration


                   idef = idef +1

                   if (idef > 20) exit cn_iteration

                end do cn_iteration

                call lsysbl_copyVector(rsolTempBlock,rsolBlock)



                ! Limit the solution vector
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)



                ! Go on to the next time step
                ttime = ttime + dt
                ! If we would go beyond the final time in our next time step,
                ! then reduce the timestep
                !if (ttfinal-ttime<dt) dt = ttfinal-ttime

                ! Leave the time stepping loop if final time is reached
                if (ttime .ge. ttfinal-0.001_DP*dt) exit timestepping3


             end do timestepping3
          end if


          ! Release the profiler and print statistics
          call profiler_release(rprofiler)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       case (4)
          ttime = 0.0_DP

          ! Initialise the profiler with 5 timers
          call profiler_init(rprofiler, 5)

          if (ttfinal > 0.0_dp)then
             timestepping4: do

                call getDtByCfl (rsolBlock, raddTriaData, dCFL, dt, dgravconst)
                ! dt = 0.5_dp*dt


                if (dt>ttfinal-ttime) dt = ttfinal-ttime

                ! Compute solution from time step t^n to time step t^{n+1}
                write(*,*)
                write(*,*)
                write(*,*) 'TIME STEP: t =', ttime,'dt =', dt
                write(*,*)

                call profiler_measure(rprofiler,1)

                call lsysbl_copyVector(rsolBlock,rsoltempBlock)
                call lsysbl_copyVector(rsolBlock,rsolOldBlock)

                ! Step 1/4

                ! Create RHS-Vector

                ! First use the dg-function for the edge terms
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
                end if

                call profiler_measure(rprofiler,4)
                ! Solve for solution update
                call linsol_solveAdaptively (p_rsolverNode,rk0,rrhsBlock,rtempBlock)

                call profiler_measure(rprofiler,1)
                ! Get new temp solution
                call lsysbl_vectorLinearComb (rsolOldBlock,rk0,1.0_DP,0.5_dp*dt,rsoltempBlock)

                call profiler_measure(rprofiler,5)

                ! Limit the solution vector
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolTempBlock, raddTriaData)


                ! Step 2/4

                ! Create RHS-Vector

                call profiler_measure(rprofiler,2)
                ! First use the dg-function for the edge terms
                rcollection%p_rvectorQuickAccess1 => rsolTempBlock
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
                end if

                ! Solve for solution update
                call profiler_measure(rprofiler,4)
                call linsol_solveAdaptively (p_rsolverNode,rk1,rrhsBlock,rtempBlock)

                call profiler_measure(rprofiler,1)
                ! Get new temp solution
                call lsysbl_vectorLinearComb (rk1,rsolOldBlock,0.5_dp*dt,1.0_DP,rsoltempBlock)

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
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolTempBlock, raddTriaData)


                ! Step 3/4

                ! Create RHS-Vector

                ! First use the dg-function for the edge terms
                rcollection%p_rvectorQuickAccess1 => rsolTempBlock
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
                end if

                ! Solve for solution update
                call profiler_measure(rprofiler,4)
                call linsol_solveAdaptively (p_rsolverNode,rk2,rrhsBlock,rtempBlock)

                ! Get new temp solution
                call profiler_measure(rprofiler,1)
                call lsysbl_vectorLinearComb (rsolOldBlock,rk2,1.0_DP,dt,rsolTempBlock)

                !       ! Limit the solution vector
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolTempBlock, raddTriaData)


                ! Step 4/4

                ! Create RHS-Vector

                ! First use the dg-function for the edge terms
                rcollection%p_rvectorQuickAccess1 => rsolTempBlock
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
                end if

                ! Solve for solution update
                call profiler_measure(rprofiler,4)
                call linsol_solveAdaptively (p_rsolverNode,rk3,rrhsBlock,rtempBlock)

                ! Get solution at next timestep
                call profiler_measure(rprofiler,1)


                call lsysbl_vectorLinearComb (rk1,rk0,2.0_dp,1.0_DP)
                call lsysbl_vectorLinearComb (rk2,rk0,2.0_dp,1.0_DP)
                call lsysbl_vectorLinearComb (rk3,rk0,1.0_dp,1.0_DP)

                call lsysbl_vectorLinearComb (rsolOldBlock,rk0,1.0_dp,1.0_dp/6.0_dp*dt,rsolBlock)

                ! Limit the solution vector
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)



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
                if (ttime .ge. ttfinal-0.001_DP*dt) exit timestepping4

                !if (dL2updnorm.le.1.0e-6) exit timestepping

             end do timestepping4
          end if

          ! Release the profiler and print statistics
          call profiler_release(rprofiler)


          ! Fully implicit linear convection
       case (5)





          ! Now as the discretisation is set up, we can start to generate
          ! the structure of the system matrix which is to solve.
          ! We create a scalar matrix, based on the discretisation structure
          ! for our one and only solution component.
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixA,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)

          ! And now to the entries of the matrix. For assembling of the entries,
          ! we need a bilinear form, which first has to be set up manually.
          ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
          ! scalar system matrix in 2D.

          rform%itermCount = 2
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_DERIV_X
          rform%Idescriptors(1,2) = DER_FUNC
          rform%Idescriptors(2,2) = DER_DERIV_Y

          ! In the standard case, we have constant coefficients:
          rform%ballCoeffConstant = .false.
          rform%BconstantCoeff = .false.

          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          !call lsyssc_clearMatrix (rmatrixMC)
          ! We need to give the solution vector to the routine (by the collection)
          rcollection%p_rvectorQuickAccess1 => rsolBlock
          call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGConvection,rcollection)


          call lsyssc_scaleMatrix (rmatrixA,-1.0_DP)

          !call matio_writeMatrixHR(rmatrixMC,'./gmv/feat_routine.txt',.true.,0,'./gmv/feat_routine.txt','(E20.10)')



          ! And now to the entries of the matrix. For assembling of the entries,
          ! we need a bilinear form, which first has to be set up manually.
          ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
          ! scalar system matrix in 2D.

          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC

          ! In the standard case, we have constant coefficients:
          rform%ballCoeffConstant = .false.
          rform%BconstantCoeff = .false.

          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          rcollection%p_rvectorQuickAccess1 => rsolBlock
          call bilf_dg_buildMatrixScEdge2D (rform, CUB_G5_1D, .false., rmatrixA,&
               rsolBlock%Rvectorblock(1), raddTriaData,&
               flux_dg_implicitConvection_sim)!,&
          !rcollection, cconstrType)

          !call matio_writeMatrixHR(rmatrixMC,'./gmv/dg_routine.txt',.true.,0,'./gmv/dg_routine.txt','(E20.10)')


          ! Release solver data and structure
          call linsol_doneData (p_rsolverNode)
          call linsol_doneStructure (p_rsolverNode)

          ! Release the solver node and all subnodes attached to it (if at all):
          call linsol_releaseSolver (p_rsolverNode)

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
          !
          ! We 'know' already (from the problem definition) that we have four boundary
          ! segments in the domain. Each of these, we want to use for enforcing
          ! some kind of boundary condition.
          !
          ! We ask the bondary routines to create a 'boundary region' - which is
          ! simply a part of the boundary corresponding to a boundary segment.
          ! A boundary region roughly contains the type, the min/max parameter value
          ! and whether the endpoints are inside the region or not.
          call boundary_createRegion(rboundary,1,1,rboundaryRegion)

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
               getBoundaryValues_2D_rad)

          ! Now to the edge 2 of boundary component 1 the domain.
          call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
               rboundaryRegion,rdiscreteBC,&
               getBoundaryValues_2D_rad)

          !    ! Edge 3 of boundary component 1.
          !    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !    
          !    ! Edge 4 of boundary component 1. That is it.
          !    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_ones)

          ! Hang the pointer into the vector and matrix. That way, these
          ! boundary conditions are always connected to that matrix and that
          ! vector.
          rmatrixBlock%p_rdiscreteBC => rdiscreteBC
          rrhsBlock%p_rdiscreteBC => rdiscreteBC
          rsolBlock%p_rdiscreteBC => rdiscreteBC
          rtempBlock%p_rdiscreteBC => rdiscreteBC
          !                             
          !    ! Now we have block vectors for the RHS and the matrix. What we
          !    ! need additionally is a block vector for the solution and
          !    ! temporary data. Create them using the RHS as template.
          !    ! Fill the solution vector with 0:
          !    call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
          !    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)

          ! Next step is to implement boundary conditions into the RHS,
          ! solution and matrix. This is done using a vector/matrix filter
          ! for discrete boundary conditions.
          ! The discrete boundary conditions are already attached to the
          ! vectors/matrix. Call the appropriate vector/matrix filter that
          ! modifies the vectors/matrix according to the boundary conditions.

          !    call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
          !    pause

          call lsyssc_duplicateMatrix (rmatrixA, &
               rmatrixBlock%Rmatrixblock(1,1), &
               LSYSSC_DUP_SHARE, &
               !LSYSSC_DUP_EMPTY)
          LSYSSC_DUP_SHARE)



          call lsysbl_clearVector (rrhsBlock)

          call vecfil_discreteBCrhs (rrhsBlock)
          call vecfil_discreteBCsol (rsolBlock)
          call matfil_discreteBC (rmatrixBlock)

          !call matio_writeMatrixHR(rmatrixMC,'./gmv/bd_cond.txt',.true.,0,'./gmv/bd_cond.txt','(E20.10)')

          !    call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
          !    pause

          ! During the linear solver, the boundary conditions are also
          ! frequently imposed to the vectors. But as the linear solver
          ! does not work with the actual solution vectors but with
          ! defect vectors instead.
          ! So, set up a filter chain that filters the defect vector
          ! during the solution process to implement discrete boundary conditions.
          RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL


          ! Create a BiCGStab-solver.
          nullify(p_rpreconditioner)
          !call linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
          !call linsol_initILU0(p_rpreconditioner)

          !call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
          call linsol_initUMFPACK4 (p_rsolverNode)
          !call linsol_initGMRES (p_rsolverNode,50,p_rpreconditioner)

          ! Set the output level of the solver to 2 for some output
          p_rsolverNode%ioutputLevel = 2

          ! The linear solver stops, when this relative or absolut norm of
          ! the residual is reached.
          p_rsolverNode%depsRel = 1.0e-14
          p_rsolverNode%depsAbs = 1.0e-14

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
          if (ierror .ne. LINSOL_ERR_NOERROR) then
             write(*,*) 'linsol_initStructure',ierror
             pause
             stop
          end if
          call linsol_initData (p_rsolverNode, ierror)
          if (ierror .ne. LINSOL_ERR_NOERROR) then
             write(*,*) 'linsol_initData', ierror
             pause
             stop
          end if






          call linsol_solveAdaptively (p_rsolverNode,rsolBlock,rrhsBlock,rtempBlock)



          ! Nonlinear implicit Burgers equation
       case (6)

          ! Now as the discretisation is set up, we can start to generate
          ! the structure of the system matrix which is to solve.
          ! We create a scalar matrix, based on the discretisation structure
          ! for our one and only solution component.
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixBlock%Rmatrixblock(1,1),&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)



          ! Now as the discretisation is set up, we can start to generate
          ! the structure of the system matrix which is to solve.
          ! We create a scalar matrix, based on the discretisation structure
          ! for our one and only solution component.
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixMC,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)

          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC

          ! In the standard case, we have constant coefficients:
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.


          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          !call lsyssc_clearMatrix (rmatrixMC)
          ! We need to give the solution vector to the routine (by the collection)
          call bilf_buildMatrixScalar (rform,.true.,rmatrixMC)

          ! Now as the discretisation is set up, we can start to generate
          ! the structure of the system matrix which is to solve.
          ! We create a scalar matrix, based on the discretisation structure
          ! for our one and only solution component.
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixA,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)


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
          !
          ! We 'know' already (from the problem definition) that we have four boundary
          ! segments in the domain. Each of these, we want to use for enforcing
          ! some kind of boundary condition.
          !
          ! We ask the bondary routines to create a 'boundary region' - which is
          ! simply a part of the boundary corresponding to a boundary segment.
          ! A boundary region roughly contains the type, the min/max parameter value
          ! and whether the endpoints are inside the region or not.
          call boundary_createRegion(rboundary,1,1,rboundaryRegion)

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
               getBoundaryValues_2D_zeros)

          ! Now to the edge 2 of boundary component 1 the domain.
          call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
               rboundaryRegion,rdiscreteBC,&
               getBoundaryValues_2D_zeros)

          ! Edge 3 of boundary component 1.
          call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
               rboundaryRegion,rdiscreteBC,&
               getBoundaryValues_2D_zeros)

          ! Edge 4 of boundary component 1. That is it.
          call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
               rboundaryRegion,rdiscreteBC,&
               getBoundaryValues_2D_zeros)

          ! Hang the pointer into the vector and matrix. That way, these
          ! boundary conditions are always connected to that matrix and that
          ! vector.
          rmatrixBlock%p_rdiscreteBC => rdiscreteBC
          rrhsBlock%p_rdiscreteBC => rdiscreteBC
          rsolBlock%p_rdiscreteBC => rdiscreteBC
          rtempBlock%p_rdiscreteBC => rdiscreteBC

          ttime = 0.0_dp

          Burgerstimestepping: do

             ! Print the time level
             write(*,*) 'Time: ', ttime

             ! Calculate the right-hand-side vector b
             call lsyssc_scalarMatVec(rMatrixMC,rsolBlock%RvectorBlock(1),rrhsBlock%RvectorBlock(1),1.0_dp,0.0_dp)


             Burgers_defcorr: do

                ! Release solver data and structure
                call linsol_doneData (p_rsolverNode)
                call linsol_doneStructure (p_rsolverNode)

                ! Release the solver node and all subnodes attached to it (if at all):
                call linsol_releaseSolver (p_rsolverNode)

                ! Calculate the defect vector

                ! And now to the entries of the matrix. For assembling of the entries,
                ! we need a bilinear form, which first has to be set up manually.
                ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
                ! scalar system matrix in 2D.

                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_DERIV_X
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_DERIV_Y

                ! In the standard case, we have constant coefficients:
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.

                ! Now we can build the matrix entries.
                ! We specify the callback function coeff_Laplace for the coefficients.
                ! As long as we use constant coefficients, this routine is not used.
                ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
                ! the framework will call the callback routine to get analytical
                ! data.
                !call lsyssc_clearMatrix (rmatrixMC)
                ! We need to give the solution vector to the routine (by the collection)
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGBurgers,rcollection)


                call lsyssc_scaleMatrix (rmatrixA,-1.0_DP)

                !call matio_writeMatrixHR(rmatrixMC,'./gmv/feat_routine.txt',.true.,0,'./gmv/feat_routine.txt','(E20.10)')



                ! And now to the entries of the matrix. For assembling of the entries,
                ! we need a bilinear form, which first has to be set up manually.
                ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
                ! scalar system matrix in 2D.

                rform%itermCount = 1
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_FUNC

                ! In the standard case, we have constant coefficients:
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.

                ! Now we can build the matrix entries.
                ! We specify the callback function coeff_Laplace for the coefficients.
                ! As long as we use constant coefficients, this routine is not used.
                ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
                ! the framework will call the callback routine to get analytical
                ! data.
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                call bilf_dg_buildMatrixScEdge2D (rform, CUB_G5_1D, .false., rmatrixA,&
                     rsolBlock%Rvectorblock(1), raddTriaData,&
                     flux_dg_implicitBurgers_sim)!,&
                !rcollection, cconstrType)


                !        call lsyssc_clearMatrix (rmatrixA)

                ! Calculate the preconditioner matrix                                         
                call lsyssc_matrixLinearComb (rmatrixMC,rmatrixA,1.0_dp,dt,&
                     .true.,.true.,.true.,.false.,&
                     rmatrixBlock%Rmatrixblock(1,1))

                !        call matio_writeMatrixHR(rmatrixMC,'./gmv/MC.txt',.true.,0,'./gmv/MC.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixA,'./gmv/A.txt',.true.,0,'./gmv/A.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,1),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        pause

                ! Calculate the defect vector
                call lsyssc_copyVector (rrhsBlock%RvectorBlock(1),rdefBlock%RvectorBlock(1))
                call lsyssc_scalarMatVec(rmatrixBlock%Rmatrixblock(1,1),rsolBlock%RvectorBlock(1),&
                     rdefBlock%RvectorBlock(1),-1.0_dp,1.0_dp)

                !       ! * If different Preconditioners should be used *        
                !       ! Use MC as preconditioner                                     
                !        call lsyssc_matrixLinearComb (rmatrixMC,rmatrixA,1.0_dp,0.0_dp,&
                !                                       .true.,.true.,.true.,.false.,&
                !                                       rmatrixBlock%Rmatrixblock(1,1)) 

                !       ! Use MC+edge as preconditioner
                !        call bilf_dg_buildMatrixScEdge2D (rform, CUB_G5_1D, .true., rmatrixA,&
                !                                                 rsolBlock%Rvectorblock(1), raddTriaData,&
                !                                                 flux_dg_implicitBurgers_sim)!,&
                !                                                 !rcollection, cconstrType)                                     
                !        call lsyssc_matrixLinearComb (rmatrixMC,rmatrixA,1.0_dp,0.0_dp,&
                !                                       .true.,.true.,.true.,.false.,&
                !                                       rmatrixBlock%Rmatrixblock(1,1))

                !       ! Use MC+cell as preconditioner
                !        rform%itermCount = 2
                !        rform%Idescriptors(1,1) = DER_FUNC
                !        rform%Idescriptors(2,1) = DER_DERIV_X
                !        rform%Idescriptors(1,2) = DER_FUNC
                !        rform%Idescriptors(2,2) = DER_DERIV_Y
                !        rform%ballCoeffConstant = .false.
                !        rform%BconstantCoeff = .false.
                !        rcollection%p_rvectorQuickAccess1 => rsolBlock
                !        call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGBurgers,rcollection)                                    
                !        call lsyssc_matrixLinearComb (rmatrixMC,rmatrixA,1.0_dp,0.0_dp,&
                !                                       .true.,.true.,.true.,.false.,&
                !                                       rmatrixBlock%Rmatrixblock(1,1))




                ! During the linear solver, the boundary conditions are also
                ! frequently imposed to the vectors. But as the linear solver
                ! does not work with the actual solution vectors but with
                ! defect vectors instead.
                ! So, set up a filter chain that filters the defect vector
                ! during the solution process to implement discrete boundary conditions.
                RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL


                ! Fill in the boundary conditions
                call vecfil_discreteBCrhs (rrhsBlock)
                call vecfil_discreteBCsol (rsolBlock)
                call matfil_discreteBC (rmatrixBlock)

                ! Create a BiCGStab-solver.
                nullify(p_rpreconditioner)
                !    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
                call linsol_initUMFPACK4 (p_rsolverNode)
                !    call linsol_initJacobi (p_rsolverNode)

                ! Set the output level of the solver to 2 for some output
                p_rsolverNode%ioutputLevel = 0

                ! The linear solver stops, when this relative or absolut norm of
                ! the residual is reached.
                p_rsolverNode%depsRel = 1.0e-4
                p_rsolverNode%depsAbs = 1.0e-4

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
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initStructure',ierror
                   pause
                   stop
                end if
                call linsol_initData (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initData', ierror
                   pause
                   stop
                end if




                ! Solve for solution updade
                ! Afterwards the defect vector contains the solution update
                call linsol_precondDefect (p_RsolverNode,rdefBlock)


                ! Update the solution
                call lsyssc_vectorLinearComb (rDefBlock%Rvectorblock(1),rSolBlock%Rvectorblock(1),1.0_dp,1.0_dp)
                call vecfil_discreteBCsol (rsolBlock)

                write(*,*) '  Defect:',lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2)

                if (lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2) < 1.0e-10) exit Burgers_defcorr

                !        write(*,*) lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2)
                !        call dg2vtk(rDefBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause


             end do Burgers_defcorr

             ttime = ttime + dt
             if (ttime.ge.(ttfinal-1.0e-12)) exit Burgerstimestepping
          end do Burgerstimestepping








          ! Fully implicit Euler for linear system (decoupled)
       case (7)

          !    call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
          !    pause

          ! Calculate matrix MC
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixMC,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)

          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.
          !call lsyssc_clearMatrix (rmatrixMC)
          call bilf_buildMatrixScalar (rform,.true.,rmatrixMC)


          ! Set up the structure of the system matrix
!!! call lsysbl_createMatBlockByDiscr (rblockDiscretisationTrial,rmatrixBlock)
          call lsysbl_createEmptyMatrix (rmatrixBlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixBlock%Rmatrixblock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixBlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixBlock)


          ! Set up the structure of the matrix A
          call lsysbl_createEmptyMatrix (rmatrixABlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixABlock%RmatrixBlock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixABlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixABlock)

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
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Now to the edge 2 of boundary component 1 the domain.
          !    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Edge 3 of boundary component 1.
          !    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !    
          !    ! Edge 4 of boundary component 1. That is it.
          !    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Hang the pointer into the vector and matrix. That way, these
          !    ! boundary conditions are always connected to that matrix and that
          !    ! vector.


          call bcasm_initDiscreteBC(rdiscreteBC)
          call boundary_createRegion(rboundary,1,1,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
               rboundaryRegion,rdiscreteBC,&
               getBoundaryValues_2D_zeros)
          call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
               rboundaryRegion,rdiscreteBC,&
               getBoundaryValues_2D_hump)

          rmatrixBlock%p_rdiscreteBC => rdiscreteBC
          rrhsBlock%p_rdiscreteBC => rdiscreteBC
          rsolBlock%p_rdiscreteBC => rdiscreteBC
          rtempBlock%p_rdiscreteBC => rdiscreteBC

          ttime = 0.0_dp

          iEulertimestepping: do

             ! Print the time level
             write(*,*) 'Time: ', ttime

             ! Calculate the right-hand-side vector b
             do i = 1, nvar2d
                call lsyssc_scalarMatVec(rMatrixMC,rsolBlock%RvectorBlock(i),rrhsBlock%RvectorBlock(i),1.0_dp,0.0_dp)
             end do



             iEuler_defcorr: do

                ! Release solver data and structure
                call linsol_doneData (p_rsolverNode)
                call linsol_doneStructure (p_rsolverNode)

                ! Release the solver node and all subnodes attached to it (if at all):
                call linsol_releaseSolver (p_rsolverNode)

                ! Calculate the defect vector

                ! First calculate the matrix for the cell terms
                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_DERIV_X
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_DERIV_Y
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                !call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGBurgers,rcollection)
                if (ielementtype.ne.EL_DG_T0_2D) then
                   call bilf_buildMatrixBlock2 (rform, .true., rmatrixABlock,&
                        fcoeff_buildMatrixBl_sim_iSystem,rcollection)!,rscalarAssemblyInfo)
                else
                   call lsysbl_clearMatrix(rmatrixABlock)
                end if

                call lsysbl_scaleMatrix (rmatrixABlock,-1.0_DP)
                !call matio_writeMatrixHR(rmatrixMC,'./gmv/feat_routine.txt',.true.,0,'./gmv/feat_routine.txt','(E20.10)')



                ! Next calculate the edge terms
                rform%itermCount = 1
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_FUNC
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock   
                call bilf_dg_buildMatrixBlEdge2D (rform, CUB_G5_1D, .false., rmatrixABlock,&
                     rsolBlock, raddTriaData,&
                     flux_dg_buildMatrixBlEdge2D_sim_iSystem)!,&
                !rcollection, cconstrType)


                ! Calculate the preconditioner matrix
                !        call lsysbl_clearMatrix(rmatrixBlock)
                call lsysbl_copyMatrix (rmatrixABlock,rmatrixBlock)
                do i = 1, nvar2d
                   call lsyssc_matrixLinearComb (rmatrixMC,rmatrixABlock%RmatrixBlock(i,i),1.0_dp,dt,&
                        .true.,.true.,.true.,.false.,&
                        rmatrixBlock%Rmatrixblock(i,i))
                   !        call lsyssc_matrixLinearComb2 (rmatrixMC,1.0_dp,rmatrixABlock%RmatrixBlock(i,j),dt,rmatrixABlock%RmatrixBlock(i,j),&
                   !                                       .false.,.false.,.true.)
                end do

                !        call matio_writeMatrixHR(rmatrixMC,'./gmv/MC.txt',.true.,0,'./gmv/MC.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixA,'./gmv/A.txt',.true.,0,'./gmv/A.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,1),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        pause


                ! Calculate the defect vector
                call lsysbl_copyVector (rrhsBlock,rdefBlock)
                call lsysbl_blockMatVec(rmatrixBlock,rsolBlock,&
                     rdefBlock,-1.0_dp,1.0_dp)


                !        do i = 1, nvar2d
                !        do j = 1, nvar2d
                !          call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(i,j),'./gmv/Pblock.txt',.true.,0,'./gmv/Pblock.txt','(E20.10)')                       
                !          pause
                !        end do 
                !        end do


                ! During the linear solver, the boundary conditions are also
                ! frequently imposed to the vectors. But as the linear solver
                ! does not work with the actual solution vectors but with
                ! defect vectors instead.
                ! So, set up a filter chain that filters the defect vector
                ! during the solution process to implement discrete boundary conditions.
                RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL


                ! Fill in the boundary conditions
                call vecfil_discreteBCrhs (rrhsBlock)
                call vecfil_discreteBCsol (rsolBlock)
                call matfil_discreteBC (rmatrixBlock)

                ! Create a BiCGStab-solver.
                nullify(p_rpreconditioner)
                !    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
                call linsol_initUMFPACK4 (p_rsolverNode)
                !    call linsol_initJacobi (p_rsolverNode)

                ! Set the output level of the solver to 2 for some output
                p_rsolverNode%ioutputLevel = 0

                ! The linear solver stops, when this relative or absolut norm of
                ! the residual is reached.
                p_rsolverNode%depsRel = 1.0e-4
                p_rsolverNode%depsAbs = 1.0e-4

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
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initStructure',ierror
                   pause
                   stop
                end if
                call linsol_initData (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initData', ierror
                   pause
                   stop
                end if




                ! Solve for solution updade
                ! Afterwards the defect vector contains the solution update
                call linsol_precondDefect (p_RsolverNode,rdefBlock)


                ! Update the solution
                call lsysbl_vectorLinearComb (rDefBlock,rSolBlock,1.0_dp,1.0_dp)
                !        call vecfil_discreteBCsol (rsolBlock)


                call vecfil_discreteBCsol (rsolBlock)

                write(*,*) '  Defect:',lsysbl_vectorNorm(rDefBlock,LINALG_NORML2)



                ! Output solution to vtk file
                sofile = './gmv/u2d_rho'
                call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                sofile = './gmv/u2d_rhou'
                call dg2vtk(rsolBlock%Rvectorblock(2),iextraPoints,sofile,ifilenumber)
                pause



                if (lsysbl_vectorNorm(rDefBlock,LINALG_NORML2) < 1.0e-12) exit iEuler_defcorr

                !        write(*,*) lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2)
                !        call dg2vtk(rDefBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause


                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause

             end do iEuler_defcorr

             ttime = ttime + dt
             if (ttime.ge.(ttfinal-1.0e-12)) exit iEulertimestepping
          end do iEulertimestepping








          ! Fully implicit Euler equation
       case (8)

          !    call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
          !    pause

          ! Calculate matrix MC
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixMC,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)

          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.
          !call lsyssc_clearMatrix (rmatrixMC)
          call bilf_buildMatrixScalar (rform,.true.,rmatrixMC)


          ! Set up the structure of the system matrix
!!! call lsysbl_createMatBlockByDiscr (rblockDiscretisationTrial,rmatrixBlock)
          call lsysbl_createEmptyMatrix (rmatrixBlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixBlock%Rmatrixblock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixBlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixBlock)


          ! Set up the structure of the matrix A
          call lsysbl_createEmptyMatrix (rmatrixABlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixABlock%RmatrixBlock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixABlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixABlock)

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
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Now to the edge 2 of boundary component 1 the domain.
          !    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Edge 3 of boundary component 1.
          !    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !    
          !    ! Edge 4 of boundary component 1. That is it.
          !    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Hang the pointer into the vector and matrix. That way, these
          !    ! boundary conditions are always connected to that matrix and that
          !    ! vector.
          !    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
          !    rrhsBlock%p_rdiscreteBC => rdiscreteBC
          !    rsolBlock%p_rdiscreteBC => rdiscreteBC
          !    rtempBlock%p_rdiscreteBC => rdiscreteBC

          ttime = 0.0_dp

          iEulertimestepping2: do

             ! Print the time level
             write(*,*) 'Time: ', ttime

             ! Calculate the right-hand-side vector b
             do i = 1, nvar2d
                call lsyssc_scalarMatVec(rMatrixMC,rsolBlock%RvectorBlock(i),rrhsBlock%RvectorBlock(i),1.0_dp,0.0_dp)
             end do



             iEuler_defcorr2: do

                ! Release solver data and structure
                call linsol_doneData (p_rsolverNode)
                call linsol_doneStructure (p_rsolverNode)

                ! Release the solver node and all subnodes attached to it (if at all):
                call linsol_releaseSolver (p_rsolverNode)

                ! Calculate the defect vector

                ! First calculate the matrix for the cell terms
                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_DERIV_X
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_DERIV_Y
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                !call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGBurgers,rcollection)
                if (ielementtype.ne.EL_DG_T0_2D) then
                   call bilf_buildMatrixBlock2 (rform, .true., rmatrixABlock,&
                        fcoeff_buildMatrixBl_sim_iEuler,rcollection)!,rscalarAssemblyInfo)
                else
                   call lsysbl_clearMatrix(rmatrixABlock)
                end if

                call lsysbl_scaleMatrix (rmatrixABlock,-1.0_DP)
                !call matio_writeMatrixHR(rmatrixMC,'./gmv/feat_routine.txt',.true.,0,'./gmv/feat_routine.txt','(E20.10)')



                ! Next calculate the edge terms
                rform%itermCount = 1
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_FUNC
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock   
                rcollection%Dquickaccess(1) = dt
                call bilf_dg_buildMatrixBlEdge2D (rform, CUB_G5_1D, .false., rmatrixABlock,&
                     rsolBlock, raddTriaData,&
                     flux_dg_buildMatrixBlEdge2D_sim_iEuler,&
                     rcollection)!, cconstrType)

                ! Calculate the preconditioner matrix
                !        call lsysbl_clearMatrix(rmatrixBlock)
                call lsysbl_copyMatrix (rmatrixABlock,rmatrixBlock)
                do i = 1, nvar2d
                   call lsyssc_matrixLinearComb (rmatrixMC,rmatrixABlock%RmatrixBlock(i,i),1.0_dp,dt,&
                        .true.,.true.,.true.,.false.,&
                        rmatrixBlock%Rmatrixblock(i,i))
                   !        call lsyssc_matrixLinearComb2 (rmatrixMC,1.0_dp,rmatrixABlock%RmatrixBlock(i,j),dt,rmatrixABlock%RmatrixBlock(i,j),&
                   !                                       .false.,.false.,.true.)
                end do

                !        call matio_writeMatrixHR(rmatrixMC,'./gmv/MC.txt',.true.,0,'./gmv/MC.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixA,'./gmv/A.txt',.true.,0,'./gmv/A.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,1),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,3),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        pause


                ! Calculate the defect vector
                call lsysbl_copyVector (rrhsBlock,rdefBlock)
                call lsysbl_blockMatVec(rmatrixBlock,rsolBlock,&
                     rdefBlock,-1.0_dp,1.0_dp)


                !        do i = 1, nvar2d
                !        do j = 1, nvar2d
                !          call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(i,j),'./gmv/Pblock.txt',.true.,0,'./gmv/Pblock.txt','(E20.10)')                       
                !          pause
                !        end do 
                !        end do


                !        ! During the linear solver, the boundary conditions are also
                !    ! frequently imposed to the vectors. But as the linear solver
                !    ! does not work with the actual solution vectors but with
                !    ! defect vectors instead.
                !    ! So, set up a filter chain that filters the defect vector
                !    ! during the solution process to implement discrete boundary conditions.
                !    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
                !    
                !    
                !    ! Fill in the boundary conditions
                !        call vecfil_discreteBCrhs (rrhsBlock)
                !        call vecfil_discreteBCsol (rsolBlock)
                !        call matfil_discreteBC (rmatrixABlock)

                ! Create a BiCGStab-solver.
                nullify(p_rpreconditioner)
                !    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
                call linsol_initUMFPACK4 (p_rsolverNode)
                !    call linsol_initJacobi (p_rsolverNode)

                ! Set the output level of the solver to 2 for some output
                p_rsolverNode%ioutputLevel = 0

                ! The linear solver stops, when this relative or absolut norm of
                ! the residual is reached.
                p_rsolverNode%depsRel = 1.0e-4
                p_rsolverNode%depsAbs = 1.0e-4

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
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initStructure',ierror
                   pause
                   stop
                end if
                call linsol_initData (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initData', ierror
                   pause
                   stop
                end if




                ! Solve for solution updade
                ! Afterwards the defect vector contains the solution update
                call linsol_precondDefect (p_RsolverNode,rdefBlock)


                ! Update the solution
                call lsysbl_vectorLinearComb (rDefBlock,rSolBlock,1.0_dp,1.0_dp)
                !        call vecfil_discreteBCsol (rsolBlock)

                write(*,*) '  Defect:',lsysbl_vectorNorm(rDefBlock,LINALG_NORML2)

                if (lsysbl_vectorNorm(rDefBlock,LINALG_NORML2) < 1.0e-12) exit iEuler_defcorr2

                !        write(*,*) lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2)
                !        call dg2vtk(rDefBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause


                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause

                !        ! Output solution to vtk file
                !        sofile = './gmv/u2d_rho'
                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhou'
                !        call dg2vtk(rsolBlock%Rvectorblock(2),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhov'
                !        call dg2vtk(rsolBlock%Rvectorblock(3),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhoE'
                !        call dg2vtk(rsolBlock%Rvectorblock(4),iextraPoints,sofile,ifilenumber)
                !        pause


             end do iEuler_defcorr2




             ttime = ttime + dt
             if (ttime.ge.(ttfinal-1.0e-12)) exit iEulertimestepping2
          end do iEulertimestepping2





          ! Fully implicit linear convection, Newton-Like implementation
       case(9)

          ! Calculate matrix MC
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixMC,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)

          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.
          !call lsyssc_clearMatrix (rmatrixMC)
          call bilf_buildMatrixScalar (rform,.true.,rmatrixMC)


          ! Set up the structure of the system matrix
!!! call lsysbl_createMatBlockByDiscr (rblockDiscretisationTrial,rmatrixBlock)
          call lsysbl_createEmptyMatrix (rmatrixBlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixBlock%Rmatrixblock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixBlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixBlock)


          ! Set up the structure of the matrix A
          call lsysbl_createEmptyMatrix (rmatrixABlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixABlock%RmatrixBlock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixABlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixABlock)

          

          ttime = 0.0_dp

          iConvNewtontimestepping2: do

             ! Print the time level
             write(*,*) 'Time: ', ttime

             ! Create RHS-Vector

             ! At first set up the corresponding linear form (f,Phi_j):
             rlinformedge%itermCount = 1
             rlinformedge%Idescriptors(1) = DER_FUNC2D
             
             ! Now use the dg-function for the edge terms
             call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G5_1D, .true.,&
                  rrhsBlock,rsolBlock,&
                  raddTriaData,&
                  iConv_flux_dg_buildVectorBlEdge2D_sim_Newton,&
                  rcollection)


             call lsysbl_scaleVector (rrhsBlock,-1.0_DP)


             ! Then add the convection terms
             if(ielementType .ne. EL_DG_T0_2D) then
                ! Set up linear form
                rlinformconv%itermCount = 2
                rlinformconv%Idescriptors(1) = DER_DERIV_X
                rlinformconv%Idescriptors(2) = DER_DERIV_Y
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                ! Call the linearformroutine
                call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                     iConv_flux_sys_block,rcollection)
             end if



!             iEuler_defcorr3: do

                ! Release solver data and structure
                call linsol_doneData (p_rsolverNode)
                call linsol_doneStructure (p_rsolverNode)

                ! Release the solver node and all subnodes attached to it (if at all):
                call linsol_releaseSolver (p_rsolverNode)

                ! Calculate the defect vector

                ! First calculate the matrix for the cell terms
                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_DERIV_X
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_DERIV_Y
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                !call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGBurgers,rcollection)
                if (ielementtype.ne.EL_DG_T0_2D) then
                   call bilf_buildMatrixBlock2 (rform, .true., rmatrixABlock,&
                        fcoeff_buildMatrixBl_sim_iConvNewton,rcollection)!,rscalarAssemblyInfo)
                else
                   call lsysbl_clearMatrix(rmatrixABlock)
                end if

                call lsysbl_scaleMatrix (rmatrixABlock,-1.0_DP)
                !call matio_writeMatrixHR(rmatrixMC,'./gmv/feat_routine.txt',.true.,0,'./gmv/feat_routine.txt','(E20.10)')



                ! Next calculate the edge terms
                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_FUNC
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_FUNC
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock   
                rcollection%Dquickaccess(1) = dt
                call bilf_dg_buildMatrixBlEdge2D_ss (rform, CUB_G5_1D, .false., rmatrixABlock,&
                     rsolBlock, raddTriaData,&
                     flux_dg_buildMatrixBlEdge2D_sim_iConv_Newton_ss,&
                     rcollection)!, cconstrType)

                ! Calculate the preconditioner matrix
                !        call lsysbl_clearMatrix(rmatrixBlock)
                !call lsysbl_copyMatrix (rmatrixABlock,rmatrixBlock)
                call lsysbl_duplicateMatrix (rmatrixABlock,rmatrixBlock,&
                                     LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
                do i = 1, nvar2d

                   call lsyssc_matrixLinearComb (rmatrixMC,rmatrixABlock%RmatrixBlock(i,i),1.0_dp/dt,1.0_dp,&
                        .false.,.false.,.true.,.false.)

                end do

                !        call matio_writeMatrixHR(rmatrixMC,'./gmv/MC.txt',.true.,0,'./gmv/MC.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixA,'./gmv/A.txt',.true.,0,'./gmv/A.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,1),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,3),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        pause


!                ! Calculate the defect vector
!                call lsysbl_copyVector (rrhsBlock,rdefBlock)
!                call lsysbl_blockMatVec(rmatrixBlock,rsolBlock,&
!                     rdefBlock,-1.0_dp,1.0_dp)


                !        do i = 1, nvar2d
                !        do j = 1, nvar2d
                !          call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(i,j),'./gmv/Pblock.txt',.true.,0,'./gmv/Pblock.txt','(E20.10)')                       
                !          pause
                !        end do 
                !        end do


                !        ! During the linear solver, the boundary conditions are also
                !    ! frequently imposed to the vectors. But as the linear solver
                !    ! does not work with the actual solution vectors but with
                !    ! defect vectors instead.
                !    ! So, set up a filter chain that filters the defect vector
                !    ! during the solution process to implement discrete boundary conditions.
                !    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
                !    
                !    
                !    ! Fill in the boundary conditions
                !        call vecfil_discreteBCrhs (rrhsBlock)
                !        call vecfil_discreteBCsol (rsolBlock)
                !        call matfil_discreteBC (rmatrixABlock)

                ! Create a BiCGStab-solver.
                nullify(p_rpreconditioner)
                CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
                !call linsol_initJacobi (p_rpreconditioner)
                call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner)!,RfilterChain)
                !call linsol_initUMFPACK4 (p_rsolverNode)
                !call linsol_initGMRES (p_rsolverNode,50,p_rpreconditioner)

                ! Set the output level of the solver to 2 for some output
                p_rsolverNode%ioutputLevel = 2

                ! The linear solver stops, when this relative or absolut norm of
                ! the residual is reached.
                p_rsolverNode%depsRel = 1.0e-8
                p_rsolverNode%depsAbs = 1.0e-10

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
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initStructure',ierror
                   pause
                   stop
                end if
                call linsol_initData (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initData', ierror
                   pause
                   stop
                end if




!                ! Solve for solution updade
!                ! Afterwards the defect vector contains the solution update
!                call linsol_precondDefect (p_RsolverNode,rdefBlock)
                
                
                ! Solve for solution update
                call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)


                ! Update the solution
                call lsysbl_vectorLinearComb (rSolUpBlock,rSolBlock,1.0_dp,1.0_dp)
                !        call vecfil_discreteBCsol (rsolBlock)
                
                ! Limit solution
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)
                
!                ! Fix negative solution values
!                call dg_makeSolPos (rSolBlock)

                write(*,*) '  Update:',lsysbl_vectorNorm(rSolUpBlock,LINALG_NORML2)
!
!                if (lsysbl_vectorNorm(rDefBlock,LINALG_NORML2) < 1.0e-12) exit iEuler_defcorr3

                !        write(*,*) lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2)
                !        call dg2vtk(rDefBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause


                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause

                !        ! Output solution to vtk file
                !        sofile = './gmv/u2d_rho'
                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhou'
                !        call dg2vtk(rsolBlock%Rvectorblock(2),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhov'
                !        call dg2vtk(rsolBlock%Rvectorblock(3),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhoE'
                !        call dg2vtk(rsolBlock%Rvectorblock(4),iextraPoints,sofile,ifilenumber)
                !        pause


!             end do iEuler_defcorr3



!             call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror,&
!                  getReferenceFunction_2D)
!             call output_line ('L1-error: ' // sys_sdEL(derror,10) )

!             call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
!             call output_line ('Mass: ' // sys_sdEL(derror,10) )


             ttime = ttime + dt
             
             call lsysbl_releaseMatrix (rmatrixBlock)
             
             
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
             
             
             
             
             if (ttime.ge.(ttfinal-1.0e-12)) exit iConvNewtontimestepping2
          end do iConvNewtontimestepping2



          ! Implicit Euler equations Newton-like implementation
       case (10)

          ! Calculate matrix MC
          call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
               LSYSSC_MATRIX9,rmatrixMC,&
               rdiscretisation%RspatialDiscr(1),&
               BILF_MATC_EDGEBASED)

          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC
          rform%ballCoeffConstant = .true.
          rform%BconstantCoeff = .true.
          !call lsyssc_clearMatrix (rmatrixMC)
          call bilf_buildMatrixScalar (rform,.true.,rmatrixMC)


          ! Set up the structure of the system matrix
!!! call lsysbl_createMatBlockByDiscr (rblockDiscretisationTrial,rmatrixBlock)
          call lsysbl_createEmptyMatrix (rmatrixBlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixBlock%Rmatrixblock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixBlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixBlock)


          ! Set up the structure of the matrix A
          call lsysbl_createEmptyMatrix (rmatrixABlock,nvar2d,nvar2d)
          do i = 1, nvar2d
             do j = 1, nvar2d
                !    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                !                                     LSYSSC_MATRIX9,rmatrixABlock%RmatrixBlock(i,j),&
                !                                     rdiscretisation%RspatialDiscr(1),&
                !                                     BILF_MATC_EDGEBASED)
                call lsyssc_copyMatrix (rmatrixMC,rmatrixABlock%RmatrixBlock(i,j))
             end do
          end do
          call lsysbl_updateMatStrucInfo (rmatrixABlock)

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
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Now to the edge 2 of boundary component 1 the domain.
          !    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Edge 3 of boundary component 1.
          !    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !    
          !    ! Edge 4 of boundary component 1. That is it.
          !    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          !    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
          !                                       rboundaryRegion,rdiscreteBC,&
          !                                       getBoundaryValues_2D_zeros)
          !                             
          !    ! Hang the pointer into the vector and matrix. That way, these
          !    ! boundary conditions are always connected to that matrix and that
          !    ! vector.
          !    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
          !    rrhsBlock%p_rdiscreteBC => rdiscreteBC
          !    rsolBlock%p_rdiscreteBC => rdiscreteBC
          !    rtempBlock%p_rdiscreteBC => rdiscreteBC

          ttime = 0.0_dp

          iEulertimestepping3: do

             ! Print the time level
             write(*,*) 'Time: ', ttime

             ! Create RHS-Vector

             ! At first set up the corresponding linear form (f,Phi_j):
             rlinformedge%itermCount = 1
             rlinformedge%Idescriptors(1) = DER_FUNC2D
             
             ! Now use the dg-function for the edge terms
             call linf_dg_buildVectorBlockEdge2d (rlinformedge, CUB_G5_1D, .true.,&
                  rrhsBlock,rsolBlock,&
                  raddTriaData,&
                  Euler_flux_dg_buildVectorBlEdge2D_sim_Newton,&
                  rcollection)


             call lsysbl_scaleVector (rrhsBlock,-1.0_DP)


             ! Then add the convection terms
             if(ielementType .ne. EL_DG_T0_2D) then
                ! Set up linear form
                rlinformconv%itermCount = 2
                rlinformconv%Idescriptors(1) = DER_DERIV_X
                rlinformconv%Idescriptors(2) = DER_DERIV_Y
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                ! Call the linearformroutine
                call linf_buildVectorBlock2 (rlinformconv, .false., rrhsBlock,&
                     Euler_flux_sys_block,rcollection)
             end if



!             iEuler_defcorr3: do

                ! Release solver data and structure
                call linsol_doneData (p_rsolverNode)
                call linsol_doneStructure (p_rsolverNode)

                ! Release the solver node and all subnodes attached to it (if at all):
                call linsol_releaseSolver (p_rsolverNode)

                ! Calculate the defect vector

                ! First calculate the matrix for the cell terms
                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_DERIV_X
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_DERIV_Y
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock
                !call bilf_buildMatrixScalar (rform,.true.,rmatrixA,coeff_implicitDGBurgers,rcollection)
                if (ielementtype.ne.EL_DG_T0_2D) then
                   call bilf_buildMatrixBlock2 (rform, .true., rmatrixABlock,&
                        fcoeff_buildMatrixBl_sim_iEuler,rcollection)!,rscalarAssemblyInfo)
                else
                   call lsysbl_clearMatrix(rmatrixABlock)
                end if

                call lsysbl_scaleMatrix (rmatrixABlock,-1.0_DP)
                !call matio_writeMatrixHR(rmatrixMC,'./gmv/feat_routine.txt',.true.,0,'./gmv/feat_routine.txt','(E20.10)')



                ! Next calculate the edge terms
                rform%itermCount = 2
                rform%Idescriptors(1,1) = DER_FUNC
                rform%Idescriptors(2,1) = DER_FUNC
                rform%Idescriptors(1,2) = DER_FUNC
                rform%Idescriptors(2,2) = DER_FUNC
                rform%ballCoeffConstant = .false.
                rform%BconstantCoeff = .false.
                rcollection%p_rvectorQuickAccess1 => rsolBlock   
                rcollection%Dquickaccess(1) = dt
                call bilf_dg_buildMatrixBlEdge2D_ss (rform, CUB_G5_1D, .false., rmatrixABlock,&
                     rsolBlock, raddTriaData,&
                     flux_dg_buildMatrixBlEdge2D_sim_iEuler_Newton_ss,&
                     rcollection)!, cconstrType)

                ! Calculate the preconditioner matrix
                !        call lsysbl_clearMatrix(rmatrixBlock)
                !call lsysbl_copyMatrix (rmatrixABlock,rmatrixBlock)
!                call lsysbl_duplicateMatrix (rmatrixABlock,rmatrixBlock,&
!                                     LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

               
                do i = 1, nvar2d

                   call lsyssc_matrixLinearComb (rmatrixMC,rmatrixABlock%RmatrixBlock(i,i),1.0_dp/dt,1.0_dp,&
                        .false.,.false.,.true.,.false.)
                        
                end do


                !        call matio_writeMatrixHR(rmatrixMC,'./gmv/MC.txt',.true.,0,'./gmv/MC.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixA,'./gmv/A.txt',.true.,0,'./gmv/A.txt','(E20.10)')
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,1),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(1,3),'./gmv/P.txt',.true.,0,'./gmv/P.txt','(E20.10)')                       
                !        pause


!                ! Calculate the defect vector
!                call lsysbl_copyVector (rrhsBlock,rdefBlock)
!                call lsysbl_blockMatVec(rmatrixBlock,rsolBlock,&
!                     rdefBlock,-1.0_dp,1.0_dp)


                !        do i = 1, nvar2d
                !        do j = 1, nvar2d
                !          call matio_writeMatrixHR(rmatrixBlock%Rmatrixblock(i,j),'./gmv/Pblock.txt',.true.,0,'./gmv/Pblock.txt','(E20.10)')                       
                !          pause
                !        end do 
                !        end do


                !        ! During the linear solver, the boundary conditions are also
                !    ! frequently imposed to the vectors. But as the linear solver
                !    ! does not work with the actual solution vectors but with
                !    ! defect vectors instead.
                !    ! So, set up a filter chain that filters the defect vector
                !    ! during the solution process to implement discrete boundary conditions.
                !    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
                !    
                !    
                !    ! Fill in the boundary conditions
                !        call vecfil_discreteBCrhs (rrhsBlock)
                !        call vecfil_discreteBCsol (rsolBlock)
                !        call matfil_discreteBC (rmatrixABlock)

                ! Create a BiCGStab-solver.
                nullify(p_rpreconditioner)
                !call linsol_initGMRES (p_rpreconditioner,1)
                !call linsol_initBiCGStab (p_rpreconditioner)
                !call linsol_initJacobi (p_rpreconditioner)
                !call linsol_initSOR (p_rpreconditioner, 0.15_dp)
                !call linsol_initILU0(p_rpreconditioner)
                call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner)
                !call linsol_initGMRES (p_rsolverNode,10,p_rpreconditioner)
                !call linsol_initUMFPACK4 (p_rsolverNode)
                

                ! Set the output level of the solver to 2 for some output
                p_rsolverNode%ioutputLevel = 2

                ! The linear solver stops, when this relative or absolut norm of
                ! the residual is reached.
                p_rsolverNode%depsRel = 1.0e-4
                p_rsolverNode%depsAbs = 1.0e-6
                
                ! Set the maximum number of iteratons
                p_rsolverNode%nmaxIterations = 5000
                
                ! Attach the system matrix to the solver.
                ! First create an array with the matrix data (on all levels, but we
                ! only have one level here), then call the initialisation 
                ! routine to attach all these matrices.
                ! Remark: Do not make a call like
                !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
                ! This does not work on all compilers, since the compiler would have
                ! to create a temp array on the stack - which does not always work!
                Rmatrices = (/rmatrixABlock/)
                call linsol_setMatrices(p_RsolverNode,Rmatrices)

                ! Initialise structure/data of the solver. This allows the
                ! solver to allocate memory / perform some precalculation
                ! to the problem.
                call linsol_initStructure (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initStructure',ierror
                   pause
                   stop
                end if
                call linsol_initData (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) then
                   write(*,*) 'linsol_initData', ierror
                   pause
                   stop
                end if




!                ! Solve for solution updade
!                ! Afterwards the defect vector contains the solution update
!                call linsol_precondDefect (p_RsolverNode,rdefBlock)
                
                
                ! Solve for solution update
                call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
                
                iiterations=iiterations+p_rsolverNode%iiterations
                inumSolverCalls = inumSolverCalls + 1


                ! Update the solution
                call lsysbl_vectorLinearComb (rSolUpBlock,rSolBlock,1.0_dp,1.0_dp)
                !        call vecfil_discreteBCsol (rsolBlock)
                
                ! Limit solution
                if (ilimiter .eq. 12) call dg_realKuzmin (rsolBlock, raddTriaData)
                
!                ! Fix negative solution values
!                call dg_makeSolPos (rSolBlock)

                write(*,*) '  Defect:',lsysbl_vectorNorm(rSolUpBlock,LINALG_NORML2)
!
!                if (lsysbl_vectorNorm(rDefBlock,LINALG_NORML2) < 1.0e-12) exit iEuler_defcorr3

                !        write(*,*) lsyssc_vectorNorm(rDefBlock%Rvectorblock(1),LINALG_NORML2)
                !        call dg2vtk(rDefBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause


                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        pause

                !        ! Output solution to vtk file
                !        sofile = './gmv/u2d_rho'
                !        call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhou'
                !        call dg2vtk(rsolBlock%Rvectorblock(2),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhov'
                !        call dg2vtk(rsolBlock%Rvectorblock(3),iextraPoints,sofile,ifilenumber)
                !        sofile = './gmv/u2d_rhoE'
                !        call dg2vtk(rsolBlock%Rvectorblock(4),iextraPoints,sofile,ifilenumber)
                !        pause


!             end do iEuler_defcorr3



!             call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror,&
!                  getReferenceFunction_2D)
!             call output_line ('L1-error: ' // sys_sdEL(derror,10) )

!             call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror)
!             call output_line ('Mass: ' // sys_sdEL(derror,10) )


             ttime = ttime + dt
             
             call lsysbl_releaseMatrix (rmatrixBlock)
             
             
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
             
             
             
             
             if (ttime.ge.(ttfinal-1.0e-12)) exit iEulertimestepping3
          end do iEulertimestepping3


     


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

    !    sofile = 'l6'
    !    call loadSolutionData(rsolBlock%Rvectorblock(1),sofile)
    !    sofile = './gmv/u2d' 

    !    !select case (ioutputtype)
    !    !  case (1)
    ! Output solution to gmv file
    call dg2gmv(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
    !    !  case (2)
    ! Output solution to vtk file
    call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
    !    !end select





    ! Output solution to vtk file
    sofile = './gmv/u2d_rho'
    call dg2vtk(rsolBlock%Rvectorblock(1),iextraPoints,sofile,ifilenumber)
!    sofile = './gmv/u2d_rhou'
!    call dg2vtk(rsolBlock%Rvectorblock(2),iextraPoints,sofile,ifilenumber)
!    sofile = './gmv/u2d_rhov'
!    call dg2vtk(rsolBlock%Rvectorblock(3),iextraPoints,sofile,ifilenumber)
!    sofile = './gmv/u2d_rhoE'
!    call dg2vtk(rsolBlock%Rvectorblock(4),iextraPoints,sofile,ifilenumber)



    !    write(*,*) 'Writing steady solution to file'
    !    ! And output the steady projection
    !     call dg_proj2steady(rsolBlock,rtriangulation, rboundary)

    !    ! Saving the solution DOFs to file
    !    write(*,*) 'Writing DOFs to file'
    !    call saveSolutionData(rsolBlock%Rvectorblock(1),sofile,ifilenumber)


    !    sofile = 'l9'
    !    call loadSolutionData(rsolBlock%Rvectorblock(1),sofile)


    ! Calculate the error to the reference function.
    rsolBlock%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(1)%ccubTypeEval=CUB_G6_2d
    call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L1ERROR,derror,&
         getReferenceFunction_2D)
    call output_line ('L1-error: ' // sys_sdEL(derror,10) )

    ! Calculate the error to the reference function.
    call pperr_scalar (rsolBlock%Rvectorblock(1),PPERR_L2ERROR,derror,&
         getReferenceFunction_2D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )   
    
    
    
!    call pperr_scalarDefault (PPERR_L1ERROR, derror, rsolBlock%Rvectorblock(1),&
!                           getReferenceFunction_2D, ffunctionWeight = getWeightFunction) 

!    call dg_pperr_scalar2d_conf (PPERR_L1ERROR, derror, rsolBlock%rvectorBlock(1)%p_rspatialDiscr,&
!                                  rsolBlock%Rvectorblock(1), getReferenceFunction_2D)
!    call output_line ('L1-error: ' // sys_sdEL(derror,10) )
!    
!    call dg_pperr_scalar2d_conf (PPERR_L2ERROR, derror, rsolBlock%rvectorBlock(1)%p_rspatialDiscr,&
!                                  rsolBlock%Rvectorblock(1), getReferenceFunction_2D)
!    call output_line ('L2-error: ' // sys_sdEL(derror,10) )
    
    write(*,*) 'Average number of iterations:', real(iiterations)/real(inumSolverCalls)


    !    call calc_error(rsolBlock%Rvectorblock(1), derror, raddtriadata)
    !    write(*,*) derror


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

    call lsysbl_releaseVector (redgeBlock)
    call lsysbl_releaseVector (rconvBlock)
    call lsysbl_releaseVector (rsolLimiterBlock)
    call lsysbl_releaseVector (rsolOldBlock)
    call lsysbl_releaseVector (rk0)
    call lsysbl_releaseVector (rk1)
    call lsysbl_releaseVector (rk2)
    call lsysbl_releaseVector (rk3)
    call lsysbl_releaseVector (rdefBlock)
    call lsysbl_releaseVector (rimf1)
    call lsysbl_releaseVector (rimf2)
    call lsysbl_releaseVector (rsolSaveBlock)
    call lsysbl_releaseVector (rsourceTermBlock)

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

  end subroutine dg2d_sys

end module dg2d_systems

!##############################################################################
!# ****************************************************************************
!# <name> chemotaxis_blowup </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This is a solver for the problem proposed by Aida, Tsujikawa, Efendiev, Yagi and Mimura in their paper
!# "LOWER ESTIMATE OF THE ATTRACTOR DIMENSION FOR A CHEMOTAXIS GROWTH SYSTEM"
!#
!#	d/dt u =  0.0625 * \Delta u - \nu \nabla * ( u * \nabla c ) + u^2 ( 1-u )
!#      d/dt c = \Delta c - 32 * c + u
!#
!#
!# 	with a proper NBC for u and c as well as an IC
!#      CHI may depend on u and c , e.g. CHI=CHI(u,c)
!#      for a nonlinear CHI a simple defect correction is proposed
!#      for this case CHI should be defined in the underlying callback file chemotaxis_callback.f90
!#
!# The equation is discretised in time by an implicit Euler Method, providing the
!# following discrete equation:
!#
!#
!#	!!!! imlicit Euler !!!!!
!#
!#	1/dt * ( u_{n+1}-u{n} ) = Laplace u_{n+1} - grad (CHI*u_{n+1} * grad c_{n+1})	/ multiplying with test func and dt, take int
!#	1/dt * ( c_{n+1}-c{n} ) = Laplace c_{n+1} - c_{n+1} + u_{n}
!#	__________________
!#
!#	(u_{n+1},phi) 		=  dt*( Laplace u_{n+1},phi ) - ( grad (CHI*u_{n+1} * grad c_{n+1} ),phi ) + ( u_{n},phi )
!#
!#      (c_{n+1},phi) 		=  dt*( Laplace c_{n+1},phi ) - dt*( c_{n+1} +u_{n},phi ) + ( c_{n},phi )
!#	__________________
!#
!#	(u_{n+1},phi) 		= - dt*( grad u_{n+1},grad phi ) + dt*( CHI * u_{n+1}* grad_c, grad_phi ) + ( u_{n},phi )
!#
!#      (c_{n+1},phi) 		= - dt*( grad c_{n+1},grad phi ) - dt*( u_{n} + c_{n+1},phi ) + ( c_{n},phi )
!#	__________________
!#
!#	[M + dt*L - dt*M_2] u_{n+1} 	= [ M ] u_{n}
!#
!#      [M + dt*L + dt*M] c_{n+1} 		= [ M ] c_{n} - dt* M u_n
!#
!#
!# The whole solution process is divided into serveral subroutines, using
!# the BiCG-solver (defect correction ) for the linear (nonlinear) subproblems in every
!# timestep
!#
!# The module bases on the standard heatconduction example and shows which
!# modifications are necessary to create a chemotaxis solver from a heatconduction
!# solver. Boundary conditions, matrices and vectors are not all reassembled in
!# every timestep, in contrast to the heatcond_method1.f90
!# </purpose>
!##############################################################################

module chemotaxis_blowup
    use fsystem
    use genoutput
    use storage
    use linearsolver
    use boundary
    use trilinearformevaluation
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
    use genoutput
    use analyticprojection
    use matrixio
    use vectorio
    use collection
    use paramlist
    use linearalgebra
    use analyticprojection
    use meshregion
    use io
    use error
    use myTestModule

    use chemotaxis_callback_blowup

    implicit none

    contains

    !##################################################################################
    !<function>
    function Fsol(x,y,a,b) result (f_result1)					!----> fct for the initial condition
	    implicit none								!	of the solutionvector
	    real(DP) :: f_result1, x, y, a, b					!	(centered)

	    f_result1 = a*dexp(-b*((x-0.5_DP)*(x-0.5_DP)+(y-0.5_DP)*(y-0.5_DP)))
    end function Fsol
    !</function>

   !<function>
    function Fchemo(x,y,a,b) result (f_result2)					!----> fct for the initial condition
	    implicit none								!	of the chemoattractant
	    real(DP) :: f_result2, x, y, a, b					!	(centered)

	    f_result2 = a*dexp(-b*((x-0.5_DP)*(x-0.5_DP)+(y-0.5_DP)*(y-0.5_DP)))
    end function Fchemo
    !</function>
  !##################################################################################

  
    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemotaxispatternFCT subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine chemotaxisBlowup
    
    !<description>
    ! The routine performs the following tasks:
    !
    ! 1.) Read in parametrisation
    ! 2.) Read in triangulation
    ! 3.) Set up RHS
    ! 4.) Set up matrix
    ! 5.) Create solver structure
    ! 6.) Solve the problem
    ! 7.) Write solution to GMV file if desired
    ! 8.) Release all variables, finish
    !</description>

 
        ! Definitions of variables.
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

        ! A trilinear , bilinear and linear form describing the analytic problem to solve
        type(t_trilinearForm) :: rtriform
        type(t_bilinearForm) :: rform
        type(t_linearForm) :: rlinform

        ! A scalar matrix and vector. The vector accepts the RHS of the problem in scalar form.
        type(t_matrixScalar)  :: rmatrix, rmatrixchemo, rsysmatrix , rlaplace, rtemp, rtempLump
        type(t_matrixScalar)  :: rmassmatrix
        type(t_matrixScalar)  :: rmatrGradX, rmatrGradY
        type(t_vectorScalar)  :: rrhscell,rcell
        type(t_vectorScalar)  :: rchemo_x,rchemo_y
        type(t_vectorScalar)  :: rrhschemo,rchemoattract,cbvector_temp, &
                                    rinitcells, rinitchemo, rdef, rfc, rfu
        type(t_vectorScalar)  :: rrhschemo_notCoupled
        type(t_vectorScalar) :: ranalyticcells , ranalyticchemo, rcold, ruold

        !A pointer to the entries of vector rchemoattract (added to access the chemoattractant)
        real(DP), dimension(:), pointer ::  p_vectordata, p_uold, p_cold,&
                                                            p_analyticcells, p_analyticchemo, p_chemodata

        ! A block matrix and a couple of block vectors. These will be filled with data for the linear solver.
        type(t_matrixBlock) :: rmatrixBlock, rmatrixBlockchemo
        type(t_vectorBlock) :: rcellBlock,rvectorBlockchemo,rrhsBlock,rrhsBlockchemo,rtempBlock, rdefBlock, rfuBlock

        ! A set of variables describing the analytic and discrete boundary conditions.
        type(t_boundaryRegion) :: rboundaryRegion
        type(t_discreteBC), target :: rdiscreteBC, rdiscreteBCchemo

        ! A solver node that accepts parameters for the linear solver
        type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner,p_rsolverNode_cells,p_rpreconditioner_cells

        ! An array for the system matrix(matrices) during the initialisation of the linear solver.
        type(t_matrixBlock), dimension(1) :: Rmatrices

        ! A filter chain that describes how to filter the matrix/vector
        ! before/during the solution process. The filters usually implement
        ! boundary conditions.
        type(t_filterChain), dimension(1), target :: RfilterChain
        type(t_filterChain), dimension(:), pointer :: p_RfilterChain

        ! NLMAX receives the level where we want to solve.
        integer :: NLMAX

        ! Error indicator during initialisation of the solver
        integer :: ierror

        ! Output block for UCD output to GMV file
        type(t_ucdExport) :: rexport
        real(DP), dimension(:), pointer :: p_Ddata

        ! Time step size, number of timesteps.(also for a fixed amount of steps, since there is
        ! a basic loop control implemented dealing w/ the steady state analysis)
        real(DP) :: dtstep
        integer :: ntimesteps, steps

        ! Auxiliary variables to determine the digits of ntimesteps
        integer :: auxStepsDigits
        real(DP) :: auxStepsDigitsReal
        character(15) :: auxStepsStr
        integer :: istartIteration

        ! Time and time step counter
        real(DP) :: dtime
        integer :: itimestep

        ! coupled step counter
        integer :: icouplestep

        ! maximal iterations for the defect correction . Since we allow possible nonlinearities
        integer :: maxiterationdef

        ! Variables for the defect correction (nonlinear loop)
        real(DP) :: defect , defectTol , negthres

        ! A pointer to the handle of a double precision array
        real(DP), dimension(:,:), pointer ::  p_DvertexCoords
        real(DP), dimension(:), pointer ::  p_cbvector_temp, p_initcells, p_initchemo

        ! Integer for loops
        integer :: i , k

        ! A collection used for the callbackfct to transfer data as additional parameter
        type (t_collection) :: chemocollection

        ! Some params used to adjust the model
        real(DP) :: CHI, D_1, D_2, A_CHEMO , A_CELLS , B_CHEMO , B_CELLS, W, SIGMA, BETA, ALPHA, R,  PHI, GAMMA, N

        ! Defining the output level
        ! =0 means no gmv files generated
        ! =1 means gmv files printed
        ! =2 means just the gmv file is printed
        ! =3 means just the ISTEP_GMV-th gmv file are printed
        integer :: OUTPUT, ISTEP_GMV

        ! Defining the gmvoutputfolder
        integer :: gmvfolder, checkneg, iResume
                    
        ! error analysis
        real(DP) :: uerror, cerror, tol, Derr_chemoL2, Derr_cellsL2, Derr_chemoH1, Derr_cellsH1

        ! error-norm constant
        ! To be set in the .dat file
        integer :: CTRLNORM
                
        ! Whether or not taking an const initial sol
        ! =0 means const initial sol
        ! =1 means exponential initial sol
        ! Whether or not we set the IC by l_2 projection ( = 1 ) or by interpolation ( =0 )
        integer :: INITSOL , L2PROJ

        ! If INITSOL = 0 , meaning const init sol, we use these variables to define the const value
        real(DP) :: C_0, U_0

        ! Scaling params for the analytic given sols
        real(DP) :: SCALE_U , SCALE_C

        ! Some variables for time stats
        real(DP):: time_min, time_max, time_accum, time_start, time_stop

        ! Some variables for iterationstats
        integer :: iteration_c_max, iteration_c_min, iteration_u_max,&
                        iteration_u_min, iteration_defcorr_max,&
                        iteration_defcorr_min
        real(DP) :: iteration_c_average, iteration_u_average, iteration_defcorr_average

        ! To detremine negative values of the solution vectors
        logical :: quit

        ! A parameter list structure accepting the parameters from the DAT file.
        type(t_parlist) :: rparams

        ! A data collection to specify the BCs lateron
        type(t_collection) :: rcollection

       ! relaxation for the convective/chemotactical term
       real(DP) ::  convecRelaxation
 
       ! for comparision with analytical solution
       type(t_vectorScalar) :: ranalytChemo, ranalytCell, rtempAnalyt
       real(DP) :: dnorm

        ! An array for saving the stats about the steady state
        real(DP), dimension(:,:), allocatable :: steady_state_array, absolute_error
    
        ! index for writing file
        integer :: iwrite, iunit
    
        ! where to write the steady-state stats
        character(LEN=30) ::  stat_ouput

        ! Whether or not a prefix is used for writing the gmv-file name
        character(LEN=10) ::  gmv_prefix
        integer :: allow_prefix
   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! Ok, let's start. !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,">>start of the subroutine chemotaxispatternFCT ... "
        call parlst_init(rparams)

        call parlst_readfromfile(rparams, 'data/chemotaxis.dat')

        ! Getting some general params for the pb
        call parlst_getvalue_int (rparams, 'GENERAL', 'NLMAX', NLMAX, 7)
        call parlst_getvalue_int (rparams, 'GENERAL' , 'OUTPUT' , OUTPUT , 0)
        call parlst_getvalue_int (rparams, 'GENERAL' , 'ISTEP_GMV' , ISTEP_GMV , 1)
        call parlst_getvalue_int (rparams, 'GENERAL' , 'INITSOL' , INITSOL , 1)
        call parlst_getvalue_int ( rparams, 'GENERAL' , 'L2PROJ' , L2PROJ , 0)
        call parlst_getvalue_int(rparams, 'GENERAL', 'GMVFOLDER', gmvfolder, 0)
        call parlst_getvalue_int(rparams, 'GENERAL', 'CHECKNEG', checkneg , 0)
        call parlst_getvalue_double(rparams, 'GENERAL', 'DEFECTTOL', defectTol , 0.0000001_DP)
        call parlst_getvalue_double(rparams, 'GENERAL', 'NEGTHRES', negthres , -0.1_DP)
        call parlst_getvalue_string(rparams, 'GENERAL', 'STAT_OUTPUT', stat_ouput, 'steady_state' )
        call parlst_getvalue_int(rparams, 'GENERAL', 'iResume', iResume , 0)
        call parlst_getvalue_int(rparams, 'GENERAL', 'ISTARTIT', istartIteration , 0)
        call parlst_getvalue_string(rparams, 'GENERAL', 'GMV_PREFIX', gmv_prefix, ' ' )
        call parlst_getvalue_int(rparams, 'GENERAL', 'ALLOW_PREFIX', allow_prefix , 0)
            
        ! We set the chemotaxis and diffusion coefficient
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'CHI', CHI, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'W', W, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'BETA', BETA, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'SIGMA', SIGMA, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'ALPHA', ALPHA, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'R', R, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'D_1', D_1, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'PHI', PHI, 1.0_DP)
        PHI = BETA
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'GAMMA', GAMMA, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'N', N, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'D_2', D_2, 1.0_DP)
        call parlst_getvalue_double (rparams, 'coefficients', 'A_CHEMO', A_CHEMO, 10.0_dp)
        call parlst_getvalue_double (rparams, 'coefficients', 'A_CELLS', A_CELLS, 0.6_dp)
        call parlst_getvalue_double (rparams, 'coefficients', 'B_CHEMO', B_CHEMO, 10.0_dp)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'B_CELLS', B_CELLS, 10.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'U_0', U_0, 2.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'C_0', C_0, 1.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'SCALE_U', SCALE_U, 2.0_DP)
        call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'SCALE_C', SCALE_C, 1.0_DP)

        ! Initialise time step size and number of timesteps
        call parlst_getvalue_int (rparams, 'TIMESTEPPING', 'NTIMESTEPS', ntimesteps, 100)
        call parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTSTEP', dtstep, 0.00001_DP)
        call parlst_getvalue_int (rparams, 'TIMESTEPPING', 'STEPS', steps, 0)
        call parlst_getvalue_int (rparams, 'TIMESTEPPING', 'MAXITERATIONDEF', maxiterationdef, 10)

        ! We start at a certain time ...
        call parlst_getvalue_double(rparams, 'TIMESTEPPING' , 'STARTTIME' , dtime , 0.0_DP)

        ! Get the errorctrl norm
        call parlst_getvalue_int (rparams, 'NORM', 'CONTROLNORM', CTRLNORM, 2)

        ! Get the relaxation for the convective/chemoattractive term
        call parlst_getvalue_double (rparams, 'TEST', 'CONVECRELAXATION', convecRelaxation, 1.0_DP)
        convecRelaxation = alpha
        ! Get the tolerance threshold value for steady state issues
        call parlst_getvalue_double (rparams, 'ERROR', 'TOL', tol, 0.0001_DP)


        ! Get the path $PREDIR from the environment, where to read .prm/.tri files from.
        ! If that does not exist, write to the directory "./pre".
        if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

        ! Reset the prefix, if we don't use it
        if ( allow_prefix == 0 ) then
            gmv_prefix = ''
        end if

       ! At first, read in the parametrisation of the boundary and save
       ! it to rboundary.
       call boundary_read_prm(rboundary, './pre/cube1.prm')
        
       ! Now read in the basic triangulation.
       call tria_readTriFile2D (rtriangulation, './pre/cube1.tri', rboundary)
    
       ! Refine it.
       call tria_quickRefine2LevelOrdering (NLMAX,rtriangulation,rboundary)
    
       ! And create information about adjacencies and everything one needs from
       ! a triangulation.
       call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

        ! Now we can start to initialise the discretisation. At first, set up
        ! a block discretisation structure that specifies the blocks in the
        ! solution vector. In this simple problem, we only have one block.
        call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                   rtriangulation, rboundary)

        ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
        ! structures for every component of the solution vector.
        ! Initialise the first element of the list to specify the element
        ! and cubature rule for this solution component:
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_E011,CUB_G2X2,rtriangulation, rboundary)
                
        ! Now as the discretisation is set up, we can start to generate the structure of the system matrix.
        ! Get some pointers  for the errorctrl
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),ruold,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rcold,.true.)
        call lsyssc_getbase_double(rcold,p_cold)
        call lsyssc_getbase_double(ruold,p_uold)

        ! allocate and create structure of some important matrices and vectors
        call chemo_creatematvec ( rmassmatrix, rsysmatrix, rlaplace, rmatrixchemo, rmatrGradX, rmatrGradY, &
                                  rchemoattract, rrhschemo, rrhschemo_notCoupled, rcell, rrhscell, &
                                  rdef, rdiscretisation, rfc, rfu, ranalytChemo, ranalytCell, rtempAnalyt, &
				  rchemo_x, rchemo_y )
                                            
        ! initialize some matrices to be used later (optimization step)
        call chemo_initmat (rmassmatrix, rsysmatrix, rlaplace, rmatrGradX, rmatrGradY, dtstep, D_1, D_2)
            
        ! set initial conditions for rcell and chemoattract.
        call chemo_initIC (rchemoattract, rcell, ranalytChemo, ranalytCell, iResume)

        ! Setting the boundary conditions for chemoattractant (filtering technique to be used later)
         call chemo_initBC ( rdiscreteBCchemo , rboundary, rdiscretisation, rtriangulation, C_0 , U_0 , INITSOL )

        ! here we set a pointer for the later filtering technique
        rvectorBlockchemo%p_rdiscreteBC => rdiscreteBCchemo

        ! Setting the Blockvectors (since the linear solver works only for block matrices/vectors)
        call lsysbl_createMatFromScalar (rsysmatrix,rmatrixBlockchemo,rdiscretisation)
        call lsysbl_createVecFromScalar (rchemoattract,rvectorBlockchemo,rdiscretisation)
        call lsysbl_createVecFromScalar (rrhschemo, rrhsBlockchemo, rdiscretisation)

        ! Setting some preliminaries for the calculations of the cell density
        call lsysbl_createVecFromScalar (rcell,rcellBlock,rdiscretisation)
        call lsysbl_createVecFromScalar (rdef,rdefBlock,rdiscretisation)
        call lsysbl_createVecFromScalar (rfu,rfuBlock,rdiscretisation)
            
        ! Setting the boundary conditions for cells (filtering technique to be used later)
        call cell_initBC ( rdiscreteBC, rboundary, rdiscretisation, rtriangulation, C_0, U_0, INITSOL )

        ! here we set a pointer for the later filtering technique
        rcellBlock%p_rdiscreteBC => rdiscreteBC
        rdefBlock%p_rdiscreteBC => rdiscreteBC

        ! printing out the initial conditions into a gmv_file
        select case ( gmvfolder )
            case ( 0 )
                call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvcpld/'//trim(gmv_prefix) // 'solution_init.gmv')
            case ( 3 )
                call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvmasscons/'//trim(gmv_prefix) // 'solution_init.gmv')
            case default
                call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvcpld2/'//trim(gmv_prefix) // 'solution_init.gmv')
        end select
        call lsyssc_getbase_double (rcell,p_Ddata)
        call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
        call lsyssc_getbase_double (rchemoattract,p_Ddata)
        call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
        !analytical part
        call lsyssc_getbase_double (ranalytCell,p_Ddata)
        call ucd_addVariableVertexBased (rexport,'analytCells',UCD_VAR_STANDARD, p_Ddata)
        call lsyssc_getbase_double (ranalytChemo,p_Ddata)
        call ucd_addVariableVertexBased (rexport,'analytChemo',UCD_VAR_STANDARD, p_Ddata)
        ! Write the file to disc, that's it.
        call ucd_write (rexport)
        call ucd_release (rexport)

        ! Setting the time init stats
        time_accum = 0.0_DP
        time_max = 0.0_DP
        time_min = 100.0_DP
        ! Resetting the iteration stats
        iteration_c_max = 0
        iteration_c_min = 500
        iteration_c_average = 0
        iteration_u_max = 0
        iteration_u_min = 500
        iteration_u_average = 0
        iteration_defcorr_max = 0
        iteration_defcorr_min = 500
        iteration_defcorr_average = 0

        ! RS: to collect some convergence stats, e.g. timederiv and abs error
        !       they will be written to certain files (specified in the .dat file)
        if (steps /=0) then
            allocate(steady_state_array(2,steps))
            allocate(absolute_error(2,steps))
        else
            allocate(steady_state_array(2,ntimesteps))
            allocate(absolute_error(2,ntimesteps))
        end if

        ! to determine negativ solution values.
        quit =.false.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! We will now start assembling the solutionvectors by first computing a new value of        !!
        !! the chemoattractant vector, e.g. c_n --->c_{n+1}. Then we use this updated value for the  !!
        !! attractant to compute the actual solution of our problem,e.g. the solution vector u_{n+1} !!
        !! ( u_n, c_{n+1} ---> u_{n+1} ). We will use in both cases am implicit Euler method.        !!
        !!                                                                                           !!
        !! STEP 1 : compute c                                                                        !!
        !! STEP 2 : compute u                                                                        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(steps /= 0) then
            ntimesteps = steps
        end if

        auxStepsDigits = 5
        if(ntimesteps+istartIteration > 99999) then
            auxStepsDigitsReal = Real( auxStepsDigits )
            auxStepsDigitsReal = log10( auxStepsDigitsReal )
            auxStepsDigits = 1+Int( auxStepsDigitsReal )
        end if
        if ( iResume /=1 ) then
            istartIteration = 0
        end if

        !create lumped mass matrix
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rtempLump)
        call lsyssc_copyMatrix(rmassmatrix, rtempLump)
        call lsyssc_lumpMatrixScalar (rtempLump,LSYSSC_LUMP_DIAG,.false.)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! Start the timeloop !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        timeloop : do itimestep=1,ntimesteps

            ! time stats
            call cpu_time(time_start)

            ! Next time step.
            dtime = dtime + dtstep

            ! Store the old sols for the errorctrl
            call lsyssc_getbase_double(rcell, p_vectordata)
            call lsyssc_getbase_double(rchemoattract, p_chemodata)
            call lalg_copyVector(p_vectordata, p_uold)
            call lalg_copyVector(p_chemodata, p_cold)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!   STEP 1: we solve for c (chemoattractant)   !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !********** we calculate the rhs of c^n = M c^n + \Delta t f_c ******!
            call chemo_initrhsC_forCoupling ( rchemoattract, rrhschemo_notCoupled, rcell, ranalytCell,&
                                     rdiscretisation, rmassmatrix, dtstep, PHI, rfc)

            !********** we calculate the rhs of u^n = M u^n + \Delta t f_u ******!
            ! Now form the actual RHS by matrix vector multiplication: rrhscell = rmatrix_lumped * u_{old}
            ! The massmatrix has already been built outside the loop
            call lsyssc_scalarMatVec(rtempLump,rcell,rrhscell,1.0_DP,0.0_DP)
!             call lsyssc_scalarMatVec(rmassmatrix,rcell,rrhscell,1.0_DP,0.0_DP)
            ! Here I embed (analytical) term F_u!
!             call collct_init (rcollection)
!             rcollection%DquickAccess(1) = dtstep
!             rcollection%DquickAccess(2) = convecRelaxation
!
!             rlinform%itermCount = 1
!             rlinform%Idescriptors(1) = DER_FUNC
!             call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1), &
!                                             rlinform, .false., rrhscell, callback_chemo_rfu, rcollection)
!
!             ! Release the collection/tempmatrix structure
!             call collct_done(rcollection)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! Start the coupledloop !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !----> seems we are ready for the coupled loop
            coupledloop : do icouplestep=1,1
                
                ! set coupled rhschemo
                call lsyssc_copyVector(rrhschemo_notCoupled, rrhschemo)
!                 call lsyssc_scalarMatVec(rtempLump, rcell, rrhschemo, dtstep*PHI, 1.0_DP)

                ! RS: trying to set it up via a linear-form
                ! e.g. invoke the linf_buildVectorScalar with rcell as the coefficient
                ! call collct_init (rcollection)
                ! call collct_setvalue_vecsca ( rcollection , 'cbvector' , rcell , .true.)
                ! rcollection%DquickAccess(1) = dtstep
                ! rcollection%DquickAccess(2) = PHI
                ! call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1), &
                !                                 rlinform, .false., rrhschemo, coeff_hillenX_RHS_c, rcollection)
                ! call collct_done(rcollection)
!                 ! <--- This actually doesn't any good :(

                call lsyssc_scalarMatVec(rtempLump, rcell, rrhschemo, dtstep*PHI, 1.0_DP)
!                 call lsyssc_scalarMatVec(rmassmatrix, rcell, rrhschemo, dtstep*PHI, 1.0_DP)
                !call lsyssc_scalarMatVec(rtempLump, ranalytCell, rrhschemo, dtstep*PHI, 1.0_DP)

                rmatrixBlockchemo%p_rdiscreteBC => rdiscreteBCchemo
                rrhsBlockchemo%p_rdiscreteBC => rdiscreteBCchemo

                ! Next step is to implement boundary conditions into the RHS, solution and matrix.
                !This is done using a vector/matrix filter for discrete boundary conditions.
                ! The discrete boundary conditions are already attached to the vectors/matrix.
                call matfil_discreteBC (rmatrixBlockchemo)
                call vecfil_discreteBCrhs (rrhsBlockchemo)
                call vecfil_discreteBCsol (rvectorBlockchemo)

                !To calculate DEFECT for chemo
                call lsysbl_copyVector( rrhsBlockchemo, rdefBlock)
                call lsysbl_blockMatVec( rmatrixBlockchemo, rvectorBlockchemo, rdefBlock, -1.0_DP, 1.0_DP )
                defect = lsysbl_vectorNorm ( rdefBlock , LINALG_NORML2 )

                !!!!!! seems everything to be ready for solving c !!!!!!

                ! Setting some resolver issues
                ! During the linear solver, the boundary conditions are also
                ! frequently imposed to the vectors. But as the linear solver
                ! does not work with the actual solution vectors but with
                ! defect vectors instead.
                RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

                ! Create a BiCGStab-solver. Attach the above filter chain
                ! to the solver, so that the solver automatically filters
                ! the vector during the solution process.
                p_RfilterChain => RfilterChain
                nullify(p_rpreconditioner)
                call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
                !other solver: call linsol_initGMRES (p_rsolverNode,64,p_rpreconditioner,p_RfilterChain)
                    
                ! Attach the system matrix to the solver.
                ! First create an array with the matrix data (on all levels, but we
                ! only have one level here), then call the initialisation
                ! routine to attach all these matrices.
                Rmatrices = (/rmatrixBlockchemo/)
                call linsol_setMatrices(p_RsolverNode,Rmatrices)
                    
                ! Initialise structure/data of the solver. This allows the
                ! solver to allocate memory / perform some precalculation to the problem.
                call linsol_initStructure (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) stop
                call linsol_initData (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) stop
                    
                ! STEP 1.6: Solve the system
                ! Now we have block vectors for the RHS and the matrix. What we
                ! need additionally is a block vector for the solution and
                ! temporary data. Create them using the RHS as template.
                call lsysbl_createVecBlockIndirect (rrhsBlockchemo, rtempBlock, .false.)
                ! Finally solve the system. As we want to solve Ax=b with
                ! b being the real RHS and x being the real solution vector,
                ! we use linsol_solveAdaptively. If b is a defect
                ! RHS and x a defect update to be added to a solution vector,
                p_rsolverNode%ioutputLevel = 2
                p_rsolverNode%depsRel=1E-15_DP
                p_rsolverNode%nminIterations=1
                p_rsolverNode%nmaxIterations=1000
                call linsol_solveAdaptively (p_rsolverNode,rvectorBlockchemo,rrhsBlockchemo,rtempBlock)

                ! Store the iterationstats
                if ( p_rsolverNode%iiterations >=  iteration_c_max ) then
                    iteration_c_max = p_rsolverNode%iiterations
                end if
                if ( p_rsolverNode%iiterations <= iteration_c_min ) then
                    iteration_c_min = p_rsolverNode%iiterations
                end if
                iteration_c_average = iteration_c_average + p_rsolverNode%iiterations

                ! We  compute the norm of the difference between two successive timesteps.
                ! (a basic error control below...)
                cerror =  lalg_errorNorm (p_cold,p_chemodata, CTRLNORM)
                            
                ! Release the block matrix/vectors
                call lsysbl_releaseVector (rtempBlock)
                ! Release solver data and structure
                call linsol_doneData (p_rsolverNode)
                call linsol_doneStructure (p_rsolverNode)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!   STEP 2: we solve for u (chemoattractant)   !!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                ! We calculate some iteratives to suit the residual tolerance
                print*,"----------------"
                print*,"timestep : ", itimestep
                print*,"----------------"

                ! defect correction loop. Since the equation for u can be non-linear, we have to perform nonlinear loops
                ! inside the defect correction subroutine.
                call  chemo_defcorr( rcell, rchemoattract, rcellBlock, rmassmatrix, rtempLump,&
                                     rlaplace, rmatrGradX, rmatrGradY, &
                                     rrhscell, rdiscretisation, rtriangulation,rdiscreteBC, dtstep, D_1, &
				     CHI, ALPHA, r, GAMMA, N, maxiterationdef, defectTol ,&
                                     iteration_u_max, iteration_u_min, iteration_u_average, iteration_defcorr_max, &
                                     iteration_defcorr_min, iteration_defcorr_average, rfu, rfuBlock, itimestep, &
				     convecRelaxation, rchemo_x, rchemo_y)
        end do coupledloop
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! end of the coupledloop  !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! STEP 2.7: Postprocessing
            ! That's it, rcellBlock now contains our solution. We can now start the postprocessing.
            if(OUTPUT .eq. 1) then
                ! Start UCD export to GMV file:
                select case (gmvfolder)
                    case (0)
                        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                    'gmvcpld/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                    case (3)
                        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                    'gmvmasscons/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                    case default
                            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                    'gmvcpld2/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                end  select
                    call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)

                    ! Write the file to disc, that's it.
                    call ucd_write (rexport)
                    call ucd_release (rexport)
            end if


            ! Now we take care of printing the gmvs of every istep_gmv-th timestep
            if(OUTPUT .eq. 3) then
                if ( mod( itimestep, istep_gmv ) == 0 ) then
                    ! Start UCD export to GMV file:
                    select case (gmvfolder)
                        case (0)
                        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                        'gmvcpld/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                        case (3)
                        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                        'gmvmasscons/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                        case default
                            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                        'gmvcpld2/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                    end  select
                    call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
                    !for test call lsyssc_getbase_double (rfu,p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    !for test call lsyssc_getbase_double (rfc,p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
                    !analytical part
                    call lsyssc_getbase_double (ranalytCell,p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'analytCells',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (ranalytChemo,p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'analytChemo',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rchemo_x,p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cx',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rchemo_y,p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cy',UCD_VAR_STANDARD, p_Ddata)

                    call ucd_write (rexport)
                    call ucd_release (rexport)
                end if
            end if

            ! We  compute the norm of the difference between two successive timesteps.
            ! (a basic error control below...)
            cerror = lalg_errorNorm (p_cold,p_chemodata, CTRLNORM )
            cerror=cerror/dtstep
            uerror = lalg_errorNorm (p_uold,p_vectordata, CTRLNORM )
            uerror=uerror/dtstep
            steady_state_array(1,itimestep) = uerror
            steady_state_array(2,itimestep) = cerror
            print *,' ############  Rel_chemo_diff=', cerror
            print *,' ############  Rel_cells_diff=', uerror

            !!!!!!! calculate the norm of |analyt - numerica|
!             ! calculate the difference: (c_analyt - c_numeric)
!             call lsyssc_vectorLinearComb (ranalytChemo,rchemoattract,1.0_DP,-1.0_DP,rtempAnalyt)
!             dnorm = lsyssc_vectorNorm (rtempAnalyt, LINALG_NORML2)
!             print*,''
!             print *,' >>>>>>>>>>  Abs_Norms_chemo_diff=', dnorm
!             ! calculate the difference: (u_analyt - u_numeric)
!             call lsyssc_vectorLinearComb (ranalytCell,rcell,1.0_DP,-1.0_DP,rtempAnalyt)
!             dnorm = lsyssc_vectorNorm (rtempAnalyt, LINALG_NORML2)
!             print *,' >>>>>>>>>>  Abs_Norms_cells_diff=', dnorm
!             print*,''
!
!             ! here I calculate the analytical L_2 error
!             call pperr_scalar (rchemoattract,PPERR_L1ERROR,Derr_chemoL2,&
!                         ffunction_Target_Chemo)
!             call pperr_scalar (rcell,PPERR_L1ERROR,Derr_cellsL2,&
!                         ffunction_Target_Cells)
!             print *,' ===========  Diff_L1_chemo=', Derr_chemoL2
!             print *,' ===========  Diff_L1_cells=', Derr_cellsL2
!
!             absolute_error(1,itimestep) = Derr_cellsL2
!             absolute_error(2,itimestep) = Derr_chemoL2
            ! here I calculate the analytical H_1 error
!             call pperr_scalar (rchemoattract,PPERR_H1ERROR,Derr_chemoH1,&
!                         ffunction_Target_ChemoH1)
!             call pperr_scalar (rcell,PPERR_H1ERROR,Derr_cellsH1,&
!                         ffunction_Target_CellsH1)
!             print *,' \\\\\\\\\\\  Diff_H1_chemo=', Derr_chemoH1
!             print *,' \\\\\\\\\\\  Diff_H1_cells=', Derr_cellsH1

            ! time statistics
            call cpu_time(time_stop)
            time_accum = time_accum+( time_stop-time_start )

            if ( time_stop-time_start .ge. time_max) then
                time_max = time_stop-time_start
            end if

            if ( time_stop-time_start .le. time_min) then
                time_min = time_stop-time_start
            end if

            ! If we compute a negative solution, we'll export the actual negativ solutions and exit the computations
            if ( checkneg .eq. 1) then
            checkneg_loop :  do i=1,rtriangulation%NVT
                if( p_vectordata(i) < negthres .AND. .NOT.quit) then
                        select case (gmvfolder)

                        case (0)
                            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                'gmvcpld/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                        case (3)
                            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                'gmvmasscons/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                        case default
                            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                                'gmvcpld2/'//trim(gmv_prefix) // 'solution.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                    end  select
                    call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
                        ! Write the file to disc, that's it.
                    call ucd_write (rexport)
                    call ucd_release (rexport)
                    quit = .true.
                    end if
            end do checkneg_loop
            end if
                                            
            ! The simulation should stop if we reached a nearly steady state. So thats whats done here...
            ! We approximate our time-derivative with the first order accurate difference quotient.
            ! d sol/ dt ~ (sol_n - sol_{n-1})/dtstep
            if(steps == 0) then
                if(tol >= uerror .AND. tol >=cerror)  then
                    print *, "tolerance threshold reached-----simulation ends up in 'nearly' steady state."
                    print *,"############ differnece to steady state in c ################"
                    print *,cerror
                    print *,"############################"
                    print *,"############## differnece to steady state in u ##############"
                    print *, uerror
                    print *,"############################"
                    quit = .true.
                end if
            end if
                                    
            if( quit ) then
                EXIT
            end if
            
        end do timeloop
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! end of the global time loop !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! RS: Here we write the steady state stats to the files
!
! ----------------------------
!
       if(steps .NE. 0 .OR. (.NOT.(quit))) then
            itimestep = itimestep -1
        end if
 ! write the steady_state_array in a file
print*,"filename = ",stat_ouput
    call io_openFileForWriting(stat_ouput, iunit, SYS_REPLACE)
        if (iunit .eq. -1) then
            call error_print(ERR_IO_FILEIO,"io_openFileForWriting", ERR_CRITICAL)
        else
            write (iunit, '(A)' ) "steady_state_array"
            write (iunit, '(A, ES8.1, A, ES8.1)') "ALPHA= ", convecRelaxation, "     BETA= ", PHI
            write (iunit, '(A, I, A, ES8.1)') "Level = ", NLMAX, "     dtstep = ", dtstep
            write (iunit, '(A, ES8.1, A, I)') "Tol = ", TOL, "     norm = ", CTRLNORM
            write (iunit,'(A, I, A, I)' ) "gmvoutput = ", gmvfolder, "ISTEP_GMV = ", ISTEP_GMV
            write (iunit, '(A, I, A, I)' ) "iResume = ", iresume,"     istartIteration = ",istartIteration
            write (iunit, '(A)') "steady state of u :"
            do iwrite=1,itimestep
                write (iunit, '(E20.10)') steady_state_array(1,iwrite )
            end do
!             if (steps ==0 .AND. quit) then
! !                 print *, steady_state_array(1,itimestep )
!                 write (iunit, '(E20.10)') steady_state_array(1,itimestep )
!             end if
            write (iunit, '(A)') "steady state of c :"
            do iwrite=1,itimestep
                write (iunit, '(E20.10)') steady_state_array(2,iwrite )
            end do
!             if (steps ==0 .AND. quit) then
!                 write (iunit, '(E20.10)') steady_state_array(2,itimestep )
!             end if

        end if

    call io_openFileForWriting(stat_ouput//"_abs_error", iunit, SYS_REPLACE)
        if (iunit .eq. -1) then
            call error_print(ERR_IO_FILEIO,"io_openFileForWriting", ERR_CRITICAL)
        else
            write (iunit, '(A)' ) "absolute_error"
            write (iunit, '(A, ES8.1, A, ES8.1)') "ALPHA= ", convecRelaxation, "     BETA= ", PHI
            write (iunit, '(A, I, A, ES8.1)') "Level = ", NLMAX, "     dtstep = ", dtstep
            write (iunit, '(A, ES8.1, A, I)') "Tol = ", TOL, "     norm = ", CTRLNORM
            write (iunit,'(A, I, A, I)' ) "gmvoutput = ", gmvfolder, "ISTEP_GMV = ", ISTEP_GMV
            write (iunit, '(A, I, A, I)' ) "iResume = ", iresume,"     istartIteration = ",istartIteration
            write (iunit, '(A)') "|| u_n - u^* || : "
            do iwrite=1,itimestep
                write (iunit, '(E20.10)') absolute_error(1,iwrite )
            end do
!             if (steps ==0 .AND. quit) then
!                 write (iunit, '(E20.10)') absolute_error(1,itimestep )
!             end if
            write (iunit, '(A)') "|| c_n - c^* || : "
            do iwrite=1,itimestep
                write (iunit, '(E20.10)') absolute_error(2,iwrite )
            end do
!             if (steps ==0 .AND. quit) then
!                 write (iunit, '(E20.10)') absolute_error(2,itimestep )
!             end if

        end if

!
! ----------------------------
!

        ! deallocate lumped mass matrix
         call lsyssc_releaseMatrix(rtempLump)

        ! If we want to export just the last vectors
        if ( OUTPUT .eq. 2 ) then
            select case (gmvfolder)
                case (0)
                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                            'gmvcpld/'//trim(gmv_prefix) // 'solutionEND.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                case (3)
                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                            'gmvmasscons/'//trim(gmv_prefix) // 'solutionEND.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
                case default
                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                            'gmvcpld2/'//trim(gmv_prefix) // 'solutionEND.gmv.'//trim(sys_si0L(itimestep+istartIteration,auxStepsDigits)))
            end select
            call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
            call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!   Printing some statistical information   !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *," patternformation_FCT:"
        print *,"----------------------------"
        print *,"Simulation parameters:"
        if(steps /= 0)then
            print*,"fixed amount of time steps = ",steps
        end if
        print*,"iResume = ", iresume
        print*, "istartIteration = ",istartIteration
        print*, "TIMESTEP = ",dtstep
        print *, "NLMAX = ",NLMAX
        print*,"def_cor Tol : ",defectTol
        print*,"def_cor iterations : ",maxiterationdef
        print*, "gmvoutput folder : ", gmvfolder
        print*,"ISTEP_GMV = ", ISTEP_GMV

        print*,"model parameters:"
        print *, "D_1 = ",D_1
        print *, "D_2 = ",D_2
        print *, "CHI = ",CHI
        print *, "ALPHA = ",ALPHA
        print *, "BETA = ",BETA
        print *, "N = ",N
        print *, "R = ",R

        print *, "---------------------------"
!         if ( quit ) then
            print *, "stopped @  timestep = ", itimestep," out of ",ntimesteps
!         else
!             print *, "stopped @  timestep = ", itimestep-1," out of ",ntimesteps
!         end if
        print*, "----------------------------"
        print*, "--------TIMESTATS-----------"
        print*, "max time per iteration : ",time_max
        print*, "min time per iteration: ", time_min
        print*, "----------------------------"
        print*, "------ITERATION STATS-------"
        print*, "SOLVING FOR C :"
        print*,"max iteration needed :",iteration_c_max
        print *,"min iteration needed :", iteration_c_min
!         if ( quit ) then
            print*, "average iterations needed :", iteration_c_average / itimestep
!         else
!             print*, "average iterations needed :", iteration_c_average / (itimestep-1)
!         end if
        print*, "SOLVING FOR U :"
        print*,"max iteration needed :", iteration_u_max
        print *,"min iteration needed :", iteration_u_min
!         if ( quit ) then
            print*, "average iterations needed :", iteration_u_average / itimestep
!         else
!             print*, "average iterations needed :", iteration_u_average / (itimestep-1)
!         end if
        print*, "SOLVING FOR U ( DEF_CORR ) :"
        print*,"max iteration needed :", iteration_defcorr_max
        print *,"min iteration needed :", iteration_defcorr_min
!         if ( quit ) then
            print*, "average iterations needed :", iteration_defcorr_average / itimestep
!         else
!             print*, "average iterations needed :", iteration_defcorr_average / (itimestep-1)
!         end if
                            
        print *,"############ differenece to steady state in c ################"
                        print *,cerror
                        print *,"############################"
                        print *,"############## differnece to steady state in u ##############"
                        print *, uerror
                        print *,"############################"
        print*, "----------------------------"

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!   Releasements of all vector.matrix and other data structures   !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        deallocate (steady_state_array)
        deallocate (absolute_error)

        call lsysbl_releaseVector (rrhsBlockchemo)

        ! releasing the chemoBlockstructures
        call lsysbl_releaseVector (rvectorBlockchemo)
        call lsysbl_releaseMatrix (rmatrixBlockchemo)
        call lsysbl_releaseVector (rdefBlock)
        ! could result in an error since template of RHS != RHS_CHEMO ?

        ! Releasing a cellblockstructure
        call lsysbl_releaseVector (rcellBlock)

        ! Release solver data and structure
        call linsol_doneData (p_rsolverNode)
        call linsol_doneStructure (p_rsolverNode)
            
        ! Release the solver node and all subnodes attached to it (if at all):
        call linsol_releaseSolver (p_rsolverNode)

        ! Release our discrete version of the boundary conditions
        call bcasm_releaseDiscreteBC (rdiscreteBCchemo)
        call bcasm_releaseDiscreteBC (rdiscreteBC)

        ! Release the preallocated matrix and the solution vector.
        call lsyssc_releaseMatrix (rmatrix)
        call lsyssc_releaseMatrix (rmatrixchemo)
        call lsyssc_releaseMatrix (rmassmatrix)
        call lsyssc_releaseMatrix (rsysmatrix)
        call lsyssc_releaseMatrix (rlaplace)
        call lsyssc_releaseMatrix (rmatrGradX)
        call lsyssc_releaseMatrix (rmatrGradY)
        call lsyssc_releaseVector (rcell)
        call lsyssc_releaseVector (rchemo_x)
        call lsyssc_releaseVector (rchemo_y)
        call lsyssc_releaseVector (rrhschemo)
        call lsyssc_releaseVector (rrhschemo_notCoupled)
        call lsyssc_releaseVector (rchemoattract)
        call lsyssc_releaseVector (rdef)
        call lsyssc_releaseVector (rrhscell)
        call lsyssc_releaseVector (rfc)
        call lsyssc_releaseVector (rfu)
        call lsysbl_releaseVector (rfuBlock)

        call lsyssc_releaseVector (rcold)
        call lsyssc_releaseVector (ruold)

        ! comparision with analytical solution
        call lsyssc_releaseVector (ranalytChemo)
        call lsyssc_releaseVector (ranalytCell)
        call lsyssc_releaseVector (rtempAnalyt)
        
        ! Release parameterlist
        call parlst_done (rparams)

        ! Release the discretisation structure and all spatial discretisation structures in it.
        call spdiscr_releaseBlockDiscr(rdiscretisation)

        ! Release the triangulation.
        call tria_done (rtriangulation)

        ! Finally release the domain, that's it.
        call boundary_release (rboundary)

        print *,">>end of the subroutine chemotaxispatternFCT ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemotaxispatternFCT subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################



    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_initmat subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Here we assemble the matrices in advance ( to save computation time )
    subroutine chemo_initmat (rmassmatrix, rsysmatrix, rlaplace,  rmatrGradX, rmatrGradY, dtstep , D_1 , D_2)

        ! This is where the matrices shoiuld be stored
        type(t_matrixScalar), intent(INOUT):: rmassmatrix , rlaplace , rsysmatrix
        type(t_matrixScalar) :: rtempmatrix

        type(t_matrixScalar), intent(INOUT) :: rmatrGradX, rmatrGradY
        
        real(DP) , intent (INOUT) :: dtstep , D_1 , D_2
        ! local bilinearform to construct the matrices
        type(t_bilinearForm) :: rform
        
        print *,">>start of the subroutine chemo_initmat ... "
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = 1.0
        call bilf_buildMatrixScalar (rform,.true.,rmassmatrix)

        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        !rform%Dcoefficients(1)  = dtstep * D_1
        rform%Dcoefficients(1)  = 1.0
        !rform%Dcoefficients(2)  = dtstep * D_1
        rform%Dcoefficients(2)  = 1.0
        call bilf_buildMatrixScalar (rform,.true.,rlaplace)

        call lsyssc_copymatrix(rmassmatrix,rsysmatrix)
        !lumping the sysmatrix TESTAS
        call lsyssc_lumpMatrixScalar (rsysmatrix,LSYSSC_LUMP_DIAG,.false.)
        call lsyssc_copymatrix(rsysmatrix,rtempmatrix)

        ! Since we want to save time, we use already constructed matrices
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        !1.0_DP is exactly my alpha
        !!!TEST just for test
        !rform%Dcoefficients(1)  = dtstep * 0.0_DP
        rform%Dcoefficients(1)  = dtstep * 1.0_DP
        !rform%Dcoefficients(1)  = dtstep * 32.0_DP  ! This is the coefficient in the paper
        !test rform%Dcoefficients(1)  = 0.0_DP
        !call bilf_buildMatrixScalar (rform,.false.,rsysmatrix)
        call lsyssc_matrixLinearComb(rtempmatrix, dtstep * 1.0_DP, rsysmatrix, 1.0_DP, rsysmatrix,.false.,.false.,.true.,.true.)

        ! Since we want to save time, we use already constructed matrices
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = dtstep * D_2
        rform%Dcoefficients(2)  = dtstep * D_2
        call bilf_buildMatrixScalar (rform,.false.,rsysmatrix)
        
        call lsyssc_releaseMatrix (rtempmatrix)

        ! we evaluate gradient c_x
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff(1) = .true.
        rform%Dcoefficients(1)  = 1.0
        call bilf_buildMatrixScalar (rform,.true.,rmatrGradX)

        ! we evaluate gradient c_y
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_DERIV_Y
        rform%Idescriptors(2,1) = DER_FUNC
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff(1) = .true.
        rform%Dcoefficients(1)  = 1.0
        call bilf_buildMatrixScalar (rform,.true.,rmatrGradY)
        
        print *,">>end of the subroutine chemo_initmat ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_initmat subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################


    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_creatematvec subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! In thie routine the basic matrices and vectors are initialized
    subroutine chemo_creatematvec (rmassmatrix, rsysmatrix, rlaplace, rmatrixchemo, rmatrGradX, rmatrGradY, &
                                   rchemoattract, rrhschemo, rrhschemo_notCoupled, rcell, rrhs, rdef, &
                                   rdiscretisation, rfc, rfu, ranalytChemo, ranalytCell, rtempAnalyt, &
				   rchemo_x, rchemo_y )

        ! This is where the matrices shoiuld be stored
        type(t_matrixScalar) , intent(INOUT) :: rmassmatrix, rsysmatrix, rlaplace, rmatrixchemo

        ! This are the vectors
        type(t_vectorScalar), intent(INOUT) :: rchemoattract, rcell, rdef, rrhschemo, rrhs, rfc, rfu,rchemo_x,rchemo_y
        type(t_vectorScalar), intent(INOUT) :: rrhschemo_notCoupled
        type(t_matrixScalar), intent(INOUT) :: rmatrGradX, rmatrGradY
        type(t_vectorScalar), intent(INOUT) :: ranalytChemo, ranalytCell, rtempAnalyt
        
        ! The underlying discretisation
        type(t_blockDiscretisation) , intent (IN) :: rdiscretisation

        print *,">>start of the subroutine chemo_creatematvec ... "
        
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9,rlaplace)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9,rmatrixchemo)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9,rsysmatrix)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9,rmassmatrix)
        ! gradient matrices
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9,rmatrGradX)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9,rmatrGradY)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rchemoattract,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rcell,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rdef,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rrhschemo, .true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rrhschemo_notCoupled, .true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rfc, .true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rfu, .true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rrhs,.true.)
        ! comparision with analytical solution
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),ranalytChemo,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),ranalytCell,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rtempAnalyt,.true.)

        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rchemo_x,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rchemo_y,.true.)

        print *,">>end of the subroutine chemo_creatematvec ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_creatematvec subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################

    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_initIC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This is called to implement the initial conditions, e.g. to set the initial sol vectors
    subroutine chemo_initIC ( rchemoattract, rcell, ranalytChemo, ranalytCell, iResume)

    !<input>
        integer , intent (IN) :: iResume

        ! sol vectors
        type(t_vectorScalar) , intent (INOUT) :: rchemoattract, rcell
        ! for analitycal solution
        type(t_vectorScalar), intent (INOUT) :: ranalytChemo, ranalytCell
    !</input>

    type(t_collection) :: rcollection

        print *,">>start of the subroutine chemo_initIC ... "

                    call collct_init (rcollection)
                    ! Setting the IC for u
                    rcollection%DquickAccess(1) = 1000.0_DP
                    rcollection%DquickAccess(2) = 100.0_DP
                    call anprj_discrDirect ( rcell , coeff_anprj_ic_exp , rcollection )
                    ! Setting the IC for c
                    rcollection%DquickAccess(1) = 500.0_DP
                    rcollection%DquickAccess(2) = 50.0_DP
                    call anprj_discrDirect ( rchemoattract , coeff_anprj_ic_exp , rcollection )

                    call collct_done (rcollection)

        ! now we set the analytical solution
!         call anprj_discrDirect ( ranalytCell, analyt_u_pattern )
!         call anprj_discrDirect ( ranalytChemo, analyt_c_pattern )


        print *,">>end of the subroutine chemo_initIC ... "

    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_initIC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################

    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of cell_initBC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize the boundary conditions (used for both, the cell density sol and the chemoattractant sol)
    subroutine cell_initBC (rdiscreteBC, rboundary, rdiscretisation, rtriangulation, C_0, U_0, INITSOL)

        ! The BC_type which should be initialized
        type(t_discreteBC) , intent (INOUT) ::rdiscreteBC

        type(t_boundary) , intent (IN) :: rboundary

        ! Underlying discretisation
        type(t_blockDiscretisation) , intent (IN) ::rdiscretisation

        ! Underlying triangulation
        type(t_triangulation) , intent (IN) :: rtriangulation

        ! BCs if we take constant initsol
        real(DP) , intent (IN) :: C_0 , U_0

        integer, intent(IN) :: INITSOL

        ! local collection
        type(t_collection) :: rcollection

        type(t_boundaryRegion) :: rboundaryRegion

        ! An object for saving the boundary mesh region
        type(t_meshregion) :: rmeshRegion

        print *,">>start of the subroutine cell_initBC ... "

        ! prepare the collection for the future work
        call collct_init(rcollection)
        rcollection%DquickAccess(1) = 0.0_DP
        rcollection%DquickAccess(2) = 0.0_DP

        ! here the Dirichlet chemo-code
        call bcasm_initDiscreteBC(rdiscreteBC)
          
        call boundary_createRegion(rboundary,1,1,rboundaryRegion)
        ! We use this boundary region and specify that we want to have Dirichlet
        ! boundary there.
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_u_callback, rcollection)
                                
        ! Now to the edge 2 of boundary component 1 the domain.
        call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_u_callback, rcollection)
                                       							
        ! Edge 3 of boundary component 1.
        call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_u_callback, rcollection)

        ! Edge 4 of boundary component 1. That's it.
        call boundary_createRegion(rboundary,1,4,rboundaryRegion)
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_u_callback, rcollection)
        
        ! deactivate collection
        call collct_done(rcollection)

        print *,">>end of the subroutine cell_initBC ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of cell_initBC subroutine  !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################


    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_initBC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize the boundary conditions (used for both, the cell density sol and the chemoattractant sol)
    subroutine chemo_initBC (rdiscreteBC, rboundary, rdiscretisation, rtriangulation, C_0, U_0, INITSOL)

        ! The BC_type which should be initialized
        type(t_discreteBC) , intent (INOUT) ::rdiscreteBC

        type(t_boundary) , intent (IN) :: rboundary

        ! Underlying discretisation
        type(t_blockDiscretisation) , intent (IN) ::rdiscretisation

        ! Underlying triangulation
        type(t_triangulation) , intent (IN) :: rtriangulation

        ! BCs if we take constant initsol
        real(DP) , intent (IN) :: C_0 , U_0

        integer, intent(IN) :: INITSOL

        ! local collection
        type(t_collection) :: rcollection

        type(t_boundaryRegion) :: rboundaryRegion

        ! An object for saving the boundary mesh region
        type(t_meshregion) :: rmeshRegion

        print *,">>start of the subroutine chemo_initBC ... "

        ! prepare the collection for the future work
        call collct_init(rcollection)
        rcollection%DquickAccess(1) = 0.0_DP
        rcollection%DquickAccess(2) = 0.0_DP

        ! here the Dirichlet chemo-code
        call bcasm_initDiscreteBC(rdiscreteBC)
          
        call boundary_createRegion(rboundary,1,1,rboundaryRegion)
        ! We use this boundary region and specify that we want to have Dirichlet
        ! boundary there.
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_c_callback, rcollection)
                                
        ! Now to the edge 2 of boundary component 1 the domain.
        call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_c_callback, rcollection)
                                       							
        ! Edge 3 of boundary component 1.
        call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_c_callback, rcollection)

        ! Edge 4 of boundary component 1. That's it.
        call boundary_createRegion(rboundary,1,4,rboundaryRegion)
!         call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                         rboundaryRegion,rdiscreteBC,&
!                                         getBoundaryValues_c_callback, rcollection)
        
        ! deactivate collection
        call collct_done(rcollection)

        print *,">>end of the subroutine chemo_initBC ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_initBC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################


    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_initrhsC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Constructing the RHS of the chemoattractant part
    subroutine chemo_initrhsC_forCoupling  ( rchemoattract, rhschemo_notCoupled, rcell, &
                                ranalytCell, rdiscretisation, rmassmatrix,  dtstep,  PHI, rfc)

        type(t_vectorScalar) , intent(INOUT) :: rchemoattract , rcell, ranalytCell

        type(t_matrixScalar) , intent(IN) :: rmassmatrix
        type(t_matrixScalar) :: rlumpedmass

        type(t_vectorScalar) , intent(INOUT) :: rhschemo_notCoupled
        type(t_vectorScalar) , intent(INOUT) :: rfc

        type(t_Blockdiscretisation) , intent(IN) :: rdiscretisation

        ! param
        real(DP) ,intent (IN) :: dtstep, PHI

        ! local linform
        type(t_linearForm) :: rlinform

        ! local collection
        type(t_collection) :: rcollection

        ! for test
        real(DP), dimension(:), pointer :: p_vector

        print *,">>start of the subroutine chemo_initrhsC ... "
        
        ! To assemble the RHS , set up the corresponding linear  form (u*g(u),Phi_j):
        !rlinform%itermCount = 1
        !rlinform%Idescriptors(1) = DER_FUNC
        !call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1), &
        !                            rlinform, .true., rhschemo, coeff_hillenX_RHS_c, rcollection)
        !
        call lsyssc_copyMatrix ( rmassmatrix, rlumpedmass )
        call lsyssc_lumpMatrixScalar (rlumpedmass,LSYSSC_LUMP_DIAG,.false.)

        ! rmassmatrix* c_{old} (the corressponding massmatrix is already built)
        call lsyssc_scalarMatVec(rlumpedmass, rchemoattract, rhschemo_notCoupled, 1.0_DP, 0.0_DP)

        ! Here I embed the (analytical) term F_c!
!         call collct_init (rcollection)
! !         call collct_setvalue_vecsca ( rcollection , 'cbvector1' , rcell , .true.)
!         call collct_setvalue_vecsca ( rcollection , 'cbvector2' , ranalytCell , .true.)
!         rcollection%DquickAccess(1) = PHI
!         rlinform%itermCount = 1
!         rlinform%Idescriptors(1) = DER_FUNC
!         call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1), &
!                                     rlinform, .true., rfc, callback_chemo_rfc, rcollection)
!         ! f_c -> rhschemo
!         call lsyssc_vectorLinearComb (rfc, rhschemo_notCoupled, dtstep, 1.0_DP)
!
!         ! Release the collection structure
!         call collct_done(rcollection)
        ! Now we should have the correct RHS = rhs like its mentioned above.
        call lsyssc_releaseMatrix (rlumpedmass)

        print *,">>end of the subroutine chemo_initrhsC ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_initrhsC subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################


    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_defcorr subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine chemo_defcorr(  rcell , rchemoattract , rcellBlock , rmassmatrix , rtempLump, &
                 rlaplace, rmatrGradX, rmatrGradY, &
                rrhscell , rdiscretisation, rtriangulation, rdiscreteBC, dtstep, D_1, CHI, ALPHA, r, &
	        GAMMA, N, maxiterationdef , defectTol,&
                iteration_u_max, iteration_u_min, iteration_u_average, iteration_defcorr_max, &
                iteration_defcorr_min, iteration_defcorr_average, rfu, rfuBlock, itimestep, convecRelaxation, &
		rchemo_x, rchemo_y)


        ! sol vectors
        type(t_vectorScalar) ,  intent(INOUT) :: rcell , rchemoattract
        type(t_vectorScalar) ,  intent(INOUT) :: rchemo_x , rchemo_y
        ! The solution block - vector
        type(t_vectorBlock) , intent(INOUT) :: rcellBlock
        ! Some matrices for computing the LHS
        type(t_matrixScalar), intent(INOUT) :: rmassmatrix, rlaplace, rtempLump
        type(t_matrixScalar), intent(INOUT) :: rmatrGradX, rmatrGradY
        ! The RHs of the init pb
        type(t_vectorScalar)  , intent(IN) :: rrhscell
        integer, intent(IN) :: itimestep

        type(t_vectorScalar)  , intent(INOUT) :: rfu
        type(t_vectorBlock)  , intent(INOUT) :: rfuBlock

        ! The underlying discretisation and BCs
        type(t_Blockdiscretisation) , intent(IN) :: rdiscretisation
        type(t_triangulation) :: rtriangulation
        type(t_discreteBC) , intent(IN) , target :: rdiscreteBC
        ! Some params (needed to derive the sys matrix)
        real(DP) , intent(IN) :: dtstep , D_1, CHI, ALPHA, r, GAMMA, N
        integer , intent(IN) :: maxiterationdef
        real(DP) , intent(IN) ::defectTol
        integer , intent(INOUT) :: iteration_u_max, iteration_u_min,&
                                            iteration_defcorr_max, iteration_defcorr_min
        real(DP), intent(INOUT) ::  iteration_u_average, iteration_defcorr_average

        ! Some local variables
        type(t_vectorBlock) :: rrhsBlock , rdefBlock
        type(t_matrixBlock) :: rmatrixBlock
        type(t_collection) :: rcollection
        type(t_bilinearForm) :: rform

        ! local linform
        type(t_linearForm) :: rlinform

        ! temporal matrix
        type(t_matrixScalar) :: rK , rmatrix, rlumpedmass, rtempLinComb
        ! defect ctrl variable
        real(DP) :: defect
        ! A solver node that accepts parameters for the linear solver
        type(t_linsolNode), pointer :: p_rsolverNode_cells,p_rpreconditioner_cells
        ! An array for the system matrix(matrices) during the initialisation of
        ! the linear solver.
        type(t_matrixBlock), dimension(1) :: Rmatrices
        ! A filter chain that describes how to filter the matrix/vector
        ! before/during the solution process. The filters usually implement
        ! boundary conditions.
        type(t_filterChain), dimension(1), target :: RfilterChain
        type(t_filterChain), dimension(:), pointer :: p_RfilterChain
        ! Error indicator during initialisation of the solver
        integer :: ierror
        integer :: k, nedge
        real(DP), dimension(:), pointer :: p_cell , p_tempVec
        real(DP) :: korrection
        real(DP), dimension(:), allocatable ::dedge
        integer, dimension(:,:), allocatable ::kedge

        ! a local average iterationcounter for  iteration_u_average
        real(DP) :: u_average_local

        ! fot test
        real(DP), dimension(:), pointer :: p_vector
                    
        real(DP) :: dnorm_x, dnorm_y
        integer :: i
        
        ! relaxation for the convective/chemotactical term
        real(DP), intent(IN) ::  convecRelaxation
        
        print *,">>start of the subroutine chemo_defcorr ... "
        
        u_average_local = 0.0_DP

        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rtempLinComb)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rK)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9,rmatrix)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!   begin nonlinear loop   !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        def_correction : do k=1,maxiterationdef

            ! STEP 2.2: Assemble the system matrix M_total
            ! Now we are constructing the Systemmatrix for the implicit scheme. @# ( [M + dt*L -dt*M_tri] u_{n+1} ) #@

            !call testPrintMatrix9(rmassmatrix)
            ! copy the already constructed matrix to rmatrix
            call lsyssc_copyMatrix(rmassmatrix, rmatrix)
            ! The only matrix-part which have to be computed within the timeloop
            ! We use a bilf to obtain a suitable matrix
            ! Since the CHI is free to be nonlinear it should be implemented directly
            ! in the cb fct. we ' ll invoke in assembling the matrix

!             call lsyssc_lumpMatrixScalar (rmatrix,LSYSSC_LUMP_DIAG,.false.)
!             call lsyssc_copyMatrix ( rmatrix, rlumpedmass )

            !!!!!!!!!!!!!!!!!!!!!! Here I evaluate the gradient of c !!!!!!!!!!!!!!!!!!!!!!!!!!
            ! calculate $\nabla_x c$
            call lsyssc_scalarMatVec(rmatrGradX,rchemoattract,rchemo_x,1.0_DP,0.0_DP)
            call multInvLumpedVector(rtempLump, rchemo_x)
            ! Calculate (for test purposes) the l2_norm of c_x
            dnorm_x = lsyssc_vectorNorm (rchemo_x, LINALG_NORML2)
            print *,'!!!!!!!!!!!!!!!!!! c_x=',dnorm_x
                      
            call lsyssc_scalarMatVec(rmatrGradY,rchemoattract,rchemo_y,1.0_DP,0.0_DP)
            !call testPrintVector(rtempVec)
            call multInvLumpedVector(rtempLump, rchemo_y)
            ! Calculate (for test purposes) the l2_norm of c_y
            dnorm_y = lsyssc_vectorNorm (rchemo_y, LINALG_NORML2)
            print *,'!!!!!!!!!!!!!!!!!! c_y=',dnorm_y
            !!!!!!!!!!!!!!!!!!!!!! finish to evaluate the gradient of c (in rtempVect) !!!!!!!!!!!!!!!!!!!!!!!!!!

	    !print *,' in defcorr subroutine itimestep=',itimestep
	    !pause
	    !if(MOD(itimestep,50).eq.0) then
	    !	print *,'writing the c_x into file (defcorr subroutine) ...'
	    !	call vector_writeToFile (rtriangulation, rchemo_x,' c_x',itimestep=itimestep)
	    !	print *,'writing the c_y into file (defcorr subroutine) ...'
	    !	call vector_writeToFile (rtriangulation, rchemo_y,' c_y',itimestep=itimestep)
	    !	!print *,' stop in the defcorr subroutine: writing gradients ...'
	    !endif

            ! Initialize the collection structure
            call collct_init (rcollection)

            ! Setting up the collection structure
            ! We're considering the case of constant chemosensitivityfct CHI
            call collct_setvalue_vecsca ( rcollection , 'cbvector1' , rchemoattract , .true.)
            call collct_setvalue_vecsca ( rcollection , 'cbvector2' , rcell , .true.)
            call collct_setvalue_vecsca ( rcollection , 'rvector_x' , rchemo_x, .true.)
            call collct_setvalue_vecsca ( rcollection , 'rvector_y' , rchemo_y, .true.)
            rcollection%DquickAccess(1) = dtstep
            rcollection%DquickAccess(2) = CHI
            rcollection%DquickAccess(3) = GAMMA
            rcollection%DquickAccess(4) = ALPHA
            rcollection%DquickAccess(5) = convecRelaxation
            rform%itermCount = 2
            rform%Idescriptors(1,1) = DER_FUNC
            rform%Idescriptors(2,1) = DER_DERIV_X
            rform%Idescriptors(1,2) = DER_FUNC
            rform%Idescriptors(2,2) = DER_DERIV_Y
            rform%ballCoeffConstant = .false.
            rform%BconstantCoeff(1) = .false.
            rform%BconstantCoeff(2) = .false.
            call bilf_buildMatrixScalar (rform,.true.,rK, coeff_cherkurbilf, rcollection)
                        
            if(k.eq.1) then
                allocate ( kedge ( 2, rK%NA ) )
                allocate ( dedge ( rK%NA ) )
            endif

            !**********************
            ! Adding some artificial diffusion to obtain positivity + smooth sol.
            ! K->K* (rK)
            call chemo_artdiff( rK, dedge, kedge, nedge )
            !!!call chemo_artdiff( rlumpedmass, rK, dedge, kedge, nedge, aedge_mass )
            !**********************

            ! Adding the logistic growth of u^2(1-u) as proposed in the pattern-paper
            ! Here we implement this term as a LHS NOT a RHS
            ! In this case (not commented) the construction of the RHS_u should be modified
            !rcollection%DquickAccess(1) = dtstep
            !rform%itermCount = 1
            !rform%Idescriptors(1,1) = DER_FUNC
            !rform%Idescriptors(2,1) = DER_FUNC
            !rform%ballCoeffConstant = .false.
            !rform%BconstantCoeff(1) = .false.
            !pattern_formation_not_test call bilf_buildMatrixScalar (rform,.false.,rK, coeff_pattern_growthterm, rcollection)

            !call chemo_artdiff( rmassmatrix, rK, dedge, kedge, nedge, aedge_mass )

            ! Since we use a 1st order approx of the time-deriv. (backward euler) we add the
            ! lumped massmatrix to the existing matrix
            call lsyssc_matrixLinearComb(rK, -dtstep, rtempLump, 1.0_DP, rtempLinComb,.true.,.true.,.true.)
!             call lsyssc_matrixLinearComb(rK, -dtstep, rmassmatrix, 1.0_DP, rtempLinComb,.true.,.true.,.true.)
            call lsyssc_copyMatrix(rtempLinComb, rmatrix)
            call lsyssc_matrixLinearComb(rlaplace , dtstep*D_1, rmatrix , 1.0_DP, rtempLinComb, .true. , .true. , .true.)
            call lsyssc_copyMatrix(rtempLinComb, rmatrix)
            
            ! Adding the diffusion-part of the model, since this can also be non-linear
            ! we use the callback function callback_defcorr_laplace
!             rcollection%DquickAccess(1) = dtstep
!             rcollection%DquickAccess(2) = D_1
!             rcollection%DquickAccess(3) = N
!             rform%itermCount = 2
!             rform%Idescriptors(1,1) = DER_DERIV_X
!             rform%Idescriptors(2,1) = DER_DERIV_X
!             rform%Idescriptors(1,2) = DER_DERIV_Y
!             rform%Idescriptors(2,2) = DER_DERIV_Y
!             rform%ballCoeffConstant = .false.
!             rform%BconstantCoeff(1) = .false.
!             rform%BconstantCoeff(2) = .false.
!             call bilf_buildMatrixScalar (rform,.false.,rmatrix, callback_defcorr_laplace, rcollection)

            ! Now rmatrix is our systemmatrix
            ! Releasing the collection
            call collct_done(rcollection)
                
            ! STEP 2.3: Create block vectors and boundary conditions.
            !
            ! The linear solver only works for block matrices/vectors - but above,
            ! we created scalar ones. So the next step is to make a 1x1 block
            ! system from the matrices/vectors above which the linear solver
            ! understands.
            call lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
            call lsysbl_createVecFromScalar (rrhscell,rrhsBlock,rdiscretisation)

            ! Hang the pointer into the vector and matrix. That way, these
            ! boundary conditions are always connected to that matrix and that
            ! vector.
            rmatrixBlock%p_rdiscreteBC => rdiscreteBC
            rrhsBlock%p_rdiscreteBC => rdiscreteBC

            ! Next step is to implement boundary conditions into the RHS,
            ! solution and matrix.
            call vecfil_discreteBCrhs (rrhsBlock)
            call vecfil_discreteBCsol (rcellBlock)
            call matfil_discreteBC (rmatrixBlock)

            ! Calculate DEFECT  for cells: rdefRHS = rrhs - rmatrix * rvector
            call lsysbl_copyVector( rrhsBlock, rdefBlock)
            call lsysbl_blockMatVec( rmatrixBlock, rcellBlock, rdefBlock, -1.0_DP, 1.0_DP )
            call vecfil_discreteBCdef (rdefBlock)

            !!extra Applying the FCT code
            call chemo_fct_lim ( rcell, kedge, dedge, nedge, &
                                 rK, dtstep, rdefBlock%RvectorBlock(1) )

	    !!extra
            call vecfil_discreteBCdef (rdefBlock)

            ! Calculate the l2_norm of the defect
            defect = lsysbl_vectorNorm ( rdefBlock , LINALG_NORML2 )

            ! Here is the time to CHECK SOME RESIDUAL CONVERGENCE
            ! Certain defect error control should be implemented here
            if ( k .gt. 1 ) then
                if ( defect .le. defectTol) then
                        ! Cells Releasements
                        call lsysbl_releaseVector(rdefBlock)
                        call lsysbl_releaseVector (rrhsBlock)
                        call lsysbl_releaseMatrix (rmatrixBlock)
!                         call lsyssc_releaseMatrix (rlumpedmass)
                        print*,"-------------------"
                        print*, k-1," iterations needed "
                        print *, " residuum : " , defect
                        exit
                end if
            end if

            ! STEP 2.6: Solve the system
            !
            ! During the linear solver, the boundary conditions are also
            ! frequently imposed to the vectors. But as the linear solver
            ! does not work with the actual solution vectors but with
            ! defect vectors instead.
            RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
            ! Create a BiCGStab-solver. Attach the above filter chain
            ! to the solver, so that the solver automatically filters
            ! the vector during the solution process.
            p_RfilterChain => RfilterChain
            nullify(p_rpreconditioner_cells)
            call linsol_initBiCGStab (p_rsolverNode_cells,p_rpreconditioner_cells,p_RfilterChain)
            !call linsol_initGMRES (p_rsolverNode_cells,64,p_rpreconditioner_cells,p_RfilterChain)

            ! Attach the system matrix to the solver.
            Rmatrices = (/rmatrixBlock/)

            call linsol_setMatrices(p_RsolverNode_cells,Rmatrices)
            
            ! Initialise structure/data of the solver. This allows the
            ! solver to allocate memory / perform some precalculation
            ! to the problem.
            call linsol_initStructure (p_rsolverNode_cells, ierror)
            if (ierror .ne. LINSOL_ERR_NOERROR) stop
            call linsol_initData (p_rsolverNode_cells, ierror)
            if (ierror .ne. LINSOL_ERR_NOERROR) stop
            ! Finally solve the system.
            ! Set the output level of the solver to 2 for some output
            p_rsolverNode_cells%ioutputLevel = 2
            p_rsolverNode_cells%depsRel=1E-15_DP
            p_rsolverNode_cells%nminIterations=1
            p_rsolverNode_cells%nmaxIterations=1000
            call linsol_precondDefect( p_rsolverNode_cells , rdefBlock )
            
            ! Store the iterationstats
            if ( p_rsolverNode_cells%iiterations >=  iteration_u_max ) then
                iteration_u_max = p_rsolverNode_cells%iiterations
            end if
            if ( p_rsolverNode_cells%iiterations <= iteration_u_min ) then
                iteration_u_min = p_rsolverNode_cells%iiterations
            end if
            u_average_local = u_average_local + p_rsolverNode_cells%iiterations

            ! compute the NEW  celldensity vector by adding rdef
            call lsysbl_vectorLinearComb ( rdefBlock , rcellBlock , 1.0_DP , 1.0_DP )
            
            ! Release solver data and structure
            call linsol_doneData (p_rsolverNode_cells)
            call linsol_doneStructure (p_rsolverNode_cells)

            ! Release the solver node and all subnodes attached to it (if at all):
            call linsol_releaseSolver (p_rsolverNode_cells)

            ! Cells Releasements
            call lsysbl_releaseVector (rrhsBlock)
            call lsysbl_releaseMatrix (rmatrixBlock)
            call lsysbl_releaseVector (rdefBlock)
!             call lsyssc_releaseMatrix (rlumpedmass)
                        
            ! test block for the linearity
            if(k.eq.3) then
                print *,' problem is nonlinear!!! '
                exit
            end if
            
        end do def_correction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!   end nonlinear loop   !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Calculate the global iteration_u_average
        iteration_u_average = iteration_u_average + ( u_average_local / ( k-1 ) )

        ! Store the iterationstats
        if ( k-1 >=  iteration_defcorr_max ) then
        iteration_defcorr_max = k-1
        end if
        if ( k-1 <= iteration_defcorr_min ) then
            iteration_defcorr_min = k-1
        end if
        iteration_defcorr_average = iteration_defcorr_average + k-1

        ! Applying the FCT code
        !call chemo_fct_lim ( rcell, kedge, dedge, nedge, &
        !                     rtempLump, rlaplace, rK, dtstep )

        deallocate ( kedge )
        deallocate ( dedge )

        call lsyssc_releaseMatrix (rK)
!         call lsyssc_releaseMatrix (rlumpedmass)
        call lsyssc_releaseMatrix (rmatrix)
        call lsyssc_releaseMatrix (rtempLinComb)
        
        print *,">>end of the subroutine chemo_defcorr ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_defcorr subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################


    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_artdiff subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine adds some artificial diffusion to the matrix rmatrix
    subroutine chemo_artdiff( rmatrix, dedge, kedge, nedge )

        ! The matrix which is to be modified
        type(t_matrixScalar) , intent(INOUT) :: rmatrix
        
        real(DP), dimension (:), intent(OUT) :: dedge
        integer, dimension(:,:), intent(OUT) :: kedge
        integer, intent(INOUT) :: nedge

        ! some local variables
        integer :: i, j, nvt, iedge
        integer :: ij_loc , ii_loc , ji_loc , jj_loc
        real(DP) :: d_ij, dk_ij, dk_ji
        integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
        integer(PREC_MATIDX), dimension(:), pointer :: p_Kld , p_Kdiagonal
        real(DP), dimension(:), pointer :: p_Da_mass, p_Da
        integer , dimension(:) , allocatable :: iaux

        print *,">>start of the subroutine chemo_artdiff ... "
        
        ! Get the structure (rmassmatrix and rmatrix have the same structure)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        call lsyssc_getbase_double (rmatrix,p_Da)
        nvt = rmatrix%NEQ

        allocate ( iaux( nvt ) )

        ! (!matrix style 9) Correction of the discrete convection operator
        iedge = 0
        iaux (:) = p_Kld(:)
        DO i = 1,nvt
            ii_loc = p_Kdiagonal ( i )

            DO ij_loc = p_Kdiagonal(i)+1, p_Kld(i+1)-1
            j = p_Kcol ( ij_loc )
            jj_loc = p_Kdiagonal ( j )
            ji_loc = iaux ( j )
            iaux ( j )=iaux ( j )+1

	    dk_ij=p_Da ( ij_loc )
	    dk_ji=p_Da ( ji_loc )
            
            ! Artificial diffusion coefficient
            d_ij = MAX( -dk_ij, 0.0_DP, -dk_ji )

            ! Elimination of negative entries
            p_Da ( ij_loc ) = p_Da ( ij_loc ) + d_ij
            p_Da ( ii_loc ) = p_Da ( ii_loc ) - d_ij
            p_Da ( ji_loc ) = p_Da ( ji_loc ) + d_ij
            p_Da ( jj_loc ) = p_Da ( jj_loc ) - d_ij
            
            iedge = iedge + 1
	    if( dk_ji.gt.dk_ij ) then
              kedge ( 1, iedge ) = i
              kedge ( 2, iedge ) = j
	      dedge ( iedge ) = MIN(d_ij, dk_ji+d_ij)
	    else
              kedge ( 1, iedge ) = j
              kedge ( 2, iedge ) = i
	      dedge ( iedge ) = MIN(d_ij, dk_ij+d_ij)
	    endif

            END DO
        END DO

        nedge = iedge

        deallocate ( iaux )

        print *,">>end of the subroutine chemo_artdiff ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_artdiff subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################

    !##################################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_fct_lim subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here we set up the antidiffusive fluxes as well as the diff.limiters
    subroutine chemo_fct_lim ( rvector, kedge, dedge, nedge, &
                               rK, dtstep, rdefVector )

        type(t_vectorScalar), intent(INOUT) :: rvector
        type(t_matrixScalar), intent(IN) :: rK
        real(DP), dimension(:), intent(IN) :: dedge
        integer, dimension(:,:), intent(IN) :: kedge
        integer, intent(IN) :: nedge
        real(DP) :: dtstep
        type(t_vectorScalar), intent(INOUT) :: rdefVector

        ! Some local variables
	integer :: i, j, iloc, ii_loc
        integer :: iedge, nvt
        real(DP) :: eps
        real(DP), dimension (: ), allocatable  :: rm, rp
        real(DP) :: pp, pm, qp, qm
        real(DP), dimension(:), pointer :: p_vectorentries
        real(DP), dimension(:), pointer :: p_defVectorentries

        integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
        integer(PREC_MATIDX), dimension(:), pointer :: p_Kld , p_Kdiagonal
        real(DP), dimension(:), pointer :: p_Da

        real(DP) :: daux
        nvt = rK%NEQ

        ! Get the structure (rmassmatrix and rmatrix have the same structure)
        call lsyssc_getbase_Kld (rK,p_Kld)
        call lsyssc_getbase_Kcol (rK,p_Kcol)
        call lsyssc_getbase_Kdiagonal(rK, p_Kdiagonal)
        call lsyssc_getbase_double (rK,p_Da)

        nvt = rK%NEQ
        eps = 1e-16

        print *,">>start of the subroutine chemo_fct_lim ... "

        allocate ( rp (nvt) )
        allocate ( rm (nvt) )

        ! Initialize rp, rm with zeros
        call chemo_initAuxArrays ( rp, rm )
        call lsyssc_getbase_double( rvector, p_vectorentries )
        call lsyssc_getbase_double( rdefVector, p_defVectorentries )

        !Constructing fluxes
        DO iedge = 1,nedge
	    i = kedge (1,iedge)
	    j = kedge (2,iedge)

            daux = dedge(iedge)*(p_vectorentries(i)-p_vectorentries(j))

	    rp(i)=rp(i)+MAX(0.0_DP,daux)
	    rm(i)=rm(i)+MIN(0.0_DP,daux)
        END DO

        DO i = 1,nvt
            
	    pp=rp(i)
	    pm=rm(i)
	
	    qp=0.0_DP
	    qm=0.0_DP

	    rp(i)=0.0_DP
	    rm(i)=0.0_DP

	    !ii_loc = p_Kdiagonal ( i )
            DO iloc = p_Kld(i)+1, p_Kld(i+1)-1

		!matrix style 9 => skip the diagonal entry
		if( iloc.eq.p_Kdiagonal (i) ) then
		  cycle
		endif

            	j = p_Kcol ( iloc )
		daux = p_Da ( iloc )*(p_vectorentries(j)-p_vectorentries(i))

		qp=qp+MAX(0.0_DP, daux)
		qm=qm+MIN(0.0_DP, daux)
	    END DO

            if ( pp > eps ) then
		rp(i) = MIN ( 1.0_DP, qp/pp )
	    endif
            if ( pm < -eps ) then
		rm(i) = MIN ( 1.0_DP, qm/pm )
	    endif
	    	
	END DO

        ! Correction of the low-order solution
        DO iedge = 1 , nedge

            ! Node numbers for the current edge
            i = kedge ( 1 , iedge )
            j = kedge ( 2 , iedge )

            ! Antidiffusive flux to be limeted
            daux = dedge(iedge)*(p_vectorentries(i)-p_vectorentries(j))

            IF ( daux > 0 ) THEN
                daux = dtstep*rp(i)*daux
            ELSE
		daux = dtstep*rm(i)*daux
            END IF
	
	    !high-order resolution
	    p_defVectorentries ( i ) = p_defVectorentries ( i ) + daux
	    p_defVectorentries ( j ) = p_defVectorentries ( j ) - daux
        END DO

        ! Cleaning up the memory
        deallocate ( rm )
        deallocate ( rp )
        
        print *,">>end of the subroutine chemo_fct_lim ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_fct_lim subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   begin of chemo_initAuxArrays subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This routine sets all entries of the arrays to ZERO
    subroutine chemo_initAuxArrays ( rp, rm )

        real(DP), dimension (:), intent(INOUT) :: rp, rm

        ! locla variables
        integer :: i, na

        print *,">>start of the subroutine chemo_initAuxArrays ... "
        
        na = size ( rp )

        DO i = 1, na
	    rp ( i ) = 0.0_DP
	    rm ( i ) = 0.0_DP
        END DO

        print *,">>end of the subroutine chemo_initAuxArrays ... "
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!   end of chemo_initAuxArrays subroutine  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !##################################################################################


    !##################################################################################
    subroutine chemo_readGMV (gmvFile, rvector1, rvector2)
    !<input>

        ! The filename of the file to read.
        character(LEN=*), intent(IN) :: gmvFile

        ! Here the data is being written in
        ! rvector1 correspond to the cell-densities (u)
        ! rvector2 correspond to the chemoattractant-concentration (c)
        type(t_vectorScalar), intent(INOUT) :: rvector1, rvector2

    !</input>

        ! local variables
        real(DP), dimension(:), pointer :: p_vectorentries1, p_vectorentries2
        integer :: ientry, iunit, ios, isbuflen, ilinenum, found, iNEQ, iNODES
        character(LEN=PARLST_LENLINEBUF) :: sdata
        character(LEN=PARLST_LENLINEBUF) :: sbuf
        real(DP) :: dvalue

        print *,">>start of the subroutine chemo_readGMV ... "

        ! get the pointers to the vector entries
        call lsyssc_getbase_double( rvector1 , p_vectorentries1 )
        call lsyssc_getbase_double( rvector2 , p_vectorentries2 )
        iNEQ = rvector1%NEQ

        ! open the file (the gmv) for reading
        call io_openFileForReading(gmvFile, iunit)

        ! Oops...
        if (iunit .eq. -1) then
        call output_line ('Error opening .gmvfile: '//trim(gmvFile), &
            OU_CLASS_ERROR,OU_MODE_STD,'chemo_readGMV')
        call sys_halt()
        end if

        ! Read all lines from the file
        ios = 0
        ilinenum = 0
        found = 0
        ientry = 1
    
            do while (ios .eq. 0)
        
                ! Read a line from the file into sbuf
                call chemo_readlinefromfile (iunit, sdata, isbuflen, ios)
                ilinenum = ilinenum + 1
    
    !             if (sdata .eq.'') then
    !                 print*,"NO DATA IN FILE"
    !                 return
    !             end if
    
                sbuf = adjustl(sdata)
    
                if ( len(sbuf) >= 5 .AND. sbuf(1:5) .eq. "nodes" ) then
                    sbuf = adjustl(sbuf(6:len(sbuf)))
                    read(sbuf,*) iNODES
                    if (iNODES /= iNEQ) then
                        ! vector and gmvFile doesn't correspond
                        ! vector does not have the same amount of entries than
                        ! the gmv-file has data
                        call output_line ('Error writing vector : ', &
                                    OU_CLASS_ERROR,OU_MODE_STD,'chemo_readGMV')
                        call sys_halt()
                    end if
                else
                    ! read till we got the right lines
                    if (len(sbuf) >= 8 .AND.sbuf(1:8) .eq. "variable") then
                        ! the following lines determine the datas (e.g. the values at the nodes)
                        found = 1
                        exit
                    end if
                end if
    
    
    
            end do
    
            if (found == 1) then
    
                found = 0
                ios = 0
                ! skip one line, for this is just the name of the following datas
                ! as we restrict our gmvs to consist of the cell and chemoattractant
                ! data, we know the structure of the data in the gmv file
                call chemo_readlinefromfile (iunit, sdata, isbuflen, ios)
                ilinenum = ilinenum + 1
                do while (ios .eq. 0)
        
                    ! Read a line from the file into sbuf
                    call chemo_readlinefromfile (iunit, sdata, isbuflen, ios)
                    ilinenum = ilinenum + 1
        
                    if (sdata .eq.'') then
                        print*,"NO DATA IN FILE"
                        return
                    end if
    
                    sbuf = adjustl(sdata)
                    if (.NOT. (len(sbuf) >= 7 .AND. sbuf(1:7).eq."endvars")) then
                        ! convert the string to the double precision vector entry
    !                     dvalue = sys_Str2Double(sbuf,"(ES27.20E3)")
                        dvalue = sys_Str2Double(sbuf,"(ES16.8E3)")
                        ! and put it into the vector
                        p_vectorentries1 ( ientry ) = dvalue
                        ientry = ientry+1
                    else
                        ! the data for the cells are finished
                        ! continue with the ones of the chemoattractant
                        if (ientry /= iNEQ+1) then
                            ! the vector was not filled up completely
                            print*,"Something went wrong here"
                        end if
                        call chemo_readlinefromfile (iunit, sdata, isbuflen, ios)
                        call chemo_readlinefromfile (iunit, sdata, isbuflen, ios)
                        ilinenum = ilinenum + 2
                        found = 1
                        exit
                    end if
        
                end do
    
                        if (found == 1) then
    
                            found = 0
                            ios = 0
                            ientry = 1
                            do while (ios .eq. 0)
                    
                                ! Read a line from the file into sbuf
                                call chemo_readlinefromfile (iunit, sdata, isbuflen, ios)
                                ilinenum = ilinenum + 1
                    
                                if (sdata .eq.'') then
                                    print*,"NO DATA IN FILE"
                                    return
                                end if
                
                                sbuf = adjustl(sdata)
                                if (.NOT. (len(sbuf) >= 7 .AND. sbuf(1:7).eq."endvars")) then
                                    ! convert the string to the double precision vector entry
    !                                 dvalue = sys_Str2Double(sbuf,"(ES27.20E3)")
                                    dvalue = sys_Str2Double(sbuf,"(ES16.8E3)")
                                    ! and put it into the vector
                                    p_vectorentries2 ( ientry ) = dvalue
                                    ientry = ientry+1
                                else
                                    if (ientry /= iNEQ+1) then
                                        ! the vector was not filled up completely
                                        print*,"Something went wrong here"
                                    end if
                                    ! all data has been read
                                    exit
                                end if

                            end do

                        end if

            end if

            ! Close the file.
            close (iunit)

            print *,">>end of the subroutine chemo_readGMV ... "

        end subroutine

    !##################################################################################

        subroutine chemo_readlinefromfile (iunit, sdata, ilinelen, ios)

        ! The unit where to read from; must be connected to a file.
        integer, intent(in) :: iunit

        ! The string where to write data to
        character(LEN=*), intent(out) :: sdata
        
        ! Length of the output
        integer, intent(out) :: ilinelen
        
        ! Status of the reading process. Set to a value <> 0 if the end
        ! of the file is reached.
        integer, intent(out) :: ios
        
        ! local variables
        integer :: eol,isize
        character :: c
        
        sdata = ''
        ilinelen = 0
        isize = 0
        
        ! Read the data - as long as the line/file does not end.
        eol = NO
        ios = 0
        do while ((ios .eq. 0) .and. (eol .eq. NO))
            
            ! Read a character.
            ! Unfortunately, Fortran forces me to use this dirty GOTO
            ! to decide processor-independently whether the line or
            ! the record ends.
            read (unit=iunit,fmt='(A1)',iostat=ios,advance='NO', end=10, eor=20,size=isize) c
            goto 30
            
        10  continue
            ! End of file.
            ios = -1
            goto 30
            
        20  continue
            ! End of record = END OF LINE.
            eol = YES
        
            ! Set error flag back to 0.
            ios = 0
            
        30  continue
            ! Do not do anything in case of an error
            if ((ios .eq. 0) .and. (isize .ge. 1)) then
            
            ilinelen = ilinelen + 1
            sdata (ilinelen:ilinelen) = c
            
            end if
        
        end do
        
        end subroutine


end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   end of chemotaxis_pattern_FCT module   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!##############################################################################
!# ****************************************************************************
!# <name> chemotaxis_pattern_FCT </name>
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

module chemotaxis_pattern_FCT
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

  use chemotaxis_callback
  
  implicit none

contains

  ! ***************************************************************************

	function Fsol(x,y,a,b) result (f_result1)					!----> fct for the initial condition
		implicit none								!	of the solutionvector
		real(DP) :: f_result1, x, y, a, b					!	(centered)

		f_result1 = a*dexp(-b*((x-0.5_DP)*(x-0.5_DP)+(y-0.5_DP)*(y-0.5_DP)))
	end function Fsol

	function Fchemo(x,y,a,b) result (f_result2)					!----> fct for the initial condition
		implicit none								!	of the chemoattractant
		real(DP) :: f_result2, x, y, a, b					!	(centered)

		f_result2 = a*dexp(-b*((x-0.5_DP)*(x-0.5_DP)+(y-0.5_DP)*(y-0.5_DP)))
	end function Fchemo
        


!<subroutine>

  subroutine chemotaxispatternFCT
  
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
    
    ! A trilinear , bilinear and linear form describing the analytic problem to solve
    type(t_trilinearForm) :: rtriform
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar)  :: rmatrix, rmatrixchemo, rsysmatrix , rlaplace, rtemp
    type(t_matrixScalar)  :: rmassmatrix
    type(t_vectorScalar)  :: rrhscell,rcell
    type(t_vectorScalar)  :: rrhschemo,rchemoattract,cbvector_temp&
                                    , rinitcells , rinitchemo , rdef
    type(t_vectorScalar) :: ranalyticcells , ranalyticchemo, rcold, ruold

    !A pointer to the entries of vector rchemoattract				!----->added to access the chemoattractant
    real(DP), dimension(:), pointer ::  p_vectordata, p_uold, p_cold,&
                                                     p_analyticcells, p_analyticchemo, p_chemodata

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock, rmatrixBlockchemo
    type(t_vectorBlock) :: rcellBlock,rvectorBlockchemo,rrhsBlock,rrhsBlockchemo,rtempBlock, rdefBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC, rdiscreteBCchemo

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner,p_rsolverNode_cells,p_rpreconditioner_cells

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
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

    ! Time and time step counter
    real(DP) :: dtime
    integer :: itimestep

    ! maximal iterations for the defect correction . Since we allow possible nonlinearities
    integer :: maxiterationdef

    ! Variables for the defect correction ( nonlinear loop )
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
    integer :: gmvfolder, checkneg
              
    ! error analysis
    real(DP) :: uerror, cerror, tol

    ! error-norm constant
    ! To be set in the .dat file
    integer :: CTRLNORM
            
    ! Whether or not taking an const initial sol
    ! =0 means const initial sol
    ! =1 means exponential initial sol
    ! Whether or not we set the IC by l_2 projection ( = 1 ) or by interpolation ( =0 )
    integer :: INITSOL , L2PROJ

    ! If INITSOL = 0 , meaning const init sol, we use these variables to
    ! define the const value
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

    ! Ok, let's start.
    !
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
      
    ! We set the chemotaxis and diffusion coefficient
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'CHI', CHI, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'W', W, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'BETA', BETA, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'SIGMA', SIGMA, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'ALPHA', ALPHA, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'R', R, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'D_1', D_1, 1.0_DP)
    D_1 = 0.0625_DP ! This is the diffusioncoefficient in the paper
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'PHI', PHI, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'GAMMA', GAMMA, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'N', N, 1.0_DP)
    call parlst_getvalue_double (rparams, 'COEFFICIENTS', 'D_2', D_2, 1.0_DP)
    D_2 = 1.0_DP  ! This is the diffusioncoefficient in the paper
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

    ! Get the tolerance threshold value for steady state issues
    call parlst_getvalue_double (rparams, 'ERROR', 'TOL', tol, 0.0001_DP)
    

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/patternQUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, './pre/patternQUAD.tri', rboundary)
    
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
                                   
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix
    ! We create the matrices and vectors for the pb
                 
    ! Get some pointers  for the errorctrl
      call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),ruold,.true.)
      call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rcold,.true.)
      call lsyssc_getbase_double(rcold,p_cold)
      call lsyssc_getbase_double(ruold,p_uold)
                   
    call chemo_creatematvec ( rmassmatrix , rsysmatrix, rlaplace , rmatrixchemo ,&
                                             rchemoattract ,rrhschemo , rcell , rrhscell , rdef , rdiscretisation )
                                     
    ! We construct some matrices in advance, so the actual loop computations
    ! will be much faster
    ! This routine is defined in this file
    call chemo_initmat ( rmassmatrix, rsysmatrix, rlaplace , dtstep , D_1 , D_2 )
        
    ! Use the discretisation to create the IC for the two solution vectors  rcell and chemoattract.
    call chemo_initIC ( rchemoattract , rcell , p_vectordata , rtriangulation , rmassmatrix , L2PROJ ,&
                                 INITSOL , U_0 , C_0 , A_CELLS , B_CELLS , A_CHEMO , B_CHEMO)


    ! Setting the boundary conditions

      ! Now we have the raw problem. What is missing is the definition of the boudary
      ! conditions.
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
    call chemo_initBC ( rdiscreteBCchemo , rboundary , rdiscretisation , C_0 , U_0 , INITSOL)




    ! Setting the Blockvectors

    ! The linear solver only works for block matrices/vectors - but above,
    ! we created scalar ones. So the next step is to make a 1x1 block
    ! system from the matrices/vectors above which the linear solver
    ! understands.
    call lsysbl_createMatFromScalar (rsysmatrix,rmatrixBlockchemo,rdiscretisation)
    call lsysbl_createVecFromScalar (rchemoattract,rvectorBlockchemo,rdiscretisation)

    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    rmatrixBlockchemo%p_rdiscreteBC => rdiscreteBCchemo
    rrhsBlockchemo%p_rdiscreteBC => rdiscreteBCchemo
    rvectorBlockchemo%p_rdiscreteBC => rdiscreteBCchemo

      


    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call matfil_discreteBC (rmatrixBlockchemo)

    ! Setting some resolver issues
      ! During the linear solver, the boundary conditions are also
      ! frequently imposed to the vectors. But as the linear solver
      ! does not work with the actual solution vectors but with
      ! defect vectors instead.
      ! So, set up a filter chain that filters the defect vector
      ! during the solution process to implement discrete boundary conditions.
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

      ! Create a BiCGStab-solver. Attach the above filter chain
      ! to the solver, so that the solver automatically filters
      ! the vector during the solution process.
      p_RfilterChain => RfilterChain
      nullify(p_rpreconditioner)
      call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
      
      ! Set the output level of the solver to 2 for some output
      p_rsolverNode%ioutputLevel = 0
      
      ! Attach the system matrix to the solver.
      ! First create an array with the matrix data (on all levels, but we
      ! only have one level here), then call the initialisation
      ! routine to attach all these matrices.
      ! Remark: Don't make a call like
      !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
      ! This doesn't work on all compilers, since the compiler would have
      ! to create a temp array on the stack - which does not always work!
      Rmatrices = (/rmatrixBlockchemo/)
      call linsol_setMatrices(p_RsolverNode,Rmatrices)
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop

!###########################################################################
    ! Setting some preliminaries for the calculations of the cell density
      call lsysbl_createVecFromScalar (rcell,rcellBlock,rdiscretisation)
      call lsysbl_createVecFromScalar (rdef,rdefBlock,rdiscretisation)

      ! Now we have the raw problem. What is missing is the definition of the boudary
      ! conditions.
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
    call chemo_initBC ( rdiscreteBC , rboundary , rdiscretisation , C_0 , U_0 , INITSOL)


      rcellBlock%p_rdiscreteBC => rdiscreteBC
    rdefBlock%p_rdiscreteBC => rdiscreteBC

    ! printing out the initial conditions into a gmv_file
    select case ( gmvfolder )

        case ( 0 )
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                    'gmvcpld/solution_init.gmv.pattern')
            call lsyssc_getbase_double (rcell,p_Ddata)
            call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
            call lsyssc_getbase_double (rchemoattract,p_Ddata)
            call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)

        case default
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                    'gmvcpld2/solution_init.gmv.pattern')
            call lsyssc_getbase_double (rcell,p_Ddata)
            call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
            call lsyssc_getbase_double (rchemoattract,p_Ddata)
            call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)

    end select
      

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

    ! to determine negativ solution values.
    quit =.false.
    ! We will now start assembling the solutionvectors by first computing a new value of
    ! the chemoattractant vector, e.g. c_n --->c_{n+1}. Then we use this updated value for the
    ! attractant to compute the actual solution of our problem,e.g. the solution vector u_{n+1}
    ! ( u_n, c_{n+1} ---> u_{n+1} ). We will use in both cases am implicit Euler method.
    !
    ! STEP 1 : compute c
    ! STEP 2 : compute u


    if(steps /= 0) then
        ntimesteps = steps
    end if

    ! Start the timeloop
    timeloop : do itimestep=1,ntimesteps
       
    ! time stats
    call cpu_time(time_start)

      ! Next time step.
      dtime = dtime + dtstep
      
                      
    ! Store the old sols for the errorctrl
    call lsyssc_getbase_double( rcell , p_vectordata )
    call lsyssc_getbase_double( rchemoattract , p_chemodata)
    call lalg_copyVector (p_vectordata,p_uold)
    call lalg_copyVector (p_chemodata,p_cold)
                   
      ! STEP 1.1: Form the right hand side for c:
      ! It consists of M c_n +dt * ( u_{n} , phi )
        call chemo_initrhsC (rrhsBlockchemo , rchemoattract ,rrhschemo , rcell , rdiscretisation , rmassmatrix , dtstep, PHI )

    call lsysbl_createVecFromScalar (rrhschemo,rrhsBlockchemo,rdiscretisation)
      ! Next step is to implement boundary conditions into the RHS,
      ! solution and matrix. This is done using a vector/matrix filter
      ! for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! vectors/matrix. Call the appropriate vector/matrix filter that
      ! modifies the vectors/matrix according to the boundary conditions.
      call vecfil_discreteBCrhs (rrhsBlockchemo)
      call vecfil_discreteBCsol (rvectorBlockchemo)

      
      ! STEP 6: Solve the system
      !
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    !CALL lsysbl_createVecBlockIndirect (rrhsBlock, rcellBlock, .TRUE.)
    call lsysbl_createVecBlockIndirect (rrhsBlockchemo, rtempBlock, .false.)
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      call linsol_solveAdaptively (p_rsolverNode,rvectorBlockchemo,rrhsBlockchemo,rtempBlock)

    ! Store the iterationstats
    if ( p_rsolverNode%iiterations >=  iteration_c_max ) then
        iteration_c_max = p_rsolverNode%iiterations
    end if
    if ( p_rsolverNode%iiterations <= iteration_c_min ) then
        iteration_c_min = p_rsolverNode%iiterations
    end if
    iteration_c_average = iteration_c_average + p_rsolverNode%iiterations

      ! STEP 7: Postprocessing
      !
      ! That's it, rcellBlock now contains our solution. We can now
      ! start the postprocessing.

!    if(OUTPUT.eq.1)then
!      ! Start UCD export to GMV file:
!    select case  (gmvfolder)
!
!        case (0)
!        	call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                        'gmvcpld/chemoattract.gmv.'//trim(sys_si0L(itimestep,5)))
!        case default
!                call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                        'gmvcpld2/chemoattract.gmv.'//trim(sys_si0L(itimestep,5)))
!
!    end select
!
!      	call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!      	call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!
!      	! Write the file to disc, that's it.
!      	call ucd_write (rexport)
!      	call ucd_release (rexport)
!
!    end if
               
    !     We  compute the norm of the difference between two successive timesteps. This will be used as a basic error control, coming up below...
    !     The same will be done for the celldensity function
    cerror =  lalg_errorNorm (p_cold,p_chemodata, CTRLNORM)
!     print*,"cerror/dtstep = ", cerror/dtstep
                
      ! Release the block matrix/vectors
      call lsysbl_releaseVector (rtempBlock)
    call lsyssc_releaseVector (rrhschemo)













            ! STEP 2.1: Form the right hand side for c:  @@@@@!!!!M u_{old}!!!!!@@@@@
            ! It consists of M u_n
!             call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rrhs,.true.)
            !!!!!!!!!#########-->there should rest a instantiated rrhs for the MatVecMultiplication
            
            call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rtemp)
            call lsyssc_copyMatrix(rmassmatrix, rtemp)
            call lsyssc_lumpMatrixScalar (rtemp,LSYSSC_LUMP_DIAG,.false.)
                        
            ! Now form the actual RHS by matrix vector multiplication!
            ! rmatrix * u_{old}
            ! The massmatrix has already been built outside the loop !!!
            call lsyssc_scalarMatVec(rtemp,rcell,rrhscell,1.0_DP,0.0_DP)		!constructs the RHS


!             call collct_init (rcollection)
!             call collct_setvalue_vecsca (rcollection, 'cbvector', rcell, .true.)
!             rlinform%itermCount = 1
!             rlinform%Idescriptors(1) = DER_FUNC
!             call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1), &
!                                         rlinform, .false., rrhscell, coeff_pattern_RHS_u, rcollection)
!
!             ! Release the collection structure
!             call collct_done(rcollection)

            call lsyssc_releaseMatrix(rtemp)
            
    ! Setting up the cb collection to hand over the params
! RS : no need for f, since it doesnt exist
!             call collct_init (rcollection)
!             rcollection%DquickAccess(1) = dtstep
!             rcollection%DquickAccess(2) = CHI
!             rcollection%DquickAccess(3) = SCALE_C
!             rcollection%DquickAccess(4) = SCALE_U
!
!             rlinform%itermCount = 1
!             rlinform%Idescriptors(1) = DER_FUNC
!             call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
!                                         rlinform,.false.,rrhscell,coeff_cherkurRHSf,rcollection)
!             call collct_done(rcollection)
                                                                                        
            ! Now we should have the correct RHS = rhs like its mentioned above.
    ! Here the nonlinearity of CHI comes into account
    ! We calculate some iteratives to suit the residual tolerance
    print*,"----------------"
    print*,"timestep : ", itimestep
    print*,"----------------"
    call  chemo_defcorr(  rcell, rchemoattract, rcellBlock, rmassmatrix, rlaplace, &
                       rrhscell, rdiscretisation, rdiscreteBC, dtstep, D_1, CHI, ALPHA, r, GAMMA, N, maxiterationdef, defectTol ,&
                        iteration_u_max, iteration_u_min, iteration_u_average, iteration_defcorr_max, &
                        iteration_defcorr_min, iteration_defcorr_average)


    
      ! STEP 2.7: Postprocessing
      !
      ! That's it, rcellBlock now contains our solution. We can now
      ! start the postprocessing.
    if(OUTPUT .eq. 1) then
      ! Start UCD export to GMV file:
    select case (gmvfolder)

        case (0)
          call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvcpld/solution.gmv.pattern.'//trim(sys_si0L(itimestep,5)))

        case default
             call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvcpld2/solution.gmv.pattern.'//trim(sys_si0L(itimestep,5)))
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
                            'gmvcpld/solution.gmv.pattern.'//trim(sys_si0L(itimestep,5)))
    
            case default
                call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                            'gmvcpld2/solution.gmv.pattern.'//trim(sys_si0L(itimestep,5)))
        end  select
    
        call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
        call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
        
        ! Write the file to disc, that's it.
        call ucd_write (rexport)
        call ucd_release (rexport)

        end if
    end if
      
           
    ! We  compute the norm of the difference between two successive timesteps
    ! (as mentioned above). This will be used as a basic error control, coming up below...
    uerror = lalg_errorNorm (p_uold,p_vectordata, CTRLNORM )
!   print*,"uerror/dtstep",uerror/dtstep


    ! time stats
    call cpu_time(time_stop)
    time_accum = time_accum+( time_stop-time_start )

    if ( time_stop-time_start .ge. time_max) then
        time_max = time_stop-time_start
    end if

    if ( time_stop-time_start .le. time_min) then
        time_min = time_stop-time_start
    end if
    ! If we want to export just the last vectors
!    if ( OUTPUT .eq. 2 .AND. itimestep.eq. steps) then
!        select case (gmvfolder)
!
!        case (0)
!            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                    'gmvcpld/cellsEND.gmv.'//trim(sys_si0L(itimestep,5)))
!            call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
!            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!            ! Write the file to disc, that's it.
!            call ucd_write (rexport)
!            call ucd_release (rexport)
!
!            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                    'gmvcpld/chemoattractEND.gmv.'//trim(sys_si0L(itimestep,5)))
!            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!            ! Write the file to disc, that's it.
!            call ucd_write (rexport)
!            call ucd_release (rexport)
!
!        case default
!            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                    'gmvcpld2/cellsEND.gmv.'//trim(sys_si0L(itimestep,5)))
!            call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
!            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!            ! Write the file to disc, that's it.
!            call ucd_write (rexport)
!            call ucd_release (rexport)
!
!            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                    'gmvcpld2/chemoattractEND.gmv.'//trim(sys_si0L(itimestep,5)))
!            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!            ! Write the file to disc, that's it.
!            call ucd_write (rexport)
!            call ucd_release (rexport)
!
!    end select
!    else
    ! If we compute a negative solution, we'll export the actual negativ
    ! solutions andexit the computations
   if ( checkneg .eq. 1) then
   checkneg_loop :  do i=1,rtriangulation%NVT
        if( p_vectordata(i) < negthres .AND. .NOT.quit) then
             select case (gmvfolder)

                case (0)
                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvcpld/solutionEND.gmv.pattern.'//trim(sys_si0L(itimestep,5)))
                    call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
                     ! Write the file to disc, that's it.
                    call ucd_write (rexport)
                    call ucd_release (rexport)

!                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                        'gmvcpld/chemoattractEND.gmv.'//trim(sys_si0L(itimestep,5)))
!                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!                    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!                     ! Write the file to disc, that's it.
!                    call ucd_write (rexport)
!                    call ucd_release (rexport)
                case default
                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmvcpld2/solutionEND.gmv.pattern.'//trim(sys_si0L(itimestep,5)))
                    call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
                     ! Write the file to disc, that's it.
                    call ucd_write (rexport)
                    call ucd_release (rexport)

!                    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                        'gmvcpld2/chemoattractEND.gmv.'//trim(sys_si0L(itimestep,5)))
!                    call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!                    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!                     ! Write the file to disc, that's it.
!                    call ucd_write (rexport)
!                    call ucd_release (rexport)

            end  select
            quit = .true.
         end if
    end do checkneg_loop
    end if
!    end if
    ! Chemoattract Releasements
      call lsysbl_releaseVector (rrhsBlockchemo)

      ! Release the scalar matrix/rhs vector which were used to create
      ! the block matrices/vectors. These must exist as long as the
      ! block matrices/vectors exist, as the block matrices/vectors are
      ! only 'copies' of the scalar ones, sharing the same handles!


      ! Release the scalar matrix/rhs vector which were used to create
      ! the block matrices/vectors. These must exist as long as the
      ! block matrices/vectors exist, as the block matrices/vectors are
      ! only 'copies' of the scalar ones, sharing the same handles!
                                     
    ! The simulation should stop if we reached a nearly steady state. So thats whats done here...
    ! We approximate our time-derivative with the first order accurate difference quotient.
    ! d sol/ dt ~ (sol_n - sol_{n-1})/dtstep
        if(steps == 0) then
         if(tol >= uerror/dtstep .AND. tol >=cerror/dtstep)  then
                print *, "tolerance threshold reached-----simulation ends up in 'nearly' steady state."
                print *,"############ differnece to steady state in c ################"
                print *,cerror / dtstep
                print *,"############################"
                print *,"############## differnece to steady state in u ##############"
                print *, uerror / dtstep
                print *,"############################"
                quit = .true.
        end if
       end if
                           
    if( quit ) then
        EXIT
    end if
    end do timeloop
    ! If we want to export just the last vectors
    if ( OUTPUT .eq. 2 ) then
        select case (gmvfolder)

        case (0)
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                    'gmvcpld/solutionEND.gmv.pattern.'//trim(sys_si0L(itimestep-1,6)))
            call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
            call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)
    
!            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                    'gmvcpld/chemoattractEND.gmv.'//trim(sys_si0L(itimestep-1,5)))
!            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!            ! Write the file to disc, that's it.
!            call ucd_write (rexport)
!            call ucd_release (rexport)

        case default
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                    'gmvcpld2/solutionEND.gmv.pattern.'//trim(sys_si0L(itimestep-1,5)))
            call lsyssc_getbase_double (rcellBlock%RvectorBlock(1),p_Ddata)
            call ucd_addVariableVertexBased (rexport,'cells',UCD_VAR_STANDARD, p_Ddata)
            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
                    call ucd_addVariableVertexBased (rexport,'chemoattractant',UCD_VAR_STANDARD, p_Ddata)
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)
    
!            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                    'gmvcpld2/chemoattractEND.gmv.'//trim(sys_si0L(itimestep-1,5)))
!            call lsyssc_getbase_double (rvectorBlockchemo%RvectorBlock(1),p_Ddata)
!            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!            ! Write the file to disc, that's it.
!            call ucd_write (rexport)
!            call ucd_release (rexport)

    end select
    end if


    ! Again if we end up in a non steady state or if we are computing certain steps, we print
    !out the steady state error as well as the
    ! absolute error of the numerical solution:
    ! norm(sol_analytic - sol_numerical)
    ! Just printing out the params we use in our computational work
    print *," patternformation_FCT:"
    print *,"----------------------------"
    print *, "computed with :"
    print *,"Simulation parameters:"
    if(steps /= 0)then
        print*,"fixed amount of time steps = ",steps
    end if
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

    print*,"Initial conditions :"
    print*,"u ( x, y, 0 ) = 1"
    print*,"c ( x, y, 0 ) = 1 + 0.1*exp ( -10 ( ( x-1 )^2 + ( y-1 )^2 )  )"

    print *,"----------------------------"

    print *, "---------------------------"
    if ( quit ) then
        print *, "stopped @  timestep = ", itimestep," out of ",ntimesteps
    else
        print *, "stopped @  timestep = ", itimestep-1," out of ",ntimesteps
    end if
    print*, "----------------------------"
    print*, "--------TIMESTATS-----------"
    print*, "max time per iteration : ",time_max
    print*, "min time per iteration: ", time_min
    if(time_accum > 10800) then
        print*, "total time for iterations : ", time_accum/3600, " h"
    else
        print *, "total time for iterations :", time_accum/60, " mn"
    end if
    print*, "----------------------------"
    print*, "------ITERATION STATS-------"
    print*, "SOLVING FOR C :"
    print*,"max iteration needed :",iteration_c_max
    print *,"min iteration needed :", iteration_c_min
    if ( quit ) then
        print*, "average iterations needed :", iteration_c_average / itimestep
    else
        print*, "average iterations needed :", iteration_c_average / (itimestep-1)
    end if
    print*, "SOLVING FOR U :"
    print*,"max iteration needed :", iteration_u_max
    print *,"min iteration needed :", iteration_u_min
    if ( quit ) then
        print*, "average iterations needed :", iteration_u_average / itimestep
    else
        print*, "average iterations needed :", iteration_u_average / (itimestep-1)
    end if
    print*, "SOLVING FOR U ( DEF_CORR ) :"
    print*,"max iteration needed :", iteration_defcorr_max
    print *,"min iteration needed :", iteration_defcorr_min
    if ( quit ) then
        print*, "average iterations needed :", iteration_defcorr_average / itimestep
    else
        print*, "average iterations needed :", iteration_defcorr_average / (itimestep-1)
    end if
                      
 print *,"############ differnece to steady state in c ################"
                    print *,cerror / dtstep
                    print *,"############################"
                    print *,"############## differnece to steady state in u ##############"
                    print *, uerror / dtstep
                    print *,"############################"
                    
    print*, "----------------------------"

    ! Now, clean up so that all the memory is available again.
      !
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

    call lsyssc_releaseVector (rcell)
    call lsyssc_releaseVector (rchemoattract)
    call lsyssc_releaseVector (rdef)
    call lsyssc_releaseVector (rrhscell)
          
    call lsyssc_releaseVector (rcold)
    call lsyssc_releaseVector (ruold)
         


    ! Release parameterlist
     call parlst_done (rparams)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

!##################################################################################
!##################################################################################

! Here we assemble the matrices in advance ( to save computation time )
    subroutine chemo_initmat (rmassmatrix, rsysmatrix, rlaplace , dtstep , D_1 , D_2)

    ! This is where the matrices shoiuld be stored
    type(t_matrixScalar), intent(INOUT):: rmassmatrix , rlaplace , rsysmatrix

    real(DP) , intent (INOUT) :: dtstep , D_1 , D_2
    ! local bilinearform to construct the matrices
    type(t_bilinearForm) :: rform

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

    ! Since we want to save time, we use already constructed matrices
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = dtstep * 32.0_DP  ! This is the coefficient in the paper
    call bilf_buildMatrixScalar (rform,.false.,rsysmatrix)
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

    end subroutine

!##################################################################################
!##################################################################################

! In thie routine the basic matrices and vectors are initialized
    subroutine chemo_creatematvec ( rmassmatrix , rsysmatrix, rlaplace , rmatrixchemo ,&
                                                        rchemoattract , rrhschemo , rcell , rrhs , rdef , rdiscretisation )

    ! This is where the matrices shoiuld be stored
    type(t_matrixScalar) , intent(INOUT)  ::rmassmatrix, rsysmatrix, rlaplace , rmatrixchemo

    ! This are the vectors
    type(t_vectorScalar) , intent(INOUT) :: rchemoattract, rcell, rdef ,rrhschemo , rrhs

    ! The underlying discretisation
    type(t_blockDiscretisation) , intent (IN) :: rdiscretisation


    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rlaplace)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rmatrixchemo)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rsysmatrix)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rmassmatrix)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rchemoattract,.true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rcell,.true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rdef,.true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rrhschemo, .true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rrhs,.true.)

    end subroutine

!##################################################################################
!##################################################################################

! This is called to implement the initial conditions, e.g. to set the initial sol vectors
    subroutine chemo_initIC ( rchemoattract , rcell , p_vectordata , rtriangulation , rmassmatrix , L2PROJ ,&
                                             INITSOL , U_0 , C_0 , A_CELLS , B_CELLS , A_CHEMO , B_CHEMO )
    
    ! massmatrix to construct an analyt proj
    type(t_matrixScalar)  ,intent (INOUT) :: rmassmatrix

    ! sol vectors
    type(t_vectorScalar) , intent (INOUT) :: rchemoattract, rcell

    real(DP), dimension(:) , pointer , intent(INOUT):: p_vectordata
    ! Underlying triangulation to derive vector-coords
    type(t_triangulation) , intent (IN) :: rtriangulation

    ! Params to set the ICs
    real(DP) , intent (IN) ::U_0 , C_0 ,  A_CELLS , B_CELLS , A_CHEMO , B_CHEMO
    integer , intent (IN) :: L2PROJ , INITSOL

    ! A local collection
    type(t_collection) :: rcollection

    ! Some pointers
    real(DP) , dimension(:,:) , pointer :: p_DvertexCoords
    real(DP), dimension(:) , pointer :: p_rchemoattract

    ! integer loop
    integer :: i

    call lsyssc_getbase_double(rcell,p_vectordata)
    call lsyssc_getbase_double(rchemoattract,p_rchemoattract)
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords,p_DvertexCoords)

    select case (L2PROJ)

        case (0)

            select case (INITSOL)

                case (0)
                    do i=1,rtriangulation%NVT
                                p_vectordata(i) = U_0
                                p_rchemoattract(i) = C_0
                    end do

                case default

                    do i=1,rtriangulation%NVT
                            p_vectordata(i) = Fsol(p_DvertexCoords(1,i),p_DvertexCoords(2,i),A_CELLS,B_CELLS)
                            p_rchemoattract(i) = Fchemo(p_DvertexCoords(1,i),p_DvertexCoords(2,i),A_CHEMO,B_CHEMO)
                    end do

            end select

        case default
            ! An alternative way to prescribe the initiaLvectors is presented below
            ! We simply use the direkt l_2 projection to approximate the analytic given init. sol
            ! in contrast to calculating an interpolation of the analytic fct.( like it's done above )
                select case (INITSOL)

                case (0)
                    ! Now we're setting up our collection
                    call collct_init (rcollection)
                    ! Setting the IC for u
                    rcollection%DquickAccess(1) = U_0
                    call anprj_discrDirect ( rcell , coeff_anprj_ic_const , rcollection )
                    ! Setting the IC for c
                    rcollection%DquickAccess(1) = C_0
                    call anprj_discrDirect ( rchemoattract , coeff_anprj_ic_const , rcollection )

                    call collct_done (rcollection)

                case default

                    ! Now we're setting up our collection
                    call collct_init (rcollection)

                    ! Setting the IC for u
                    call anprj_discrDirect ( rcell , coeff_anprj_ic_pattern )

                    rcollection%DquickAccess(1) = 1.0_DP / 32.0_DP
                    ! Setting the IC for c
                    call anprj_discrDirect ( rchemoattract , coeff_anprj_ic_const, rcollection )

                    ! Releasing the collection
                    call collct_done (rcollection)

            end select


    end select

    end subroutine

!##################################################################################
!##################################################################################

! Initialize the boundary conditions (used for both, the cell density sol and the chemoattractant sol)
    subroutine chemo_initBC ( rdiscreteBC , rboundary , rdiscretisation , C_0 , U_0 , INITSOL)

    ! The BC_type which should be initialized
    type(t_discreteBC) , intent (INOUT) ::rdiscreteBC

    type(t_boundary) , intent (IN) :: rboundary

    ! Underlying discretisation
    type(t_blockDiscretisation) , intent (IN) ::rdiscretisation

    ! BCs if we take constant initsol
    real(DP) , intent (IN) :: C_0 , U_0

    integer , intent(IN) :: INITSOL

    ! local collection
    type(t_collection) :: rcollection

    type(t_boundaryRegion) :: rboundaryRegion

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
!     if(INITSOL .eq. 0) then
    
    call collct_init(rcollection)
        rcollection%DquickAccess(1) = 0.0_DP
        rcollection%DquickAccess(2) = 0.0_DP
      
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC's for the first (here: one and only) component in the
      !   solution vector.
      ! - Discretise the boundary condition so that the BC's can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC's to rdiscreteBC for later use.
!       call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                        rboundaryRegion,rdiscreteBC,&
!                                        getBoundaryValues_constC,rcollection)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!       call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_constC,rcollection)
                               									
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!       call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                        rboundaryRegion,rdiscreteBC,&
!                                        getBoundaryValues_constC,rcollection)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
!       call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                        rboundaryRegion,rdiscreteBC,&
!                                        getBoundaryValues_constC,rcollection)
    call collct_done(rcollection)
!     end if

    end subroutine

!##################################################################################
!##################################################################################

! Constructing the RHS of the chemoattractant part
    subroutine chemo_initrhsC (rrhsBlockchemo , rchemoattract ,rhschemo , rcell ,&
                                                 rdiscretisation , rmassmatrix , dtstep,  PHI)

    type(t_vectorScalar) , intent(IN) :: rchemoattract , rcell

    type(t_matrixScalar) , intent(IN) :: rmassmatrix

    type(t_vectorBlock) , intent(INOUT) :: rrhsBlockchemo

    type(t_vectorScalar) , intent(INOUT) :: rhschemo

    type(t_Blockdiscretisation) , intent(IN) :: rdiscretisation

    ! param
    real(DP) ,intent (IN) :: dtstep, PHI

    ! local linform
    type(t_linearForm) :: rlinform

    ! local collection
    type(t_collection) :: rcollection


    call collct_init (rcollection)
    call collct_setvalue_vecsca (rcollection, 'cbvector', rcell, .true.)
    rcollection%DquickAccess(1) = PHI
    ! To assemble the RHS , set up the corresponding linear  form (u*g(u),Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1), &
                                rlinform, .true., rhschemo, coeff_hillenX_RHS_c, rcollection)

    ! Since we approx. the 1st deriv. by the backward euler, there's still one component
    !missing in the RHS:
    ! rmassmatrix* c_{old}
    ! The corressponding massmatrix is already built ( just before the loop starts )
    call lsyssc_scalarMatVec(rmassmatrix,rchemoattract,rhschemo,1.0_DP,dtstep)

      ! Release the collection structure
      call collct_done(rcollection)
      ! Now we should have the correct RHS = rhs like its mentioned above.

    end subroutine

!##################################################################################
!##################################################################################


! This subroutine is called to implement a basic defect correction to deal with possible
! nonlinearities of CHI
    subroutine chemo_defcorr(  rcell , rchemoattract , rcellBlock , rmassmatrix , rlaplace ,&
               rrhscell , rdiscretisation , rdiscreteBC , dtstep , D_1, CHI, ALPHA, r, GAMMA, N, maxiterationdef , defectTol,&
                iteration_u_max, iteration_u_min, iteration_u_average, iteration_defcorr_max, &
                iteration_defcorr_min, iteration_defcorr_average)
    

    ! sol vectors
    type(t_vectorScalar) ,  intent(INOUT) :: rcell , rchemoattract
    ! The solution block -vector
    type(t_vectorBlock) , intent(INOUT) :: rcellBlock
    ! Some matrices for computing the LHS
    type(t_matrixScalar) , intent(INOUT) :: rmassmatrix , rlaplace
    ! The RHs of the init pb
    type(t_vectorScalar)  , intent(IN) :: rrhscell
    ! The underlying discretisation and BCs
    type(t_Blockdiscretisation) , intent(IN) :: rdiscretisation
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
    ! temporal matrix
    type(t_matrixScalar) :: rK , rmatrix, rtempLinComb, rlumpedmass
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
    real(DP), dimension(:), allocatable ::dedge, aedge_mass
    integer, dimension(:,:), allocatable ::kedge
    ! a local average iterationcounter for  iteration_u_average
    real(DP) :: u_average_local


    u_average_local = 0.0_DP

    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rtempLinComb)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rK)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9,rmatrix)
    
    def_correction : do k=1,maxiterationdef

            ! STEP 2.2: Assemble the system matrix M_total
            ! Now we are constructing the Systemmatrix for the implicit scheme. @# ( [M + dt*L -dt*M_tri] u_{n+1} ) #@

            ! copy the already constructed matrix to rmatrix
            call lsyssc_copyMatrix(rmassmatrix, rmatrix)
            ! The only matrix-part which have to be computed within the timeloop
            ! We use a bilf to obtain a suitable matrix
            ! Since the CHI is free to be nonlinear it should be implemented directly
            ! in the cb fct. we ' ll invoke in assembling the matrix
 
            call  lsyssc_lumpMatrixScalar (rmatrix,LSYSSC_LUMP_DIAG,.false.)
            call lsyssc_copyMatrix ( rmatrix, rlumpedmass )

            ! Initialize the collection structure
            call collct_init (rcollection)

            ! Setting up the collection structure
            ! We're considering the case of constant chemosensitivityfct CHI
            call collct_setvalue_vecsca ( rcollection , 'cbvector1' , rchemoattract , .true.)
            call collct_setvalue_vecsca ( rcollection , 'cbvector2' , rcell , .true.)
            rcollection%DquickAccess(1) = dtstep
            rcollection%DquickAccess(2) = CHI
            rcollection%DquickAccess(3) = GAMMA
            rcollection%DquickAccess(4) = ALPHA
            rform%itermCount = 2
            rform%Idescriptors(1,1) = DER_FUNC
            rform%Idescriptors(2,1) = DER_DERIV_X
            rform%Idescriptors(1,2) = DER_FUNC
            rform%Idescriptors(2,2) = DER_DERIV_Y
            rform%ballCoeffConstant = .false.
            rform%BconstantCoeff(1) = .false.
            rform%BconstantCoeff(2) = .false.
            call bilf_buildMatrixScalar (rform,.true.,rK, coeff_hillenX, rcollection)



            if(k.eq.1) then
                allocate ( kedge ( 2, rK%NA ) )
                allocate ( dedge ( rK%NA ) )
                allocate ( aedge_mass ( rK%NA ) )
            endif

            !**********************
            ! Adding some artificial diffusion to obtain positivity + smooth sol.
            ! K->K* (rK)
!             call chemo_artdiff( rmassmatrix, rK, dedge, kedge, nedge, aedge_mass )
            !**********************

! This is commented, since we implemente4d this part as a fixed RHS to avoid big
! computations
! It's implemented approx. at line 617

!             ! Adding the logistic growth that Hillen and Painter introduced in their paper, (M8)
!             rcollection%DquickAccess(1) = dtstep
!             rcollection%DquickAccess(2) = r
!             rform%itermCount = 1
!             rform%Idescriptors(1,1) = DER_FUNC
!             rform%Idescriptors(2,1) = DER_FUNC
!             rform%ballCoeffConstant = .false.
!             rform%BconstantCoeff(1) = .false.
!             call bilf_buildMatrixScalar (rform,.false.,rK, coeff_m8bilf, rcollection)

            ! Adding the logistic growth of u^2(1-u) as proposed in the pattern-paper
            ! Here we implement this term as a LHS NOT a RHS
            ! In this case (not commented) the construction of the RHS_u should be modified
             rcollection%DquickAccess(1) = dtstep
            rform%itermCount = 1
            rform%Idescriptors(1,1) = DER_FUNC
            rform%Idescriptors(2,1) = DER_FUNC
            rform%ballCoeffConstant = .false.
            rform%BconstantCoeff(1) = .false.
            call bilf_buildMatrixScalar (rform,.false.,rK, coeff_pattern_growthterm, rcollection)

            !call chemo_artdiff( rmassmatrix, rK, dedge, kedge, nedge, aedge_mass )

            ! Since we use a 1st order approx of the time-deriv. (backward euler) we add the
            ! lumped massmatrix to the existing matrix
            call lsyssc_matrixLinearComb(rmatrix, 1.0_DP,rK , -dtstep, rtempLinComb, .true. , .true. , .true.)
            call lsyssc_copyMatrix(rtempLinComb, rmatrix)


            ! Adding the diffusion-part of the model, since this can also be non-linear
            ! we use the callback function coeff_hillen_laplace
            rcollection%DquickAccess(2) = D_1
            rcollection%DquickAccess(3) = N
            rform%itermCount = 2
            rform%Idescriptors(1,1) = DER_DERIV_X
            rform%Idescriptors(2,1) = DER_DERIV_X
            rform%Idescriptors(1,2) = DER_DERIV_Y
            rform%Idescriptors(2,2) = DER_DERIV_Y
            rform%ballCoeffConstant = .false.
            rform%BconstantCoeff(1) = .false.
            rform%BconstantCoeff(2) = .false.
            call bilf_buildMatrixScalar (rform,.false.,rmatrix, coeff_hillen_laplace, rcollection)


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
        
            
            ! Else we continue with deriving the defect correction solution...

            ! Hang the pointer into the vector and matrix. That way, these
            ! boundary conditions are always connected to that matrix and that
            ! vector.
            rmatrixBlock%p_rdiscreteBC => rdiscreteBC
            rrhsBlock%p_rdiscreteBC => rdiscreteBC
            
            ! Now we have block vectors for the RHS and the matrix. What we
            ! need additionally is a block vector for the solution and
            ! temporary data. Create them using the RHS as template.
            ! Fill the solution vector with 0:
            ! CALL lsysbl_createVecBlockIndirect (rrhsBlock, rcellBlock, .TRUE.)
            ! call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
            
            ! Next step is to implement boundary conditions into the RHS,
            ! solution and matrix. This is done using a vector/matrix filter
            ! for discrete boundary conditions.
            ! The discrete boundary conditions are already attached to the
            ! vectors/matrix. Call the appropriate vector/matrix filter that
            ! modifies the vectors/matrix according to the boundary conditions.
            call vecfil_discreteBCrhs (rrhsBlock)
            call vecfil_discreteBCsol (rcellBlock)
            call matfil_discreteBC (rmatrixBlock)

            ! Since we use defect correction for the nonlinearness, we should compute the defect
            ! e.g.  rdefRHS (= rrhs ) = rrhs - rmatrix * rvector
            call lsysbl_copyVector ( rrhsBlock , rdefBlock )
            call lsysbl_blockMatVec ( rmatrixBlock , rcellBlock , rdefBlock , - 1.0_DP , 1.0_DP )

            call vecfil_discreteBCdef ( rdefBlock )
            
            ! Calculate the l2_norm of the defect
            defect = lsysbl_vectorNorm ( rdefBlock , LINALG_NORML2 )

            !print *, " residuum : " , defect

            ! Here is the time to CHECK SOME RESIDUAL CONVERGENCE to exit the iteration of the
            ! nonlinear pb
            ! But since we 're testing just constant CHI there's nether the need to iterate nor to check
            ! residual issues.
            ! Certain defect error control should be implemented here
            if ( k .gt. 1 ) then
                if ( defect .le. defectTol) then
                        ! Cells Releasements
                        call lsysbl_releaseVector(rdefBlock)
                        call lsysbl_releaseVector (rrhsBlock)
                        call lsysbl_releaseMatrix (rmatrixBlock)
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
            ! So, set up a filter chain that filters the defect vector
            ! during the solution process to implement discrete boundary conditions.
            RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
            ! Create a BiCGStab-solver. Attach the above filter chain
            ! to the solver, so that the solver automatically filters
            ! the vector during the solution process.
            p_RfilterChain => RfilterChain
            nullify(p_rpreconditioner_cells)
            call linsol_initBiCGStab (p_rsolverNode_cells,p_rpreconditioner_cells,p_RfilterChain)
        
            ! Attach the system matrix to the solver.
            ! First create an array with the matrix data (on all levels, but we
            ! only have one level here), then call the initialisation
            ! routine to attach all these matrices.
            ! Remark: Don't make a call like
            !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
            ! This doesn't work on all compilers, since the compiler would have
            ! to create a temp array on the stack - which does not always work!
            Rmatrices = (/rmatrixBlock/)
            
            call linsol_setMatrices(p_RsolverNode_cells,Rmatrices)
            
            ! Initialise structure/data of the solver. This allows the
            ! solver to allocate memory / perform some precalculation
            ! to the problem.
            call linsol_initStructure (p_rsolverNode_cells, ierror)
            if (ierror .ne. LINSOL_ERR_NOERROR) stop
            call linsol_initData (p_rsolverNode_cells, ierror)
            if (ierror .ne. LINSOL_ERR_NOERROR) stop
            ! Finally solve the system. As we want to solve Ax=b with
            ! b being the real RHS and x being the real solution vector,
            ! we use linsol_solveAdaptively. If b is a defect
            ! RHS and x a defect update to be added to a solution vector,
            ! we would have to use linsol_precondDefect instead.
!             call linsol_solveAdaptively (p_rsolverNode_cells,rcellBlock,rrhsBlock,rtempBlock)

            ! Set the output level of the solver to 2 for some output
            p_rsolverNode_cells%ioutputLevel = 0
            p_rsolverNode_cells%depsRel=1E-16_DP
            p_rsolverNode_cells%nminIterations=1
            p_rsolverNode_cells%nmaxIterations=500
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
            
            ! test block for the linearity
            if(k.eq.3) then
                print *,' problem is nonlinear!!! '
!                 stop
            end if
            
    end do def_correction

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
    !call chemo_fct_lim ( rdiscretisation, rcell, kedge, dedge, nedge, aedge_mass,&
    !                     rlumpedmass, rlaplace, rK, dtstep )

    deallocate ( kedge )
    deallocate ( dedge )
    deallocate ( aedge_mass )

    call lsyssc_releaseMatrix (rK)
    call lsyssc_releaseMatrix (rlumpedmass)
    call lsyssc_releaseMatrix (rtempLinComb)
    call lsyssc_releaseMatrix (rmatrix)
    
    end subroutine


!##################################################################################
!##################################################################################

    ! This subroutine adds some artificial diffusion to the matrix rmatrix
    subroutine chemo_artdiff( rmassmatrix, rmatrix, dedge, kedge, nedge, aedge_mass )

    ! The matrix wich is to be modified
    type(t_matrixScalar) , intent(INOUT) :: rmassmatrix, rmatrix

    real(DP), dimension (:), intent(OUT) :: dedge
    integer, dimension(:,:), intent(OUT) :: kedge
    integer, intent(INOUT) :: nedge
    real(DP), dimension (:), intent(OUT) :: aedge_mass

    ! some local variables
    integer :: i , j , k, nvt ,counter, iedge
    integer :: ij_loc , ii_loc , ji_loc , jj_loc
    real(DP) :: a,b, d_ij
    integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
    integer(PREC_MATIDX), dimension(:), pointer :: p_Kld , p_Kdiagonal
    real(DP), dimension(:), pointer :: p_Da_mass, p_Da
    integer , dimension(:) , allocatable :: iaux

    ! Get the structure (rmassmatrix and rmatrix have the same structure)
    call lsyssc_getbase_double (rmassmatrix,p_Da_mass)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
    call lsyssc_getbase_double (rmatrix,p_Da)
    nvt = rmatrix%NEQ

    allocate ( iaux( nvt ) )

    ! Correction of the discrete convection operator
    iedge = 0
    iaux ( : ) = p_Kld( : )
    DO i = 1,nvt
        ii_loc = p_Kdiagonal ( i )

        DO ij_loc = iaux ( i )+1 ,p_Kld( i+1 ) - 1
        j = p_Kcol ( ij_loc )
        jj_loc = p_Kdiagonal ( j )
        ji_loc = iaux ( j )
        iaux ( j )=iaux ( j )+1
        
        ! Artificial diffusion coefficient
        d_ij = MAX( -p_Da ( ij_loc ), 0.0_DP, -p_Da ( ji_loc ) )

        ! Elimination of negative entries
        p_Da ( ij_loc ) = p_Da ( ij_loc ) + d_ij
        p_Da ( ii_loc ) = p_Da ( ii_loc ) - d_ij
        p_Da ( ji_loc ) = p_Da ( ji_loc ) + d_ij
        p_Da ( jj_loc ) = p_Da ( jj_loc ) - d_ij
        
        iedge = iedge + 1
        dedge ( iedge ) = d_ij
        kedge ( 1, iedge ) = i
        kedge ( 2, iedge ) = j
        aedge_mass( iedge ) =  p_Da_mass( ij_loc )
        
        END DO
    END DO

    nedge = iedge
    
    deallocate ( iaux )

    end subroutine

!##################################################################################
!##################################################################################

    ! Here we set up the antidiffusive fluxes as well as the diff.limiters
    ! This routine should follow the ideas of D.Kuzmins and M.Moellers "Algebraic Flux Correction I. Scalar Conservation Laws", March 2004
    subroutine chemo_fct_lim ( rdiscretisation, rvector, kedge, dedge, nedge, aedge_mass,&
                               rlumpedmatrix, rlaplace, rK, dtstep )
                         
        type(t_vectorScalar), intent(INOUT) :: rvector
        type(t_matrixScalar), intent(IN) :: rlumpedmatrix, rlaplace, rK
        type(t_Blockdiscretisation), intent(IN) :: rdiscretisation
        real(DP), dimension(:), intent(IN) :: dedge, aedge_mass
        integer, dimension(:,:), intent(IN) :: kedge
        integer, intent(IN) :: nedge
        real(DP) :: dtstep

        ! Some local variables
        integer :: i, j, iedge, nedge_mass, nedge_D, neq
        real(DP) :: f_ij, du, eps, counter
        real(DP), dimension (: ), allocatable  :: pp, pm, qp, qm, rm, rp, f, flux
        type(t_vectorScalar) :: ru_dot, rvector_temp
        type(t_matrixScalar) :: rD, rtemp
        real(DP), dimension(:), pointer :: p_Da, p_vectorentries, p_udot
        integer , dimension(:), pointer :: p_Kdiagonal
        neq = rK%NEQ
        eps = 1e-16

        allocate ( pp (neq) )
        allocate ( pm (neq) )
        allocate ( qp (neq) )
        allocate ( qm (neq) )
        allocate ( rp (neq) )
        allocate ( rm (neq) )
        allocate ( f ( neq ) )

        ! Initialize the auxiliary arrays
        call chemo_initAuxArrays ( pp, pm, qp, qm, rp, rm, f )

        call lsyssc_getbase_double( rvector , p_vectorentries )
        call lsyssc_getbase_Kdiagonal ( rlumpedmatrix , p_Kdiagonal )
        call lsyssc_getbase_double(rlumpedmatrix, p_Da)

        ! Building of the raw antidiffusive fluxes f_ij
        allocate ( flux ( nedge ) )

        !Constructing apprx. derivative
        !initialization
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),ru_dot,.true.)
        call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rvector_temp,.true.)
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,rtemp)
        !construction
        call lsyssc_matrixLinearComb(rK, 1.0_DP, rlaplace, -1.0_DP, rtemp, .true. , .true. , .true.)
        call lsyssc_scalarMatVec(rtemp, rvector, rvector_temp, 1.0_DP, 0.0_DP)
        call lsyssc_invertedDiagMatVec ( rlumpedmatrix, rvector_temp, 1.0_DP, ru_dot)
        !get ru_dot data
        call lsyssc_getbase_double( ru_dot, p_udot )
  
        !Constructing fluxes
        DO iedge = 1,nedge
            flux ( iedge ) = aedge_mass (iedge) * ( p_udot ( kedge (1,iedge) ) - &
                                     p_udot ( kedge (2,iedge) ) ) + &
                             dedge ( iedge ) * ( p_vectorentries ( kedge (1,iedge) ) - &
                                     p_vectorentries ( kedge (2,iedge) ) )
            !flux ( iedge ) = dedge ( iedge ) * ( p_vectorentries ( kedge (1,iedge) ) - &
            !                         p_vectorentries ( kedge (2,iedge) ) )
        END DO


        ! Now the actual implementation of the FCT stabilization for the Keller-Segel model starts
        DO iedge = 1 , nedge

            i = kedge ( 1, iedge )
            j = kedge ( 2, iedge )

            ! Antidiffusive fluxes to be limeted
            du = p_vectorentries ( j ) - p_vectorentries ( i )
            f_ij = flux ( iedge )

            ! Prelimiting of antidiffusive fluxes
            if ( f_ij * du > 0 ) then
                f_ij = 0 ; flux ( iedge ) = f_ij
            end if

            ! Positive/negative edge contributions
            pp ( i ) = pp ( i ) + MAX ( 0.0_DP , f_ij )
            pp ( j ) = pp ( j ) + MAX ( 0.0_DP , -f_ij )
            pm ( i ) = pm ( i ) + MIN ( 0.0_DP , f_ij )
            pm ( j ) = pm ( j ) + MIN ( 0.0_DP , -f_ij )

            ! Maximum / minimum solution increments
            qp ( i ) = MAX ( qp ( i ) , du )
            qp ( j ) = MAX ( qp ( j ) , -du )
            qm ( i ) = MIN ( qm ( i ) , du )
            qm ( j ) = MIN ( qm ( j ) , -du )

        END DO

        ! Computation of nodal correction factors
        WHERE ( pp > eps )  rp = MIN ( 1.0_DP, p_Da ( p_Kdiagonal ( : ) ) * qp / (dtstep*pp ) )
        WHERE ( pm < -eps ) rm = MIN ( 1.0_DP, p_Da ( p_Kdiagonal ( : ) ) * qm / (dtstep*pm ) )
             
        ! Correction of the right-hand side
        DO iedge = 1 , nedge

            ! Node numbers for the current edge
            i = kedge ( 1 , iedge )
            j = kedge ( 2 , iedge )

            ! Antidiffusive flux to be limeted
            f_ij = flux ( iedge )

            ! Multiplication by alpha_ij
            ! f_ij = alpha_ij * f_ij
            IF ( f_ij  > 0 ) THEN
                f_ij = MIN ( rp ( i ) , rm ( j ) ) * f_ij
            ELSE
                f_ij = MIN ( rm ( i ) , rp ( j ) ) * f_ij
            END IF

            ! Insertion into the global vector
            f ( i ) = f ( i ) + f_ij
            f ( j ) = f ( j ) - f_ij

        END DO

        DO i = 1, neq
            ! The actual reconstruction to the "mid-order" solution
            ! u_i ^{n+1} = u_i ^L + \Delta t * f_i (u^L , u^n) / m_i
            p_vectorentries ( i ) = p_vectorentries ( i ) + dtstep * ( f ( i ) / p_Da ( p_Kdiagonal ( i ) ) )
        END DO

        ! Cleaning up the memory
        deallocate ( pp )
        deallocate ( pm )
        deallocate ( qp )
        deallocate ( qm )
        deallocate ( rm )
        deallocate ( rp )
        deallocate ( f )
        deallocate ( flux )
        call lsyssc_releaseMatrix (rtemp)
        call lsyssc_releaseVector (ru_dot)
        call lsyssc_releaseVector (rvector_temp)
        
    end subroutine


!##################################################################################
!##################################################################################
    ! This routine sets all entries of the arrays to ZERO
    subroutine chemo_initAuxArrays ( pp, pm, qp, qm, rp, rm, f )

        real(DP), dimension (:), intent(INOUT) :: pp, pm, qp, qm, rp, rm, f

        ! locla variables
        integer :: i, na

        na = size ( pp )

        DO i = 1, na

            pp ( i ) = 0.0_DP
            pm ( i ) = 0.0_DP
            qp ( i ) = 0.0_DP
            qm ( i ) = 0.0_DP
            rp ( i ) = 1.0_DP
            rm ( i ) = 1.0_DP
            f ( i ) = 0.0_DP

        END DO

    end subroutine

!##################################################################################
!##################################################################################

!!##################################################################################
!
!    ! For testing purpose only
!    ! This subroutine should construct the arrays which are used in the CES matrix format
!    ! the arrays should be deallocated,  as soon this routine is called
!    subroutine chemo_makeCES(rdiscretisation, rmatrix, kedge , aedge , nedge)
!
!        ! The matrix which should be considered
!        type(t_matrixScalar), intent(IN) :: rmatrix
!
!        type(t_Blockdiscretisation), intent(IN) :: rdiscretisation
!
!        ! Arrays which will contain the CES structure
!        real(DP), dimension(:,:), intent(INOUT) :: aedge
!        integer, dimension(:,:), intent(INOUT) :: kedge
!
!        integer, intent( OUT ) :: nedge
!
!        ! local variables
!        type(t_matrixScalar) :: rmatrixT
!        integer :: i, j, ji_loc, ij_loc, nvt, na, bias, iedge
!        integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol, p_KcolT
!        integer(PREC_MATIDX), dimension(:), pointer :: p_Kld, p_Kdiagonal, p_KldT, p_KdiagonalT
!        real(DP), dimension(:), pointer :: p_Da, p_DaT
!
!
!        ! Get the transposed matrix
!        call lsyssc_transposeMatrix ( rmatrix, rmatrixT )
!
!        ! Get the current matrix structure
!        call lsyssc_getbase_Kld (rmatrix,p_Kld)
!        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
!        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
!        call lsyssc_getbase_double (rmatrix,p_Da)
!        call lsyssc_getbase_Kld (rmatrixT,p_KldT)
!        call lsyssc_getbase_Kcol (rmatrixT,p_KcolT)
!        call lsyssc_getbase_Kdiagonal(rmatrixT, p_KdiagonalT)
!        call lsyssc_getbase_double (rmatrixT,p_DaT)
!
!
!
!
!        nvt = rmatrix%NEQ
!        na = rmatrix%NA
!
!        ! Calculating the maximal number of nonzero off-diagonal entries , e.g. nedge
!        nedge = (na-nvt) / 2
!
!        ! Allocate some mem for the CES-Arrays
!!         allocate ( kedge( 2,nedge ) )
!!         allocate ( aedge( 2,nedge ) )
!
!        iedge = 1
!
!
!        ! build the arrays
!        DO i = 1, nvt-1
!
!        ! resetting of the bias, since we start a new "row"
!        bias = 0
!
!            ! loop through the entry-array of rmatrix. We're concernd about the
!            ! upper triangular entries which are non-zero
!            DO ij_loc = p_Kdiagonal( i ) +1 , p_Kld ( i+1 ) -1
!
!                ! loop through the entry-array of rmatrixT. we 're concerned about the
!                ! lower triangular entries which are non-zero
!
!                DO ji_loc = p_KdiagonalT ( i ) +1 +bias , p_KldT ( i+1 ) -1
!                    IF ( p_Kcol ( ij_loc ) > p_KcolT ( ji_loc ) ) THEN
!                        ! edge is induced by a lower triangular entry
!                        j = p_Kcol ( ji_loc )
!                        kedge ( 1 , iedge ) = i
!                        kedge ( 2 , iedge ) = j
!                        aedge ( 1 , iedge ) = 0.0_DP
!                        aedge ( 2 , iedge ) = p_DaT ( ji_loc )
!
!                        ! Increment edgeindex
!                        iedge = iedge +1
!                    ELSE
!                    IF ( p_Kcol ( ij_loc ) == p_KcolT ( ji_loc ) ) THEN
!                        ! edge is induced by both lower  and upper triangular entry
!                        j = p_Kcol ( ji_loc )
!                        kedge ( 1 , iedge ) = i
!                        kedge ( 2 , iedge ) = j
!                        aedge ( 1 , iedge ) = p_Da ( ij_loc )
!                        aedge ( 2 , iedge ) = p_DaT ( ji_loc )
!
!                        ! Increment edgeindex
!                        iedge = iedge +1
!                        ! skipping the concerning entry for the next comparison
!                        bias = ji_loc - p_KdiagonalT ( i )
!                        EXIT
!                    ELSE
!                    IF ( p_Kcol ( ij_loc ) < p_KcolT ( ji_loc ) ) THEN
!                        ! edge is induced by a upper triangular entry
!                        j = p_Kcol ( ij_loc )
!                        kedge ( 1 , iedge ) = i
!                        kedge ( 2 , iedge ) = j
!                        aedge ( 1 , iedge ) = p_Da ( ij_loc )
!                        aedge ( 2 , iedge ) = 0.0_DP
!
!                        ! Increment edgeindex
!                        iedge = iedge +1
!                        ! skipping the concerning entries for the next comparison
!                        bias = ji_loc - p_KdiagonalT ( i ) -1
!                        EXIT
!                    END IF
!                    END IF
!                    END IF
!
!                END DO
!
!                IF ( ji_loc >= p_KldT ( i+1 ) ) THEN
!                        ! the upper triangular-part is not finished yet....
!                        j = p_Kcol ( ij_loc )
!                        kedge ( 1 , iedge ) = i
!                        kedge ( 2 , iedge ) = j
!                        aedge ( 1 , iedge ) = p_Da ( ij_loc )
!                        aedge ( 2 , iedge ) = 0.0_DP
!
!                        ! Increment edgeindex
!                        iedge = iedge +1
!                END IF
!
!
!            END DO
!
!                IF ( ij_loc == p_Kld ( i+1 ) .AND.  ji_loc < p_KldT ( i+1 )) THEN
!                    ! the lower triangular-part is not finished yet....
!                    DO ji_loc = p_KdiagonalT ( i ) +1 +bias , p_KldT ( i+1 ) -1
!
!                        j = p_Kcol ( ji_loc )
!                        kedge ( 1 , iedge ) = i
!                        kedge ( 2 , iedge ) = j
!                        aedge ( 1 , iedge ) = 0.0_DP
!                        aedge ( 2 , iedge ) = p_DaT ( ji_loc )
!
!                        ! Increment edgeindex
!                        iedge = iedge +1
!                    END DO
!                END IF
!
!
!        END DO
!
!        ! Actualize the number of edges
!        nedge = iedge -1
!
!
!        call lsyssc_releaseMatrix(rmatrixT)
!
!
!    end subroutine
!
!
!!##################################################################################

end module

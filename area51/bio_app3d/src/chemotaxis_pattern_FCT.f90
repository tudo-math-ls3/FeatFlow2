
!##############################################################################
!# ****************************************************************************
!# <name> chemotaxis_pattern (FCT) 3D </name>
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
!#   [M + dt*L + dt*M] c_{n+1} 		= [ M ] c_{n} + dt* M u_n
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Definitions of variables. !!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! A parameter list structure accepting the parameters from the DAT file.
    type(t_parlist) :: rparams
    !
    !!!!! parameters' variables !!!!! 
    ! maximum level of refinment
    integer :: NLMAX    
    !
    ! Defining the output level
    ! =0 means no gmv files generated
    ! =1 means gmv files printed
    ! =2 means just the gmv file is printed
    ! =3 means just the ISTEP_GMV-th gmv file are printed
    integer :: OUTPUT, ISTEP_GMV
    !
    ! Whether or not taking an const initial sol
    ! =0 means const initial sol
    ! =1 means exponential initial sol
    ! Whether or not we set the IC by l_2 projection ( = 1 ) or by interpolation ( =0 )
    integer :: INITSOL , L2PROJ
    !
    ! Defining the gmvoutputfolder
    integer :: gmvfolder, checkneg
    !
    ! Variables for the defect correction ( nonlinear loop )
    real(DP) :: defect , defectTol , negthres
    !
    ! Some params used to adjust the model
    real(DP) :: CHI, D_1, D_2, A_CHEMO , A_CELLS , B_CHEMO , B_CELLS, W, SIGMA, BETA, ALPHA, R,  PHI, GAMMA, N
    !
    ! If INITSOL = 0 , meaning const init sol, we use these variables to
    ! define the const value
    real(DP) :: C_0, U_0
    !
    ! Scaling params for the analytic given sols
    real(DP) :: SCALE_U , SCALE_C
    !
    ! Time step size, number of timesteps.(also for a fixed amount of steps, since there is
    ! a basic loop control implemented dealing w/ the steady state analysis)
    real(DP) :: dtstep
    integer :: ntimesteps, steps
    !
    ! Time and time step counter
    real(DP) :: dtime
    integer :: itimestep
    !
    ! maximal iterations for the defect correction . Since we allow possible nonlinearities
    integer :: maxiterationdef
    !
    ! error-norm constant
    ! To be set in the .dat file
    integer :: CTRLNORM
    !
    ! error analysis
    real(DP) :: uerror, cerror, tol

    
    !!!!! triangulation variables !!!!! 
    ! maximum level of refinment
    ! 
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation
    !
    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir
    !
    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation


    !!!!! solution + rhs + def vectors !!!!!
    type(t_vectorScalar) :: ranalyticcells , ranalyticchemo, rcold, ruold
    !A pointer to the entries of vector rchemoattract				
    real(DP), dimension(:), pointer ::  p_analyticcells, p_analyticchemo, p_cold, p_uold    
    real(DP), dimension(:), pointer ::  p_vectordata, p_chemodata
    !
    ! cell and chemoattractant solution vectors
    type(t_vectorScalar)  :: rchemoattract, rrhschemo, rcell, rrhscell
    !
    ! defect vectors
    type(t_vectorScalar)  :: rdef
    !
        
    
    !!!!! matrices !!!!!
    type(t_matrixScalar)  :: rmassmatrix, rsysmatrix, rlaplace, rmatrixchemo, rtemp
    
    
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Ok, let's start. !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!  Reading the mesh file + triangulation !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get the path $PREDIR from the environment, where to read .prm/.tri files 
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './tri'

    ! At first, read in the basic triangulation.
    call tria_readTriFile3D (rtriangulation, trim(spredir)//'/cCube3.tri')
   
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation)

    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_Q1_3D,CUB_G3_3D,rtriangulation)
    
    !!!!! Missed the boundary conditions ??? !!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!  Creating matrices !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    call chemo_creatematvec ( rmassmatrix , rsysmatrix, rlaplace , rmatrixchemo ,&
                                             rchemoattract ,rrhschemo , rcell , rrhscell , rdef , rdiscretisation )    
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!  Working with initial and boundary conditions !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get some pointers  for the errorctrl
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),ruold,.true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rcold,.true.)
    call lsyssc_getbase_double(rcold,p_cold) 
    call lsyssc_getbase_double(ruold,p_uold)
    !
    
    




    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!  Releasing data  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! releasing old vectors
    call lsyssc_releaseVector (rcold)
    call lsyssc_releaseVector (ruold)

    ! releasing matrices and vectors
     call chemo_releasematvec ( rmassmatrix , rsysmatrix, rlaplace , rmatrixchemo ,&
                                             rchemoattract ,rrhschemo , rcell , rrhscell , rdef)    

    ! Release parameterlist
     call parlst_done (rparams)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    call tria_done (rtriangulation)
    
    
    ! Definitions of variables.
    print *,' 3D version responds ... '
    print *,''
    
  end subroutine
  !!! end of chemotaxispatternFCT !!!




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!! Auxiliary subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    ! In this subroutine the basic matrices and vectors are initialized
    subroutine chemo_creatematvec ( rmassmatrix , rsysmatrix, rlaplace , rmatrixchemo ,&
                                                        rchemoattract , rrhschemo , rcell , rrhscell , rdef , rdiscretisation )

    ! This is where the matrices shoiuld be stored
    type(t_matrixScalar) , intent(INOUT)  ::rmassmatrix, rsysmatrix, rlaplace , rmatrixchemo

    ! This are the vectors
    type(t_vectorScalar) , intent(INOUT) :: rchemoattract, rcell, rdef ,rrhschemo , rrhscell

    ! The underlying discretisation
    type(t_blockDiscretisation) , intent (IN) :: rdiscretisation


    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rmassmatrix)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rsysmatrix)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rlaplace)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9,rmatrixchemo)
    
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rchemoattract,.true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rrhschemo, .true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rcell,.true.)
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rrhscell,.true.) 
    call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1),rdef,.true.)

    end subroutine








    ! In this subroutine the basic matrices and vectors are initialized
    subroutine chemo_releasematvec ( rmassmatrix , rsysmatrix, rlaplace , rmatrixchemo ,&
                                                        rchemoattract , rrhschemo , rcell , rrhscell , rdef)

    ! This is where the matrices shoiuld be stored
    type(t_matrixScalar) , intent(INOUT)  ::rmassmatrix, rsysmatrix, rlaplace , rmatrixchemo

    ! This are the vectors
    type(t_vectorScalar) , intent(INOUT) :: rchemoattract, rcell, rdef ,rrhschemo , rrhscell


    ! Release the preallocated matrix and the solution vector.
    call lsyssc_releaseMatrix (rmassmatrix)
    call lsyssc_releaseMatrix (rsysmatrix)
    call lsyssc_releaseMatrix (rlaplace)
    call lsyssc_releaseMatrix (rmatrixchemo)

    call lsyssc_releaseVector (rchemoattract)
    call lsyssc_releaseVector (rrhschemo)
    call lsyssc_releaseVector (rcell)
    call lsyssc_releaseVector (rrhscell)
    call lsyssc_releaseVector (rdef)

    end subroutine


end module
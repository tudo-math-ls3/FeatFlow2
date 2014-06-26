!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method1_deflgmres </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!#
!# This module is based on poisson2d_method0_simple, but using a multi-grid
!# solver.
!# </purpose>
!##############################################################################

module anisotropicdiffusion_methodMG2

  use fsystem
  use genoutput
  use storage
  use linearsolver2
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
  use bcassemblybase
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use ucd
  use pprocerror
  use collection
!  use convection
  use spdiscprojection
  use statistics 
  use sortstrategy
  use sortstrategybase
  use paramlist
  use jumpstabilisation  
  
  use matrixio
  use anisotropicdiffusion_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC

!****! Declaration for using sorting algorithms
     ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy
  
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine anisotropicdiffusionMG2
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
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A variable that is used to specify a region on the boundary.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother,p_rprecond
    type(t_linsolNode), pointer :: p_rspecDefl

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! One level of Deflated GMRES
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo

 ! Collection structure for setting up the RHS
    type(t_collection) :: rcollection

!****! Declare a timer structure
    type(t_timer) :: rtimer

    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
!     ! Maximum Eigenvalue of preconditioned matrix
!     real(DP) :: dmaxeigvalue
!     ! Number of GMRES Iterations at each level
!     integer :: iiterations
!     ! Relaxation Parameter
!     real(DP) :: drelaxation


    ! Maximum number of solver iterations                          
    integer :: imaxiter
    ! Relative stopping criterea
    real(DP) :: dRelError
    
    ! Error indicator during initialisation of the solver
    integer :: ierror

    ! Anisotropic diffusion coefficients
    real(DP) :: Dalpha,Dbeta
    ! Difussion coefficient vector
    real(DP) , dimension(2) :: Dvec
    ! temp vect component
    real (DP):: v1, v2

    ! Final diffusion matrix after rotation
    real(DP), dimension(2,2) :: DdiffusionMatrix

    ! Error of FE function to reference function
    real(DP) :: derror,dtemp,dtemp2
    
    type(t_matrixScalar) :: rtempMatrix
    
!    type(t_convStreamDiff2) :: rsdconfig
!    type(t_vectorBlock) :: rsdconvection
!    
!    type(t_blockDiscretisation) :: rvelDiscr
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir, smaster
    real(DP), dimension(:), pointer :: p_Ddata
!    integer, dimension(:), pointer :: p_Kld

! PRM/TRI file
    character(LEN=SYS_STRLEN) :: sstring,sfilePRM,sfileTRI, sfileGMV,sfileLOG

!   Structure for saving parameters from the DAT file
    type (t_parlist) :: rparams
    
    ! Element type of the discretisation
    integer :: ieltype
    
    ! Whether to convert to triangle mesh
    integer :: iconvertToTriangleMesh
    
     ! Mesh distortion
    real(DP) :: dmeshDistortion
    
    ! Stabilisation
    real(DP) :: dgamma
    integer :: istabilisation
    
    ! Type of solution
    integer :: isolution
    
    ! Type of smoother
     integer :: ismoothertype
    ! Type of MG Cycle
    integer :: icycletype

    ! Sorting strategy
     integer :: isorttype

    ! Some temporary variables
    integer :: i,j,k,l
    integer :: h_Ipermutation
    


    ! Ok, let us start.
    !
    ! +------------------------------------------------------------------------
    ! | READ DAT FILE PARAMETERS
    ! +------------------------------------------------------------------------
    !
    ! Initialise the parameter structure and read the DAT file.
    call parlst_init(rparams)
    
    ! Get the data file.
    call sys_getcommandLineArg(1,smaster,sdefault='./dat/mg_anisotropicdiffusion.dat')
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

    ! GMV Directory
    call parlst_getvalue_string (rparams, '', &
                                 'sucddir', sstring)
    read(sstring,*) sucddir
    ! GMV file
    call parlst_getvalue_string (rparams, '', &
                                 'sfileGMV', sstring)
    read(sstring,*) sfileGMV

    ! LOG file
    call parlst_getvalue_string (rparams, '', &
                                 'sfileLOG', sstring)
    read(sstring,*) sfileLOG

    ! Whether to convert to a triangle mesh
    call parlst_getvalue_int (rparams, '', &
      'iconvertToTriangleMesh', iconvertToTriangleMesh, 0)		! calls the subroutine parlst_getvalue_int_direct

    ! Mesh distortion
    call parlst_getvalue_double (rparams, '', &
      'dmeshDistortion', dmeshDistortion, 0.0_DP)

    ! Element type
    call parlst_getvalue_int (rparams, '', &
                              'ielementType', ieltype, 7)
                              
    ! Type of stabilisation
    call parlst_getvalue_int (rparams, '', &
      'istabilisation', istabilisation, 0)
!    call parlst_getvalue_double (rparams, '', &
!      'dgamma', dgamma, 0.01_DP)

    ! Minimum level
    call parlst_getvalue_int (rparams, '', &
                                     'NLMIN', NLMIN, 2)
    
    ! Maximum level
    call parlst_getvalue_int (rparams, '', &
                                     'NLMAX', NLMAX, 7)

    ! Anisotropic diffusion coefficients
    call parlst_getvalue_double (rparams, '', 'Dalpha', Dalpha, 1.0_DP)
    
    ! Anisotropic diffusion coefficients
    call parlst_getvalue_double (rparams, '', 'Dbeta', Dbeta, 1.0_DP)
    
    ! Diffusion vector
    call parlst_getvalue_double (rparams, '', 'v1', v1, 0.0_DP)
    call parlst_getvalue_double (rparams, '', 'v2', v2, 0.0_DP)
        
    ! solution
    call parlst_getvalue_int (rparams, '', &
                              'isolution', isolution, 0)

    ! Smoother type
    call parlst_getvalue_int (rparams, '', &
                              'ismoothertype', ismoothertype, 0)

    ! Sorting Strategy
    call parlst_getvalue_int (rparams, '', &
                              'isorttype', isorttype, 0)   

    ! Type of Multigrid cycle used. 1 := V cycle , 2:= W cycle. Default is 0 := F cycle
    call parlst_getvalue_int (rparams, '', &
                              'icycletype', icycletype, 0)    
!     ! Max Eigenvalue of preconditioned matrix
!     call parlst_getvalue_double (rparams, '', 'dmaxeigvalue', dmaxeigvalue, 0.0_DP)
!     
!     ! Number of GMRES Iterations at each level
!     call parlst_getvalue_int (rparams, '', 'iiterations', iiterations, -1)
!     
!     ! Relaxation Parameter
!     call parlst_getvalue_double (rparams, '', 'drelaxation', drelaxation, 0.0_DP)

    ! Maximum number of solver iterations                      
    call parlst_getvalue_int (rparams, '', 'imaxiter', imaxiter, 50)

    ! Relative stopping criterea
    call parlst_getvalue_double (rparams, '', 'dRelError', dRelError, 1E-5_DP)
    
    ! Difussion coefficient vector
    Dvec(1) = v1/sqrt(v1*v1+v2*v2)
    Dvec(2) = v2/sqrt(v1*v1+v2*v2)
    
    ! Create the actual diffusion matrix:
    !
    ! A = (  ((alpha - beta)*v1*v1+ beta) 	((alpha - beta)*v1*v2)	   )
    !     (  ((alpha - beta)*v1*v2)         	((alpha - beta)*v2*v2+ beta)  ) 
    !
    DdiffusionMatrix(1,1) = &
      ((Dalpha - Dbeta)*Dvec(1)*Dvec(1)+ Dbeta)

    DdiffusionMatrix(1,2) = &
      ((Dalpha - Dbeta)*Dvec(1)*Dvec(2))  

    DdiffusionMatrix(2,1) = &
       ((Dalpha - Dbeta)*Dvec(1)*Dvec(2))

    DdiffusionMatrix(2,2) = &
        ((Dalpha - Dbeta)*Dvec(2)*Dvec(2)+ Dbeta)        

!****!**********************************************************************
!	Open Data File in the same folder				   *
!***************************************************************************

    OPEN (unit = 7, file = trim(sfileLOG)) ! If file not present, first creates and then opens
    WRITE(7,*) 'Level   ', 'DOF   ','#Iterations   ', 'L2-error   ', 'H1-error'


!        do NLMAX = NLMAX,9
      
        ! Allocate memory for all levels
        allocate(Rlevels(NLMIN:NLMAX))

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Read the domain, read the mesh, refine
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        ! At first, read in the parametrisation of the boundary and save
        ! it to rboundary.
        call boundary_read_prm(rboundary, sfilePRM)
            
        ! Now read in the basic triangulation into our coarse level.
        call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                                sfileTRI, rboundary)
        
        ! Refine it.
        call tria_quickRefine2LevelOrdering (NLMIN-1,&
            Rlevels(NLMIN)%rtriangulation,rboundary)
        
        ! And create information about adjacencies and everything one needs
        ! from a triangulation.
        call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
            rboundary)
        
        ! Now refine the grid for the fine levels.
        do i = NLMIN+1, NLMAX

          ! Refine the grid using the 2-Level-Ordering algorithm
          call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
              Rlevels(i)%rtriangulation,rboundary)
          
          ! Create a standard mesh
          call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
              rboundary)
        
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Set up discretisation structures which tells the code which
      ! finite element to use
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      ! Do this for all levels
      do i = NLMIN, NLMAX
        call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                    Rlevels(i)%rtriangulation, rboundary)
      end do
        
      ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
      ! structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! for this solution component:
      do i = NLMIN, NLMAX
        
      select case (ieltype)
      case (1)
          call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_E001,Rlevels(i)%rtriangulation, rboundary)
      case (11)
          call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_E011,Rlevels(i)%rtriangulation, rboundary)
      case (13)
          call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_E013,Rlevels(i)%rtriangulation, rboundary)
      case (-30)
          call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_EM30,Rlevels(i)%rtriangulation, rboundary)
      case (-1)
          call spdiscr_initDiscr_triquad (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_E001,EL_E011,Rlevels(i)%rtriangulation, rboundary)
      case (-2)
          call spdiscr_initDiscr_triquad (Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          EL_E002,EL_E013,Rlevels(i)%rtriangulation, rboundary)
      end select
   
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Create an assembly information structure on each level which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      do i = NLMIN, NLMAX
        call spdiscr_createDefCubStructure(&  
            Rlevels(i)%rdiscretisation%RspatialDiscr(1),Rlevels(i)%rcubatureInfo,&
            CUB_GEN_AUTO_G3)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Create a 1x1 block matrix with the operator on every level
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create a scalar matrix, based on the discretisation structure
      ! for our one and only solution component.
      do i = NLMIN, NLMAX

      ! Initialise the block matrix with default values based on
      ! the discretisation.
        call lsysbl_createMatBlockByDiscr (&
            Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)
              
              
      select case (istabilisation)
      case (1)
      ! Jump stabilisation. Needs an extended matrix stencil.
          call bilf_createMatrixStructure (&
              Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
              LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1),cconstrType=BILF_MATC_EDGEBASED)
      case DEFAULT
      ! No stabilisation
      call bilf_createMatrixStructure ( &
           Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
           LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))
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
        ! Jump stabilisation. Create the Laplace matrix and add the
        ! stabilisation
        call bilf_buildMatrixScalar (rform,.true.,&
              Rlevels(i)%rmatrix%RmatrixBlock(1,1),Rlevels(i)%rcubatureInfo)
        call jstab_calcUEOJumpStabilisation (&
            Rlevels(i)%rmatrix%RmatrixBlock(1,1),dgamma,dgamma,2.0_DP,1.0_DP,CUB_G3_1D,1.0_DP)
      case DEFAULT
        ! No stabilisation. Create the Laplace matrix directly.
        call bilf_buildMatrixScalar (rform,.true.,&
              Rlevels(i)%rmatrix%RmatrixBlock(1,1),Rlevels(i)%rcubatureInfo)
      end select
               
!     call sstrat_calcCuthillMcKee (Rlevels(i)%rmatrix%RmatrixBlock(1,1), h_Ipermutation)
!          call lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.true.,&
!              SSTRAT_CM,h_Ipermutation)
          
    !      call spdiscr_initBlockDiscr (rvelDiscr, 2, &
    !                                   Rlevels(i)%rtriangulation, rboundary)
    !      call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
    !          rvelDiscr%RspatialDiscr(1), .true.)
    !      call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
    !          rvelDiscr%RspatialDiscr(2), .true.)
    !
    !      call lsysbl_createVectorBlock(rvelDiscr,rsdconvection,.true.)
    !      call lsyssc_clearVector(rsdconvection%RvectorBlock(2),1.0_DP)
    !      rsdconfig%dupsam = 1.0_DP
    !      rsdconfig%dnu = rform%Dcoefficients(1)
    !      rsdconfig%ddelta = 1.0_DP
    !      call conv_streamDiff2Blk2dMat (rsdconfig,Rlevels(i)%rmatrix,rsdconvection,&
    !          rcubatureInfo=Rlevels(i)%rcubatureInfo)
    !      call lsysbl_releaseVector(rsdconvection)
    !      call spdiscr_releaseBlockDiscr(rvelDiscr)
               
        end do
          
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create RHS and solution vectors
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Next step: Create a RHS vector, a solution vector and a temporary
        ! vector. All are filled with zero.
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rrhsBlock,.true.)
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rvectorBlock,.true.)
        call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,rtempBlock,.true.)

!****! Print DoF 
!  	print *,  rvectorBlock%neq
! 	call output_line ('# Dof :', rvectorBlock%neq))
          
        ! Create a collection; used for passing parameters to the RHS.
      	call collct_init (rcollection)
	
		! Put the parameters to the quick-access array
	rcollection%IquickAccess(1) = isolution
	rcollection%DquickAccess(1) = DdiffusionMatrix(1,1)
	rcollection%DquickAccess(2) = DdiffusionMatrix(1,2)
	rcollection%DquickAccess(3) = DdiffusionMatrix(2,1)
	rcollection%DquickAccess(4) = DdiffusionMatrix(2,2)
        !
        ! At first set up the corresponding linear form (f,Phi_j):
        rlinform%itermCount = 1
        rlinform%Idescriptors(1) = DER_FUNC
        
        ! ... and then discretise the RHS to get a discrete version of it.
        ! Again we simply create a scalar vector based on the one and only
        ! discretisation structure.
        ! This scalar vector will later be used as the one and only first
        ! component in a block vector.
        call linf_buildVectorScalar (&
            rlinform,.true.,rrhsBlock%RvectorBlock(1),Rlevels(NLMAX)%rcubatureInfo,coeff_RHS,rcollection)

        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        do i = NLMIN, NLMAX
        
          ! Initialise the discrete BC structure
          call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

          ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.
          call boundary_createRegion(rboundary,1,1,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues)
                                   
          ! Now to the edge 2 of boundary component 1 the domain.
          call boundary_createRegion(rboundary,1,2,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues)
                                   
          ! Edge 3 of boundary component 1.
          call boundary_createRegion(rboundary,1,3,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues)
          
          ! Edge 4 of boundary component 1. That is it.
          call boundary_createRegion(rboundary,1,4,rboundaryRegion)
          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                            rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                            getBoundaryValues)
          
          ! Assign the BC`s to the matrix. That way, these
          ! boundary conditions are always connected to that matrix.
          call lsysbl_assignDiscreteBC(Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
      
          ! Also implement the boundary conditions into the matrix.
          call matfil_discreteBC (Rlevels(i)%rmatrix)



!****! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Impose a resorting into matrix at each level
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	  select case (isorttype)
!****! ROW-WISE Sort Strategy
	    case (1)
!****! Create a sort strategy structure for our discretisation
	      call sstrat_initBlockSorting (Rlevels(i)%rsortStrategy,Rlevels(i)%rdiscretisation)

!****! Calculate the resorting strategy.
	      call sstrat_initXYZsorting (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
				Rlevels(i)%rdiscretisation%RspatialDiscr(1),0)

!****! Impose sorting strategy on matrix at each level
	      call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,Rlevels(i)%rsortStrategy, &
					Rlevels(i)%rsortStrategy)
!****! Resort the matrix at each level
	      call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)

!****! COLUMN-WISE Sort Strategy
	    case(2)
!****! Create a sort strategy structure for our discretisation
	      call sstrat_initBlockSorting (Rlevels(i)%rsortStrategy,Rlevels(i)%rdiscretisation)

!****! Calculate the resorting strategy.
	      call sstrat_initXYZsorting (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
		Rlevels(i)%rdiscretisation%RspatialDiscr(1),1)

!****! Impose sorting strategy on matrix at each level
	      call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,Rlevels(i)%rsortStrategy, &
		Rlevels(i)%rsortStrategy)
!****! Resort the matrix at each level
	      call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
	 
	  end select
          
	end do


        ! Our right-hand-side/solution/temp vectors also needs to 
        ! know the boundary conditions.
        call lsysbl_assignDiscreteBC(rrhsBlock,Rlevels(NLMAX)%rdiscreteBC)
        call lsysbl_assignDiscreteBC(rvectorBlock,Rlevels(NLMAX)%rdiscreteBC)
        call lsysbl_assignDiscreteBC(rtempBlock,Rlevels(NLMAX)%rdiscreteBC)

        ! Next step is to implement boundary conditions into the RHS,
        ! solution and matrix. This is done using a vector/matrix filter
        ! for discrete boundary conditions.
        ! The discrete boundary conditions are already attached to the
        ! vectors. Call the appropriate vector filter that
        ! modifies the vectors according to the boundary conditions.
        call vecfil_discreteBCrhs (rrhsBlock)
        call vecfil_discreteBCsol (rvectorBlock)
        
!        call lsyssc_sortVectorInSitu (rrhsBlock%RvectorBlock(1),&
!            rtempBlock%RvectorBlock(1),-SSTRAT_CM,h_Ipermutation)
!        call lsyssc_sortVectorInSitu (rvectorBlock%RvectorBlock(1),&
!            rtempBlock%RvectorBlock(1),-SSTRAT_CM,h_Ipermutation)
 
!****! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Impose a resorting into vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	select case (isorttype)
	  case (1)
   
!****! Impose the sorting strategy to the vectors.
	    call lsysbl_setSortStrategy (rrhsBlock,Rlevels(NLMAX)%rsortStrategy)
	    call lsysbl_setSortStrategy (rvectorBlock,Rlevels(NLMAX)%rsortStrategy)
    
    
!****! Resort the vectors.
	    call lsysbl_sortVector (rrhsBlock,.true.)
	    call lsysbl_sortVector (rvectorBlock,.true.)

	  case (2)
   
!****! Impose the sorting strategy to the vectors.
	    call lsysbl_setSortStrategy (rrhsBlock,Rlevels(NLMAX)%rsortStrategy)
	    call lsysbl_setSortStrategy (rvectorBlock,Rlevels(NLMAX)%rsortStrategy)
    
    
!****! Resort the vectors.
	    call lsysbl_sortVector (rrhsBlock,.true.)
	    call lsysbl_sortVector (rvectorBlock,.true.)

	end select  
     
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up a linear solver
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! During the linear solver, the boundary conditions are also
        ! frequently imposed to the vectors. But as the linear solver
        ! does not work with the actual solution vectors but with
        ! defect vectors instead.
        ! So, set up a filter chain that filters the defect vector
        ! during the solution process to implement discrete boundary conditions.
        RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

        ! Now we have to build up the level information for multigrid.
        !
        ! Create a Multigrid-solver. Attach the above filter chain
        ! to the solver, so that the solver automatically filters
        ! the vector during the solution process.
        call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)
        
        ! Set up a coarse grid solver.
        ! The coarse grid in multigrid is always grid 1!
        call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
        call linsol_initUMFPACK4 (p_rcoarseGridSolver)
        p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
        
        ! Now set up the other levels...
        do i = NLMIN+1, NLMAX
        
	  select case (ismoothertype)
	    case (0)
          ! Create a Jacobi smoother
	      call linsol_initJacobi(p_rsmoother)
	    case (1)
	    ! Create a Gauss Seidel smoother
	      call linsol_initSOR(p_rsmoother,1.0_DP)
	    case (2)    
	    ! Create an ILU(0) smoother
	      call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
	    case(3)
	      call linsol_initBlockSOR (p_rsmoother, 1.0_DP)
	    end select

          ! We will use 4 smoothing steps with damping parameter 0.7
          call linsol_convertToSmoother(p_rsmoother, 4, 0.7_DP)
          
          ! And add this multi-grid level. We will use the same smoother
          ! for pre- and post-smoothing.
          call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
          p_rlevelInfo%p_rpresmoother => p_rsmoother
          p_rlevelInfo%p_rpostsmoother => p_rsmoother
          
        end do
         
        ! Set the output level of the solver to 2 for some output
        p_rsolvernode%ioutputLevel = 2
        p_rsolvernode%nmaxIterations = imaxiter
        p_rsolvernode%depsRel = dRelError
	p_rsolvernode%p_rsubnodeMultigrid2%icycle = icycletype
	!p_rsolverNode%domega = 0.7_DP
	        
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

        ! We can release Rmatrices immediately -- as long as we do not
        ! release Rlevels(i)%rmatrix!
        do i=NLMIN,NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

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

!****! Initialise the timer by zero:
    call stat_clearTimer(rtimer)
        
!****! Start the timer
	  call stat_startTimer(rtimer)

        ! Finally solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)

!****! Stop the timer
	  call stat_stopTimer(rtimer)

!****! Print the wall clock time that was necessary for the computation  
	WRITE(*,*) "Time for computation: ",rtimer%delapsedReal
!****! Print computation time to output file
	 call output_line ('Time for computation:' // sys_sdEL(rtimer%delapsedReal,10))
        
        print *,"Level=",NLMAX,", ite=",p_rsolverNode%iiterations
        
!        call lsyssc_sortVectorInSitu (rvectorBlock%RvectorBlock(1),rtempBlock%RvectorBlock(1),-SSTRAT_CM)

!****! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Recover the unsorted solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
!****! Unsort the solution
	select case(isorttype) 
	  case(1,2)
	  call lsysbl_sortVector (rvectorBlock,.false.)
        end select
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Postprocessing of the solution
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        ! That is it, rvectorBlock now contains our solution. We can now
        ! start the postprocessing.
        !
        ! Get the path for writing postprocessing files from the environment variable
        ! $UCDDIR. If that does not exist, write to the directory "./gmv".
        if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'
    
        ! Start UCD export to GMV file:
        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
            Rlevels(NLMAX)%rtriangulation,trim(sucddir)//trim(sfileGMV))

	nullify(p_Ddata)
        call spdp_projectToVertices (rvectorBlock%RvectorBlock(1), p_Ddata, DER_FUNC)   
!         call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
        
        ! Write the file to disc, that is it.
        call ucd_write (rexport)
        call ucd_release (rexport)
	WRITE(7,*) ' ',NLMAX, ' ',rvectorBlock%neq, ' ', p_rsolverNode%iiterations
	select case(isolution)
	  case(0,2,3)		
        ! Calculate the error to the reference function.
        call pperr_scalar (PPERR_L2ERROR,derror,rvectorBlock%RvectorBlock(1),&
                           getReferenceFunction,rcollection=rcollection)
        call output_line ('L2-error: ' // sys_sdEL(derror,10) )
	WRITE(7,*) ' ',derror
        call pperr_scalar (PPERR_H1ERROR,derror,rvectorBlock%RvectorBlock(1),&
                           getReferenceFunction,rcollection=rcollection)
        call output_line ('H1-error: ' // sys_sdEL(derror,10) )
        WRITE(7,*) ' ',derror
	end select

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
        
        ! Release the block matrices/vectors
        call lsysbl_releaseVector (rtempBlock)
        call lsysbl_releaseVector (rvectorBlock)
        call lsysbl_releaseVector (rrhsBlock)
        do i = NLMAX, NLMIN, -1
          call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
        end do

        ! Release our discrete version of the boundary conditions
        do i = NLMAX, NLMIN, -1
          call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
        end do

!****! Release the sorting strategy.
	select case (isorttype) 
	  case(1,2)
	    do i = NLMAX, NLMIN, -1
	      call sstrat_doneBlockSorting (Rlevels(i)%rsortStrategy)
	    end do
	end select

        ! Release the cubature info structures
        do i=NLMAX,NLMIN,-1
          call spdiscr_releaseCubStructure(Rlevels(i)%rcubatureInfo)
        end do

        ! Release the discretisation structure and all spatial discretisation
        ! structures in it.
        do i = NLMAX, NLMIN, -1
          call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
        end do
        
        ! Release the triangulation.
        do i = NLMAX, NLMIN, -1
          call tria_done (Rlevels(i)%rtriangulation)
        end do
        
	! Release the collection
     	call collct_done (rcollection)
     
        deallocate(Rlevels)

        ! Finally release the domain, that is it.
        call boundary_release (rboundary)
        
!        end do
      
       print *
!
!****! Close the data file
    CLOSE(7)     


  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> solidmech2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple linear elasticity
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

module solidmech2d_method0_simple

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
  use spdiscprojection
  use scalarpde
  use ucd
  use paramlist
  use pprocerror
  use collection
  use feevaluation
    
  use solidmech2d_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine solidmech2d_0_simple
  
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
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
!     ! An object for saving the triangulation on the domain
!     type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

!     ! An object specifying the discretisation.
!     ! This contains also information about trial/test functions,...
!     type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
!     type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock

    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
!     type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo

    ! A list of points where to evaluate.
    ! DIMENSION(1..ndim,1..npoints)
    real(DP), dimension(2,1)  :: Dpoints

    ! Values of the FE function at the points specified by Dpoints.
    real(DP), dimension(1)    :: Dvalues
    
     ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror1, derror2
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir, smaster,sstring
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
    real(DP), Dimension(:,:), pointer :: p_DvertexCoords

    ! some temporary variables
    integer :: i, j, k, imax
    real(DP) ::  U1, U2, U1_X, U2_X, U1_Y, U2_Y, sigma11, sigma12, sigma22, sigma33
    real(DP) :: devsigma11, devsigma12, devsigma22, devsigma33
    ! Structure for saving parameters from the DAT file
    type (t_parlist) :: rparams

    
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection) :: rcollection
    
    ! Ok, let's start.
    !
    ! +------------------------------------------------------------------------
    ! | READ DAT FILE PARAMETERS
    ! +------------------------------------------------------------------------
    !
    ! Initialise the parameter structure and read the DAT file.
    call parlst_init(rparams)
 
    ! Get the data file.
    call sys_getcommandLineArg(1,smaster,sdefault='./dat/solidmech2d.dat')
    call parlst_readfromfile (rparams, smaster)
    
    ! Get the parameters...
    !
    ! PRM file
    call parlst_getvalue_string (rparams, '', &
                                 'sgridFilePRM', sstring)
    read(sstring,*) rproblem%sgridFilePRM
                                 
    ! TRI file
    call parlst_getvalue_string (rparams, '', &
                                 'sgridFileTRI', sstring)
    read(sstring,*) rproblem%sgridFileTRI
       
    ! Element type
    call parlst_getvalue_string (rparams, '', &
                              'selementType', sstring)
    read(sstring,*) rproblem%selement
                              
    ! number of boundaries (has to be set manually by the user who has to know
    ! the number of boundaries in the current grid)
    call parlst_getvalue_int (rparams, '', &
                                     'nboundaries', rproblem%nboundaries)
                       
    ! number of boundary segments (here, one has to know the number of
    ! boundary segments in the current grid)
    allocate(rproblem%NboundarySegments(rproblem%nboundaries))
    do i = 1,rproblem%nboundaries
      call parlst_getvalue_int (rparams, '', &
                 'NboundarySegments', rproblem%NboundarySegments(i), iarrayindex = i)
    end do

    ! for max boundarysegment
    imax = -1
    do i = 1,rproblem%nboundaries
	    if (rproblem%NboundarySegments(i) .gt. imax) then
        imax = rproblem%NboundarySegments(i)
      end if
    end do


    ! boundary conditions
    allocate(rproblem%Sbc(rproblem%nboundaries,imax))
    do i = 1, rproblem%nboundaries
      do j = 1,rproblem%NboundarySegments(i)
        call parlst_getvalue_string (rparams, '', 'Sbc'//trim(sys_siL(i,3)), sstring, isubstring = j)
	      read(sstring,*) rproblem%Sbc(i,j)
      end do
    end do
!     do i = 1,rproblem%nboundarySegments
! 	      call parlst_getvalue_string (rparams, '', 'Sbc'//trim(sys_siL(i,3)), sstring)
! 	      read(sstring,*) rproblem%Sbc(i)
!     end do

    !  kind of equation (possible values: POISSON, ELASTICITY)
    call parlst_getvalue_string (rparams, '', &
                                     'ctypeOfEquation',sstring)

    if(trim(sstring) .eq. 'POISSON') then
      rproblem%ctypeOfEquation = POISSON
      print *, 'POISSON'
    else if(trim(sstring) .eq. 'ELASTICITY') then
      rproblem%ctypeOfEquation = ELASTICITY
      print *, 'ELASTICITY'
    else
      print *,'Invalid sstring for equation type:', trim(sstring)
    Stop
    end if
    
    ! type of configuration (possible values: SIMUL_REAL, SIMUL_ANALYTICAL)
    call parlst_getvalue_string (rparams, '', &
                                     'ctypeOfSimulation',sstring)

    if(trim(sstring) .eq. 'ANALYTICAL') then
	    rproblem%ctypeOfSimulation = SIMUL_ANALYTICAL
      print *, 'ANALYTICAL SIMULATION'
    else if(trim(sstring) .eq. 'REAL') then
	    rproblem%ctypeOfSimulation = SIMUL_REAL
      print *, 'REAL SIMULATION'
    else
	    print *,'Invalid sstring for simulation type:', trim(sstring)
    Stop
    end if

    ! type of smoother (possible values: J_SMOOTHER, ILU_SMOOTHER)
    call parlst_getvalue_string (rparams, '', &
                                     'ctypeOfSmoother',sstring)

    if(trim(sstring) .eq. 'JACOBI') then
      rproblem%ctypeOfSmoother = J_SMOOTHER
      print *, 'JACOBI SMOOTHER'
    else if(trim(sstring) .eq. 'ILU') then
      rproblem%ctypeOfSmoother = ILU_SMOOTHER
      print *, 'ILU SMOOTHER'
    else
      print *,'Invalid sstring for smoother type:', trim(sstring)
    Stop
    end if

    !     kind of solver (possible values: DIRECT_SOLVER,BICGSTAB_SOLVER,MG_SOLVER)
    call parlst_getvalue_string (rparams, '', &
                                     'ctypeOfSolver',sstring)

    if(trim(sstring) .eq. 'DIRECT_SOLVER') then
      rproblem%ctypeOfSolver = DIRECT_SOLVER
      print *, 'DIRECT SOLVER'
    else if(trim(sstring) .eq. 'MG_SOLVER') then
      rproblem%ctypeOfSolver = MG_SOLVER
      print *, 'MG_SOLVER'
    else if(trim(sstring) .eq. 'CG_MG_SOLVER') then
      rproblem%ctypeOfSolver = CG_MG_SOLVER
      print *, 'CG_MG_SOLVER'
    else if(trim(sstring) .eq. 'CG_SOLVER') then
      rproblem%ctypeOfSolver = CG_SOLVER
      print *, 'CG_SOLVER'
    else if(trim(sstring) .eq. 'MG_CG_SOLVER') then
      rproblem%ctypeOfSolver = MG_CG_SOLVER
      print *, 'MG_CG_SOLVER'
    else if(trim(sstring) .eq. 'MG_BICG_SOLVER') then
      rproblem%ctypeOfSolver = MG_BICG_SOLVER
      print *, 'MG_BICG_SOLVER'
    else if(trim(sstring) .eq. 'BICGSTAB_SOLVER') then
      rproblem%ctypeOfSolver = BICGSTAB_SOLVER
      print *, 'BICGSTAB_SOLVER'
    else
      print *,'Invalid sstring for solver type:', trim(sstring)
    Stop
    end if

    ! Minimum level
    call parlst_getvalue_int (rparams, '', &
                                     'NLMIN', rproblem%NLMIN)
    ! Maximum level
    call parlst_getvalue_int (rparams, '', &
                                     'NLMAX', rproblem%NLMAX)

    if (rproblem%ctypeOfSolver /= MG_SOLVER .and. &
        rproblem%ctypeOfSolver /= CG_MG_SOLVER .and. &
        rproblem%ctypeOfSolver /= MG_BICG_SOLVER .and. &
        rproblem%ctypeOfSolver /= MG_CG_SOLVER) then
      rproblem%NLMIN = rproblem%NLMAX
    end if

    ! Allocate memory for all levels
    allocate(Rlevels(rproblem%NLMIN:rproblem%NLMAX))

    ! INPUT PARAMETER: Cycle identifier.
    !  0=F-cycle,
    !  1=V-cycle,
    !  2=W-cycle.
    call parlst_getvalue_int (rparams, '', &
                                    'ccycle', rproblem%ccycle)
    print *,'ccycle:', rproblem%ccycle
        
    ! function IDs (only needed in case of ctypeOfSimulation .eq. SIMUL_ANALYTICAL)
    call parlst_getvalue_int (rparams, '', &
                                     'cfuncID_u1', rproblem%cfuncID_u1)

    call parlst_getvalue_int (rparams, '', &
                                     'cfuncID_u2', rproblem%cfuncID_u2)
         
    ! deformation in gmv(possible values: Y (YES), N (NO))
    call parlst_getvalue_string (rparams, '', &
                                     'Deformation', sstring)
    if(trim(sstring) .eq. 'Y') then
	    rproblem%DEFORMATION = Y
      print *, 'Deformation Y'
    else if(trim(sstring) .eq. 'N') then
	    rproblem%DEFORMATION = N
      print *, 'Deformation N'
    else
	    print *,'Invalid sstring for Deformation type:', trim(sstring)
    Stop
    end if

    ! calculate sol on a point(possible values: Y (YES), N (NO))
    call parlst_getvalue_string (rparams, '', &
                                     'inquirePoint', sstring)
    if(trim(sstring) .eq. 'Y') then
      rproblem%inquirePoint = Y
      print *, 'inquirePoint Y'
      ! PointX
      call parlst_getvalue_double (rparams, '', &
                                     'inquirePointX', rproblem%inquirePointX)
      ! PointY
      call parlst_getvalue_double (rparams, '', &
                                     'inquirePointY', rproblem%inquirePointY)
      ! Reference sol for U1
      call parlst_getvalue_double (rparams, '', &
                                     'refSolU1', rproblem%refSolU1)
      ! Reference sol for U2
      call parlst_getvalue_double (rparams, '', &
                                     'refSolU2', rproblem%refSolU2)
    else if(trim(sstring) .eq. 'N') then
      rproblem%inquirePoint = N
      print *, 'inquirePoint N'
    else
      print *,'Invalid sstring for inquirePoint type:', trim(sstring)
    Stop
    end if

    Dpoints(1,1) = rproblem%inquirePointX
    Dpoints(2,1) = rproblem%inquirePointY

    !    max number of iterations
    call parlst_getvalue_int (rparams, '', &
                                     'niterations', rproblem%niterations)

    !    number of smoothing steps
    call parlst_getvalue_int (rparams, '', &
                                     'nsmoothingSteps', rproblem%nsmoothingSteps)

    !  tolerance
    call parlst_getvalue_double (rparams, '', &
                                     'dtolerance', rproblem%dtolerance)

    !  Damping parameter
    call parlst_getvalue_double (rparams, '', &
                                     'ddamp', rproblem%ddamp)
         
    ! material parameters (Poisson ratio nu and shear modulus mu)

    call parlst_getvalue_double (rparams, '', &
                                     'dnu', rproblem%dnu)

    call parlst_getvalue_double (rparams, '', &
                                     'dmu', rproblem%dmu)
  
    rproblem%dlambda = 2.0_DP * rproblem%dmu * rproblem%dnu/(1 - 2.0_DP * rproblem%dnu)
                     
    ! set constant RHS values (only needed in case of ctypeOfSimulation .eq. SIMUL_REAL)
    call parlst_getvalue_double (rparams, '', &
                                     'drhsVol1', rproblem%drhsVol1)

    call parlst_getvalue_double (rparams, '', &
                                     'drhsVol2', rproblem%drhsVol2)

    allocate(rproblem%DrhsBoundx(rproblem%nboundaries,imax))
    do i = 1, rproblem%nboundaries
      do j = 1,rproblem%NboundarySegments(i)
        call parlst_getvalue_string (rparams, '', 'DrhsBoundx'//trim(sys_siL(i,3)), sstring, isubstring = j)
        read(sstring,*) rproblem%DrhsBoundx(i,j)
      end do
    end do

    allocate(rproblem%DrhsBoundy(rproblem%nboundaries,imax))
    do i = 1, rproblem%nboundaries
      do j = 1,rproblem%NboundarySegments(i)
        call parlst_getvalue_string (rparams, '', 'DrhsBoundy'//trim(sys_siL(i,3)), sstring, isubstring = j)
        read(sstring,*) rproblem%DrhsBoundy(i,j)
      end do
    end do

                 
    print *,'NLMIN = ', rproblem%NLMIN
    print *,'NLMAX = ', rproblem%NLMAX
    print *,'cfuncID_u1 = ', rproblem%cfuncID_u1
    print *,'cfuncID_u2 = ', rproblem%cfuncID_u2
    print *,'dnu = ', rproblem%dnu
    print *,'dmu = ', rproblem%dmu
    print *,'dlambda = ', rproblem%dlambda
    print *,'selementType = ', rproblem%selement
    print *,'inquirePointX = ', rproblem%inquirePointX
    print *,'inquirePointY = ', rproblem%inquirePointY
   
    do i = 1, rproblem%nboundaries
      do j = 1,rproblem%NboundarySegments(i)
        call parlst_getvalue_string (rparams, '', 'Sbc'//trim(sys_siL(i,3)), sstring, isubstring = j)
        call output_line ('Sbc = ' // trim(rproblem%Sbc(i,j)) )
      end do
    end do

    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------
    !
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, rproblem%sgridFilePRM)
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (Rlevels(rproblem%NLMIN)%rtriangulation, rproblem%sgridFileTRI, rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering(rproblem%NLMIN-1, &
                                        Rlevels(rproblem%NLMIN)%rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (Rlevels(rproblem%NLMIN)%rtriangulation,rboundary)

    ! Now refine the grid for the fine levels.
    do i = rproblem%NLMIN+1, rproblem%NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
          Rlevels(i)%rtriangulation,rboundary)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
        rboundary)
    
    end do

    ! poisson equation with multigrid
    if (rproblem%ctypeOfEquation .eq. POISSON) then

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      ! Do this for all levels
      do i = rproblem%NLMIN, rproblem%NLMAX
        call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                    Rlevels(i)%rtriangulation, rboundary)
      end do
      
      ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
      ! structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      do i = rproblem%NLMIN, rproblem%NLMAX
        if (rproblem%selement .EQ. 'Q2') then
    
          call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
            EL_Q2, CUB_G3X3, Rlevels(i)%rtriangulation, rboundary)
        
        else if (rproblem%selement .EQ. 'Q1') then
 
          call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
              EL_Q1, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
      
        end if
      end do
                  
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create a scalar matrix, based on the discretisation structure
      ! for our one and only solution component.
      do i = rproblem%NLMIN, rproblem%NLMAX
  
        ! Initialise the block matrix with default values based on
        ! the discretisation.
        call lsysbl_createMatBlockByDiscr (&
            Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)
  
        ! Now as the discretisation is set up, we can start to generate
        ! the structure of the system matrix which is to solve.
        ! We create that directly in the block (1,1) of the block matrix
        ! using the discretisation structure of the first block.
        call bilf_createMatrixStructure ( &
            Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
            LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))
        
        ! And now to the entries of the matrix. For assembling of the entries,
        ! we need a bilinear form, which first has to be set up manually.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 2D.
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
  
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
        call bilf_buildMatrixScalar (rform,.true.,&
            Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D)
      
      end do
        
      ! Although we could manually create the solution/RHS vector,
      ! the easiest way to set up the vector structure is
      ! to create it by using our matrix as template:
      call lsysbl_createVecBlockIndMat (Rlevels(rproblem%NLMAX)%rmatrix,rrhs, .false.)
      call lsysbl_createVecBlockIndMat (Rlevels(rproblem%NLMAX)%rmatrix,rvector, .false.)
  
      ! The vector structure is ready but the entries are missing.
      ! So the next thing is to calculate the content of that vector.
      !
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! ... and then discretise the RHS to get a discrete version of it.
      ! Again we simply create a scalar vector based on the one and only
      ! discretisation structure.
      ! This scalar vector will later be used as the one and only first
      ! component in a block vector.
      call linf_buildVectorScalar (Rlevels(rproblem%NLMAX)%rdiscretisation%RspatialDiscr(1),&
                    rlinform,.true.,rrhs%RvectorBlock(1),coeff_RHS_Vol_u1_2D)

      ! print # of DOF
      call lsysbl_getbase_double (rVector,p_Ddata)
  
      call output_line ('Number of DOF: ' // trim(sys_siL(size(p_Ddata),12)) )
  
      ! Clear the solution vector on the finest level.
      call lsysbl_clearVector(rvector)
  
      
  
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponsolidmechding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      do i = rproblem%NLMIN, rproblem%NLMAX
        call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)
      end do


      do j = 1, rproblem%nboundaries
        do k = 1,rproblem%NboundarySegments(j)
          call boundary_createRegion(rboundary,j,k,rboundaryRegion)
          if (rproblem%Sbc(j,k) .eq. 'D') then
            ! The endpoint of this segment should also be Dirichlet. We set this by
            ! changing the region properties in rboundaryRegion.
            rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
            print *,'D',j,k
            do i = rproblem%NLMIN, rproblem%NLMAX
              ! We use this boundary region and specify that we want to have Dirichlet
              ! boundary there. The following call does the following:
              ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
              ! We specify icomponent='1' to indicate that we set up the
              ! Dirichlet BC's for the first (here: one and only) component in the
              ! solution vector.
              ! - Discretise the boundary condition so that the BC's can be applied
              ! to matrices and vectorssolidmech
              ! - Add the calculated discrete BC's to rdiscreteBC for later use.
              call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                getBoundaryValues_2D)

            end do
          else if (rproblem%Sbc(j,k) .eq. 'N') then
            rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
            print *,'N',j,k
            rcollection%IquickAccess(1) = k
            ! For Element Q1
            if (rproblem%selement .EQ. 'Q1') then
  
              ! We use this boundary region and specify that we want to have Neumann boundary there
              ! and add this to volumetric part(the neumann boundary part on rhs) for U1
              call linf_buildVectorScalarBdr2d(rlinform,CUB_G2_1D,.false.,&
                rrhs%RvectorBlock(1),coeff_RHS_neumBdr_u1_2D,rboundaryRegion,rcollection)
  
            ! For Element Q2
            else if (rproblem%selement .EQ. 'Q2') then
  
              ! For U1
              call linf_buildVectorScalarBdr2d(rlinform,CUB_G3_1D,.false.,&
                rrhs%RvectorBlock(1),coeff_RHS_neumBdr_u1_2D,rboundaryRegion,rcollection)

            end if
          else
            print *, ' Invalid Input for Boundary Condition'
            stop
          end if
        end do ! end segments
      end do ! end boundaries
      

    ! elsticity equations
    else if (rproblem%ctypeOfEquation .eq. ELASTICITY) then
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we have two blocks.
      do i = rproblem%NLMIN, rproblem%NLMAX
        call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 2, &
                                    Rlevels(i)%rtriangulation, rboundary)
      end do
      
      ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
      ! structures for every component of the solution vector.
      ! We have a solution vector with two components:
      do i = rproblem%NLMIN, rproblem%NLMAX
        if (rproblem%selement .EQ. 'Q2') then
    
	        !  Component 1 = X-velocity
	        call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
			      EL_Q2, CUB_G3X3, Rlevels(i)%rtriangulation, rboundary)
	      
	        ! Component 2 = Y-velocity
	        call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(2),&
			        EL_Q2, CUB_G3X3, Rlevels(i)%rtriangulation, rboundary)
	      
	      else if (rproblem%selement .EQ. 'Q1') then
	        ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
	        ! structures for every component of the solution vector.
	        ! We have a solution vector with two components:
	        !  Component 1 = X-velocity
	        call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
			        EL_Q1, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
	    
	        ! Component 2 = Y-velocity
	        call spdiscr_initDiscr_simple (Rlevels(i)%rdiscretisation%RspatialDiscr(2),&
			        EL_Q1, CUB_G2X2, Rlevels(i)%rtriangulation, rboundary)
        end if
      end do
  
      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      do i = rproblem%NLMIN, rproblem%NLMAX
  
        ! Initialise the block matrix with default values based on
        ! the discretisation.
        call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)
        
        ! Now as the discretisation is set up, we can start to generate
        ! the structure of the system matrix which is to solve.
        ! We create that directly in the block (1,1) of the block matrix
        ! using the discretisation structure of the first block.
        !
        ! Create the matrix structure of the X-velocity.
        call bilf_createMatrixStructure (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9, Rlevels(i)%rmatrix%RmatrixBlock(1,1))
        
        ! And now to the entries of the matrix. For assembling of the entries,
        ! we need a bilinear form, which first has to be set up manually.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 2D.
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
    
        ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = 2*rproblem%dmu + rproblem%dlambda
        rform%Dcoefficients(2)  = rproblem%dmu
                  
        ! Now we can build the matrix entries.
        ! We specify the callback function coeff_Laplace for the coefficients.
        ! As long as we use constant coefficients, this routine is not used.
        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
        ! the framework will call the callback routine to get analytical
        ! data.
        call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                                    coeff_Laplace_2D)
                  
        ! Now We create the block (1,2) of the block matrix
        ! using the discretisation structure of the first block.
        call bilf_createMatrixStructure (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9, Rlevels(i)%rmatrix%RmatrixBlock(1,2))
        
        ! And now to the entries of the matrix. For assembling of the entries,
        ! we need a bilinear form, which first has to be set up manually.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 2D.
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_Y
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_X
        rform%Idescriptors(2,2) = DER_DERIV_Y
    
        ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = rproblem%dlambda
        rform%Dcoefficients(2)  = rproblem%dmu
        ! Now we can build the matrix entries.
        ! We specify the callback function coeff_Laplace for the coefficients.
        ! As long as we use constant coefficients, this routine is not used.
        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
        ! the framework will call the callback routine to get analytical
        ! data.
        call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrix%RmatrixBlock(1,2),&
                                    coeff_Laplace_2D)
                  
        ! Now We create the block (2,1) of the block matrix
        ! using the discretisation structure of the first block.
        call bilf_createMatrixStructure (Rlevels(i)%rdiscretisation%RspatialDiscr(2),&
                                        LSYSSC_MATRIX9, Rlevels(i)%rmatrix%RmatrixBlock(2,1))
        
        ! And now to the entries of the matrix. For assembling of the entries,
        ! we need a bilinear form, which first has to be set up manually.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 2D.
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_Y
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_X
    
        ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = rproblem%dlambda
        rform%Dcoefficients(2)  = rproblem%dmu
                  
        ! Now we can build the matrix entries.
        ! We specify the callback function coeff_Laplace for the coefficients.
        ! As long as we use constant coefficients, this routine is not used.
        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
        ! the framework will call the callback routine to get analytical
        ! data.
        call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrix%RmatrixBlock(2,1),&
                                    coeff_Laplace_2D)
                    
        ! Now We create the block (2,2) of the block matrix
        ! using the discretisation structure of the first block.
        call bilf_createMatrixStructure (Rlevels(i)%rdiscretisation%RspatialDiscr(2),&
                                        LSYSSC_MATRIX9, Rlevels(i)%rmatrix%RmatrixBlock(2,2))
        
        ! And now to the entries of the matrix. For assembling of the entries,
        ! we need a bilinear form, which first has to be set up manually.
        ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
        ! scalar system matrix in 2D.
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
    
        ! In the standard case, we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = rproblem%dmu
        rform%Dcoefficients(2)  = 2*rproblem%dmu + rproblem%dlambda
                  
        ! Now we can build the matrix entries.
        ! We specify the callback function coeff_Laplace for the coefficients.
        ! As long as we use constant coefficients, this routine is not used.
        ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
        ! the framework will call the callback routine to get analytical
        ! data.
        call bilf_buildMatrixScalar (rform,.true.,Rlevels(i)%rmatrix%RmatrixBlock(2,2),&
                                    coeff_Laplace_2D)
  
      end do
                
      ! Although we could manually create the solution/RHS vector,
      ! the easiest way to set up the vector structure is
      ! to create it by using our matrix as template:
      call lsysbl_createVecBlockIndMat (Rlevels(rproblem%NLMAX)%rmatrix,rrhs, .false.)
      call lsysbl_createVecBlockIndMat (Rlevels(rproblem%NLMAX)%rmatrix,rvector, .false.)
  
      ! The vector structure is ready but the entries are missing.
      ! So the next thing is to calculate the content of that vector.
      ! At first set up the corresponding linear form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! Clear the solution vector on the finest level.
      call lsysbl_clearVector(rrhs)
  
      ! ... and then discretise the RHS to the first two subvectors of
      ! the block vector using the discretisation structure of the
      ! corresponding blocks.
      !
      ! Note that the vector is unsorted after calling this routine!
      ! For U1
      call linf_buildVectorScalar (Rlevels(rproblem%NLMAX)%rdiscretisation%RspatialDiscr(1),&
                    rlinform,.true.,rrhs%RvectorBlock(1),coeff_RHS_Vol_u1_2D)
  
      ! For U2
      call linf_buildVectorScalar (Rlevels(rproblem%NLMAX)%rdiscretisation%RspatialDiscr(2),&
                    rlinform,.true.,rrhs%RvectorBlock(2),coeff_RHS_Vol_u2_2D)
  
      ! print # of DOF
      call lsysbl_getbase_double (rVector,p_Ddata)
  
      call output_line ('Number of DOF: ' // trim(sys_siL(size(p_Ddata),12)) )
  
      ! Clear the solution vector on the finest level.
      call lsysbl_clearVector(rvector)
  
      
  
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponsolidmechding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      do i = rproblem%NLMIN, rproblem%NLMAX
        call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)
      end do
    
      ! We first set up the boundary conditions for the U1-velocity, then those
      ! of the U2-velocity.
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the bondary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      do j = 1, rproblem%nboundaries
        do k = 1,rproblem%NboundarySegments(j)
	        call boundary_createRegion(rboundary,j,k,rboundaryRegion)
	        if (rproblem%Sbc(j,k) .eq. 'D') then
            ! The endpoint of this segment should also be Dirichlet. We set this by
            ! changing the region properties in rboundaryRegion.
		        rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
            print *,'D',j,k
            do i = rproblem%NLMIN, rproblem%NLMAX
              ! We use this boundary region and specify that we want to have Dirichlet
              ! boundary there. The following call does the following:
              ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
              ! We specify icomponent='1' to indicate that we set up the
              ! Dirichlet BC's for the first (here: one and only) component in the
              ! solution vector.
              ! - Discretise the boundary condition so that the BC's can be applied
              ! to matrices and vectorssolidmech
              ! - Add the calculated discrete BC's to rdiscreteBC for later use.
		          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
						    rboundaryRegion,Rlevels(i)%rdiscreteBC,&
						    getBoundaryValues_2D)
    
              ! Now continue with defining the boundary conditions of the U2-velocity:
		          call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,2,&
						    rboundaryRegion,Rlevels(i)%rdiscreteBC,&
						    getBoundaryValues_2D)
	          end do
          else if (rproblem%Sbc(j,k) .eq. 'N') then
		        rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
            print *,'N',j,k
            rcollection%IquickAccess(1) = k
            ! For Element Q1
		        if (rproblem%selement .EQ. 'Q1') then
  
              ! We use this boundary region and specify that we want to have Neumann boundary there
              ! and add this to volumetric part(the neumann boundary part on rhs) for U1
			        call linf_buildVectorScalarBdr2d(rlinform,CUB_G2_1D,.false.,&
					      rrhs%RvectorBlock(1),coeff_RHS_neumBdr_u1_2D,rboundaryRegion,rcollection)
  
              ! Now its done for U2
			        call linf_buildVectorScalarBdr2d(rlinform,CUB_G2_1D,.false.,&
					      rrhs%RvectorBlock(2),coeff_RHS_neumBdr_u2_2D,rboundaryRegion,rcollection)
  
            ! For Element Q2
		        else if (rproblem%selement .EQ. 'Q2') then
  
              ! For U1
			        call linf_buildVectorScalarBdr2d(rlinform,CUB_G3_1D,.false.,&
					      rrhs%RvectorBlock(1),coeff_RHS_neumBdr_u1_2D,rboundaryRegion,rcollection)
  
              ! For U2
			        call linf_buildVectorScalarBdr2d(rlinform,CUB_G3_1D,.false.,&
					      rrhs%RvectorBlock(2),coeff_RHS_neumBdr_u2_2D,rboundaryRegion,rcollection)
		        end if
	        else
	          print *, ' Invalid Input for Boundary Condition'
	          stop
	        end if
        end do ! end segments
      end do ! end boundaries
    end if ! elasticity equation

    do i = rproblem%NLMIN, rproblem%NLMAX
      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix.
      Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
  
      ! Also implement the boundary conditions into the matrix.
      call matfil_discreteBC (Rlevels(i)%rmatrix)

    end do

 
    ! Hang the pointer into the vector. That way, these
    ! boundary conditions are always connected to that
    ! vector.!
    rrhs%p_rdiscreteBC => Rlevels(rproblem%NLMAX)%rdiscreteBC
    rvector%p_rdiscreteBC => Rlevels(rproblem%NLMAX)%rdiscreteBC

!     if (rproblem%ctypeOfEquation .eq. POISSON) then
!       ! Now we have block vectors for the RHS and the matrix. What we
!       ! need additionally is a block vector for the solution and
!       ! temporary data. Create them using the RHS as template.
!       ! Fill the solution vector with 0:
!       call lsysbl_createVecBlockIndirect (rrhs, rvector, .true.)
!       call lsysbl_createVecBlockIndirect (rrhs, rtempBlock, .false.)
!     end if
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhs)
    call vecfil_discreteBCsol (rvector)

    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver.
    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
    if (rproblem%ctypeOfSolver .eq. BICGSTAB_SOLVER) then
      p_RfilterChain => RfilterChain
      nullify(p_rpreconditioner)
      call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)

    ! Create Conjugate solver
    else if (rproblem%ctypeOfSolver .eq. CG_SOLVER) then
      p_RfilterChain => RfilterChain
      nullify(p_rpreconditioner)
      call linsol_initCG (p_rsolverNode,p_rpreconditioner,p_RfilterChain)

    ! For Direct solver we use following
    else if (rproblem%ctypeOfSolver .eq. DIRECT_SOLVER) then
      call linsol_initUMFPACK4 (p_rsolverNode)
      p_rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput = 0

    ! Now we have to build up the level information for multigrid.

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    else if (rproblem%ctypeOfSolver .eq. MG_SOLVER) then
      p_RfilterChain => RfilterChain
      call linsol_initMultigrid2 (p_rsolverNode,rproblem%NLMAX-rproblem%NLMIN+1,p_RfilterChain)

      ! Set up a coarse grid solver.
      ! The coarse grid in multigrid is always grid 1!
      call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
      call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
      
      ! Now set up the other levels...
      do i = rproblem%NLMIN+1, rproblem%NLMAX
      
        ! Create a Jacobi smoother
        if (rproblem%ctypeOfSmoother .eq. J_SMOOTHER) then
          call linsol_initJacobi(p_rsmoother)
        end if

        ! Create an ILU(0) smoother
        if (rproblem%ctypeOfSmoother .eq. ILU_SMOOTHER) then
          call linsol_initILU0(p_rsmoother)
        end if

        ! We will use nsmoothingSteps smoothing steps with damping parameter ddamp
        call linsol_convertToSmoother(p_rsmoother, rproblem%nsmoothingSteps, rproblem%ddamp)
        
        ! And add this multi-grid level. We will use the same smoother
        ! for pre- and post-smoothing.
        call linsol_getMultigrid2Level (p_rsolverNode,i-rproblem%NLMIN+1,p_rlevelInfo)
        p_rlevelInfo%p_rpresmoother => null()
        p_rlevelInfo%p_rpostsmoother => p_rsmoother
        
      end do

      ! INPUT PARAMETER: Cycle identifier.
      !  0=F-cycle,
      !  1=V-cycle,
      !  2=W-cycle.
      p_rsolverNode%p_rsubnodeMultigrid2%icycle = rproblem%ccycle
      
    else if (rproblem%ctypeOfSolver .eq. CG_MG_SOLVER) then
      p_RfilterChain => RfilterChain
      call linsol_initMultigrid2 (p_rpreconditioner,rproblem%NLMAX-rproblem%NLMIN+1,p_RfilterChain)

      ! Set up a coarse grid solver.
      ! The coarse grid in multigrid is always grid 1!
      call linsol_getMultigrid2Level (p_rpreconditioner,1,p_rlevelInfo)
      call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
      
      ! Now set up the other levels...
      do i = rproblem%NLMIN+1, rproblem%NLMAX
      
        ! Create a Jacobi smoother
        if (rproblem%ctypeOfSmoother .eq. J_SMOOTHER) then
          call linsol_initJacobi(p_rsmoother)
        end if

        ! Create an ILU(0) smoother
        if (rproblem%ctypeOfSmoother .eq. ILU_SMOOTHER) then
          call linsol_initILU0(p_rsmoother)
        end if

        ! We will use nsmoothingSteps smoothing steps with damping parameter ddamp
        call linsol_convertToSmoother(p_rsmoother, rproblem%nsmoothingSteps, rproblem%ddamp)
        
        ! And add this multi-grid level. We will use the same smoother
        ! for pre- and post-smoothing.
        call linsol_getMultigrid2Level (p_rpreconditioner,i-rproblem%NLMIN+1,p_rlevelInfo)
        p_rlevelInfo%p_rpresmoother => p_rsmoother
        p_rlevelInfo%p_rpostsmoother => p_rsmoother
        
      end do

      ! INPUT PARAMETER: Cycle identifier.
      !  0=F-cycle,
      !  1=V-cycle,
      !  2=W-cycle.
      p_rpreconditioner%p_rsubnodeMultigrid2%icycle = rproblem%ccycle

      ! Set the output level of the solver to 2 for some output
      p_rpreconditioner%ioutputLevel = 1
  
      ! We will allow the solver to perform 1 iteration
      p_rpreconditioner%nmaxIterations = 1
  
      p_rpreconditioner%depsRel = 1.0E-99_DP

      call linsol_initCG (p_rsolverNode,p_rpreconditioner,p_RfilterChain)

    else if (rproblem%ctypeOfSolver .eq. MG_CG_SOLVER .or. &
             rproblem%ctypeOfSolver .eq. MG_BICG_SOLVER ) then
      p_RfilterChain => RfilterChain
      nullify(p_rpreconditioner)
      ! Create a Multigrid-solver. Attach the above filter chain
      ! to the solver, so that the solver automatically filters
      ! the vector during the solution process.
      call linsol_initMultigrid2 (p_rsolverNode,rproblem%NLMAX-rproblem%NLMIN+1,p_RfilterChain)

      ! Set up a coarse grid solver.
      ! The coarse grid in multigrid is always grid 1!
      call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
      call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
      
      ! Now set up the other levels...
      do i = rproblem%NLMIN+1, rproblem%NLMAX
      
        ! Create a Jacobi smoother
        if (rproblem%ctypeOfSmoother .eq. J_SMOOTHER) then
          call linsol_initJacobi(p_rpreconditioner)
        end if

        ! Create an ILU(0) smoother
        if (rproblem%ctypeOfSmoother .eq. ILU_SMOOTHER) then
          call linsol_initILU0(p_rpreconditioner)
        end if

        ! Create Conjugate solver
        if (rproblem%ctypeOfSolver .eq. MG_CG_SOLVER) then
          call linsol_initCG (p_rsmoother,p_rpreconditioner,p_RfilterChain)
        ! Create BiConjugate solver
        else if (rproblem%ctypeOfSolver .eq. MG_BICG_SOLVER) then
          call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,p_RfilterChain)
        end if

        ! We will use nsmoothingSteps smoothing steps with damping parameter ddamp
        call linsol_convertToSmoother(p_rsmoother, rproblem%nsmoothingSteps, rproblem%ddamp)
 
        ! And add this multi-grid level. We will use the same smoother
        ! for pre- and post-smoothing.
        call linsol_getMultigrid2Level (p_rsolverNode,i-rproblem%NLMIN+1,p_rlevelInfo)
        p_rlevelInfo%p_rpresmoother => p_rsmoother
        p_rlevelInfo%p_rpostsmoother => p_rsmoother

        if (i .eq. rproblem%NLMAX) then
          ! Set the output level of the solver to 1 for some output
          p_rsmoother%ioutputLevel = 1
        else
          p_rsmoother%ioutputLevel = 0
        end if
        
      end do

      ! INPUT PARAMETER: Cycle identifier.
      !  0=F-cycle,
      !  1=V-cycle,
      !  2=W-cycle.
      p_rsolverNode%p_rsubnodeMultigrid2%icycle = rproblem%ccycle

    else
      print *, ' Invalid Input for Solver type'
    stop
    end if
   
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! set last 3 residuals for asymptotic rate of convergence
    p_rsolverNode%niteAsymptoticCVR = 3

    ! We will allow the solver to perform 200 iterations
    p_rsolverNode%nmaxIterations = rproblem%niterations

    p_rsolverNode%depsRel = rproblem%dtolerance
!     p_rsolverNode%depsAbs = 1E-10_DP
    print *,'tolerance',p_rsolverNode%depsRel

    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    allocate(Rmatrices(rproblem%NLMIN:rproblem%NLMAX))
    do i = rproblem%NLMIN, rproblem%NLMAX
      call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(rproblem%NLMIN:rproblem%NLMAX))

    ! We can release Rmatrices immediately -- as long as we don't
    ! release Rlevels(i)%rmatrix!
    do i = rproblem%NLMIN,rproblem%NLMAX
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
    call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)

    ! Rate of convergence
    call output_line ('Rate of Convergence: ' // sys_sdEL(p_rsolverNode%dconvergenceRate,10) )

    ! Rate of asymptotic convergence
    call output_line ('Rate of Asymptotic Convergence: ' // &
                      sys_sdEL(p_rsolverNode%dasymptoticConvergenceRate,10) )

    ! Calculate sol.U1(x,y) & U2(x,y) on a point (x,y) and their
    ! derivatives U1_x(x,y), U1_y(x,y) & U2_x(x,y), U2_y(x,y)
    ! and absolute error , Strain, Stress , plain srain
    if (rproblem%inquirePoint .eq. Y) then
      call fevl_evaluate (DER_FUNC, Dvalues, rvector%RvectorBlock(1), Dpoints)
      U1 = Dvalues(1)
      call output_line ('U1(X,Y): ' // sys_sdEL(U1,10) )
      call fevl_evaluate (DER_FUNC, Dvalues, rvector%RvectorBlock(2), Dpoints)
      U2 = Dvalues(1)
      call output_line ('U2(X,Y): ' // sys_sdEL(U2,10) )
      call fevl_evaluate (DER_DERIV_X, Dvalues, rvector%RvectorBlock(1), Dpoints)
      U1_X = Dvalues(1)
      call output_line ('U1_X(X,Y): ' // sys_sdEL(U1_X,10) )
      call fevl_evaluate (DER_DERIV_X, Dvalues, rvector%RvectorBlock(2), Dpoints)
      U2_X = Dvalues(1)
      call output_line ('U2_X(X,Y): ' // sys_sdEL(U2_X,10) )
      call fevl_evaluate (DER_DERIV_Y, Dvalues, rvector%RvectorBlock(1), Dpoints)
      U1_Y = Dvalues(1)
      call output_line ('U1_Y(X,Y): ' // sys_sdEL(U1_Y,10) )
      call fevl_evaluate (DER_DERIV_Y, Dvalues, rvector%RvectorBlock(2), Dpoints)
      U2_Y = Dvalues(1)
      call output_line ('U2_Y(X,Y): ' // sys_sdEL(U2_Y,10) )
      call output_line ('eps11: ' // sys_sdEL(U1_X,10) )
      call output_line ('eps22: ' // sys_sdEL(U2_Y,10) )
      call output_line ('eps12: ' // sys_sdEL(0.5_DP*(U2_X+U1_Y),10) )
      sigma11 = 2.0_DP*rproblem%dmu*U1_X+rproblem%dlambda*(U1_X+U2_Y)
      sigma22 = 2.0_DP*rproblem%dmu*U2_Y+rproblem%dlambda*(U1_X+U2_Y)
      sigma12 = rproblem%dmu*(U2_X+U1_Y)
      sigma33 = rproblem%dlambda*(U1_X+U2_Y)
      devsigma11 = 1.0_DP/3.0_DP*(2.0_DP*sigma11-sigma22-sigma33)
      devsigma22 = 1.0_DP/3.0_DP*(2.0_DP*sigma22-sigma11-sigma33)
      devsigma12 = sigma12
      devsigma33 = 1.0_DP/3.0_DP*(2.0_DP*sigma33-sigma22-sigma11)
      call output_line ('sigma11: ' // sys_sdEL(sigma11,10) )
      call output_line ('sigma22: ' // sys_sdEL(sigma22,10) )
      call output_line ('sigma12: ' // sys_sdEL(sigma12,10) )
      call output_line ('sigma33: ' // sys_sdEL(sigma33,10) )
      call output_line ('|devsigma|: ' // sys_sdEL(sqrt(devsigma11**2+devsigma22**2+ &
                     devsigma33**2+2.0_DP*devsigma12**2),10) )
      call output_line ('Abs-error for U: ' // sys_sdEL(sqrt((rproblem%refSolU1-U1)**2+ &
                    (rproblem%refSolU2-U2)**2)/sqrt((rproblem%refSolU1)**2+(rproblem%refSolU2)**2),10) )
    end if

    ! Calculate the error to the reference function.

    if (rproblem%ctypeOfSimulation .eq. SIMUL_ANALYTICAL) then
	    call pperr_scalar (rVector%RvectorBlock(1),PPERR_L2ERROR,derror1,&
			  getReferenceFunction_u1_2D)
! 	call output_line ('L2-error for U1: ' // sys_sdEL(derror1,10) )
	
	    call pperr_scalar (rVector%RvectorBlock(2),PPERR_L2ERROR,derror2,&
			  getReferenceFunction_u2_2D)
! 	call output_line ('L2-error for U2: ' // sys_sdEL(derror2,10) )

	    call output_line ('L2-error for U: ' // sys_sdEL(sqrt(derror1*derror1 + derror2*derror2),10) )
	
	    call pperr_scalar (rVector%RvectorBlock(1),PPERR_H1ERROR,derror1,&
			  getReferenceFunction_u1_2D)
! 	call output_line ('H1-error for U1: ' // sys_sdEL(derror1,10) )
	
	    call pperr_scalar (rVector%RvectorBlock(2),PPERR_H1ERROR,derror2,&
			  getReferenceFunction_u2_2D)
! 	call output_line ('H1-error for U2: ' // sys_sdEL(derror2,10) )

	    call output_line ('H1-error for U: ' // sys_sdEL(sqrt(derror1*derror1 + derror2*derror2),10) )
    end if

    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! For Bending in the gmv
    call storage_getbase_double2D(Rlevels(rproblem%NLMAX)%rtriangulation%h_DvertexCoords, &
                                  p_DvertexCoords)

    ! Write velocity field
    call lsyssc_getbase_double (rVector%RvectorBlock(1),p_Ddata)
    if (rproblem%ctypeOfEquation .eq. ELASTICITY) then
      call lsyssc_getbase_double (rVector%RvectorBlock(2),p_Ddata2)
    end if
  
    ! we add sol. vector to coordinates to see the bending
    if (rproblem%DEFORMATION .eq. Y) then
	    do i = 1,Rlevels(rproblem%NLMAX)%rtriangulation%NVT
		    p_Dvertexcoords(1,i) = p_Dvertexcoords(1,i) + p_Ddata(i)
        if (rproblem%ctypeOfEquation .eq. ELASTICITY) then
  		    p_Dvertexcoords(2,i) = p_dvertexCoords(2,i) + p_Ddata2(i)
        end if
	    end do
    end if

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing.
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,Rlevels(rproblem%NLMAX)%rtriangulation,&
        trim(sucddir)//'/u2d_0_simple.gmv')

    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    if (rproblem%ctypeOfEquation .eq. ELASTICITY) then
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
    else if (rproblem%ctypeOfEquation .eq. POISSON) then
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata)
    end if
    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
        
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rrhs)
    do i = rproblem%NLMAX, rproblem%NLMIN, -1
      call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    end do

    ! Release our discrete version of the boundary conditions
    do i = rproblem%NLMAX, rproblem%NLMIN, -1
      call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    do i = rproblem%NLMAX, rproblem%NLMIN, -1
      call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
    end do
    
    ! Release the triangulation.
    do i = rproblem%NLMAX, rproblem%NLMIN, -1
      call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    deallocate(Rlevels)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)

  end subroutine

end module

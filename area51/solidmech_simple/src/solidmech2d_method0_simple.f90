
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
    
  use solidmech2d_callback
  
  implicit none

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
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock

    ! A set of variables describing the discrete boundary conditions.    
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

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
    
    ! Error of FE function to reference function
    real(DP) :: derror1, derror2
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir, smaster,sstring
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
    real(DP), Dimension(:,:), pointer :: p_DvertexCoords
    integer :: i

    ! Structure for saving parameters from the DAT file
    type (t_parlist) :: rparams
    
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
                              
!       boundary conditions
    allocate(rproblem%Sbc(rproblem%nboundarySegments))
    call parlst_getvalue_string (rparams, '', &
                                 'Sbc1', sstring)
    read(sstring,*) rproblem%Sbc(1)
    call parlst_getvalue_string (rparams, '', &
                                 'Sbc2', sstring)
    read(sstring,*) rproblem%Sbc(2)
    call parlst_getvalue_string (rparams, '', &
                                 'Sbc3', sstring)
    read(sstring,*) rproblem%Sbc(3)
    call parlst_getvalue_string (rparams, '', &
                                 'Sbc4', sstring)
    read(sstring,*) rproblem%Sbc(4)
       
    ! Maximum level
    call parlst_getvalue_int (rparams, '', &
                                     'NLMAX', rproblem%NLMAX)
       
!       number of boundary segments (here, one has to know the number of
!       boundary segments in the current grid)
    call parlst_getvalue_int (rparams, '', &
                                     'nboundarySegments', rproblem%nboundarySegments)
    
!       type of configuration (possible values: SIMUL_REAL, SIMUL_ANALYTICAL)
    call parlst_getvalue_string (rparams, '', &
                                     'ctypeOfSimulation',sstring)
   if(trim(sstring) .eq. 'ANALYTICAL') then
	rproblem%ctypeOfSimulation = SIMUL_ANALYTICAL
   else if(trim(sstring) .eq. 'REAL') then
	rproblem%ctypeOfSimulation = SIMUL_REAL
   else
	print *,'Invalid sstring for simulation type'
   end if
        
!      function IDs (only needed in case of ctypeOfSimulation .eq. SIMUL_ANALYTICAL)
    call parlst_getvalue_int (rparams, '', &
                                     'cfuncID_u1', rproblem%cfuncID_u1)

    call parlst_getvalue_int (rparams, '', &
                                     'cfuncID_u2', rproblem%cfuncID_u2)
         
!      deformation in gmv(possible values: ON, OFF)
    call parlst_getvalue_string (rparams, '', &
                                     'Deformation', sstring)
   if(trim(sstring) .eq. 'ON') then
	rproblem%DEFORMATION = ON
   else if(trim(sstring) .eq. 'OFF') then
	rproblem%DEFORMATION = OFF
   else
	print *,'Invalid sstring for simulation type'
   end if
         
!       material parameters (Poisson ratio nu and shear modulus mu)
    call parlst_getvalue_double (rparams, '', &
                                     'dnu', rproblem%dnu)

    call parlst_getvalue_double (rparams, '', &
                                     'dmu', rproblem%dmu)
  rproblem%dlambda = 2.0_DP * rproblem%dmu * rproblem%dnu/(1 - 2.0_DP * rproblem%dnu)
                     
!      set constant RHS values (only needed in case of ctypeOfSimulation .eq. SIMUL_REAL)
    call parlst_getvalue_double (rparams, '', &
                                     'drhsVol1', rproblem%drhsVol1)

    call parlst_getvalue_double (rparams, '', &
                                     'drhsVol2', rproblem%drhsVol2)

    call parlst_getvalue_double (rparams, '', &
                                     'drhsBound1', rproblem%drhsBound1)

    call parlst_getvalue_double (rparams, '', &
                                     'drhsBound2', rproblem%drhsBound2)
                 
! !     set grid file
!   rproblem%sgridFileTri = 'cook.tri'
!   rproblem%sgridFilePrm = 'cook.prm'
! 
! !       set number of boundary segments (here, one has to know the number of
! !       boundary segments in the current grid)
!   rproblem%nboundarySegments = 4
! 
! !       set boundary conditions
!   allocate(rproblem%Sbc(rproblem%nboundarySegments))
!   rproblem%Sbc(1) = 'N'
!   rproblem%Sbc(2) = 'N'
!   rproblem%Sbc(3) = 'N'
!   rproblem%Sbc(4) = 'D'
! 
! !      set element(possible values: Q1, Q2)
!   rproblem%selement = 'Q1'
! 
! !      set refinement level
!   rproblem%NLMAX = 1
! 
! !       material parameters (Poisson ratio nu and shear modulus mu)
!   rproblem%dnu     = 0.3_DP
!   rproblem%dmu     = 0.5_DP
!   rproblem%dlambda = 2.0_DP * rproblem%dmu * rproblem%dnu/(1 - 2.0_DP * rproblem%dnu)
! 
! !       set type of configuration (possible values: SIMUL_REAL, SIMUL_ANALYTICAL)
!   rproblem%ctypeOfSimulation = SIMUL_ANALYTICAL
! 
! !      set function IDs (only needed in case of ctypeOfSimulation .eq. SIMUL_ANALYTICAL)
!   rproblem%cfuncID_u1 = 53
!   rproblem%cfuncID_u2 = 52
! 
! !      set constant RHS values (only needed in case of ctypeOfSimulation .eq. SIMUL_REAL)
!   rproblem%drhsVol1   = 0.0_DP
!   rproblem%drhsVol2   = -0.001_DP
!   rproblem%drhsBound1 = 0.0_DP
!   rproblem%drhsBound2 = 0.0_DP
! 
! !      set deformation(possible values: ON, OFF)
!   rproblem%DEFORMATION = ON

    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------
    !
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, rproblem%sgridFilePRM)
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, rproblem%sgridFileTRI, rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (rproblem%NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we have two blocks.
    call spdiscr_initBlockDiscr (rdiscretisation,2,&
                                 rtriangulation, rboundary)
    
    if (rproblem%selement .EQ. 'Q2') then
	! rdiscretisation%Rdiscretisations is a list of scalar discretisation
	! structures for every component of the solution vector.
	! We have a solution vector with two components:
	!  Component 1 = X-velocity
	call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
			EL_Q2, CUB_G3X3, rtriangulation, rboundary)
	
	! Component 2 = Y-velocity
	call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(2),&
			EL_Q2, CUB_G3X3, rtriangulation, rboundary)
	
	else if (rproblem%selement .EQ. 'Q1') then
	! rdiscretisation%Rdiscretisations is a list of scalar discretisation
	! structures for every component of the solution vector.
	! We have a solution vector with two components:
	!  Component 1 = X-velocity
	call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
			EL_Q1, CUB_G2X2, rtriangulation, rboundary)
	
	! Component 2 = Y-velocity
	call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(2),&
			EL_Q1, CUB_G2X2, rtriangulation, rboundary)
    end if

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)
    
    ! Inform the matrix that we build a saddle-point problem.
    ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
    ! but probably some solvers can use the special structure later.
    rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
    
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    !
    ! Create the matrix structure of the X-velocity.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,1))
    
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
   call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(1,1),&
                                 coeff_Laplace_2D)
              
    ! Now We create the block (1,2) of the block matrix
    ! using the discretisation structure of the first block.
     call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,2))
    
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
   call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(1,2),&
                                 coeff_Laplace_2D)
              
    ! Now We create the block (2,1) of the block matrix
    ! using the discretisation structure of the first block.
     call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(2),&
                                     LSYSSC_MATRIX9, rmatrix%RmatrixBlock(2,1))
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_Y
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_X

!     In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = rproblem%dlambda
    rform%Dcoefficients(2)  = rproblem%dmu
              
!     Now we can build the matrix entries.
!     We specify the callback function coeff_Laplace for the coefficients.
!     As long as we use constant coefficients, this routine is not used.
!     By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!     the framework will call the callback routine to get analytical
!     data.
   call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(2,1),&
                                 coeff_Laplace_2D)
                
!     Now We create the block (2,2) of the block matrix
!     using the discretisation structure of the first block.
     call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(2),&
                                     LSYSSC_MATRIX9, rmatrix%RmatrixBlock(2,2))
    
!     And now to the entries of the matrix. For assembling of the entries,
!     we need a bilinear form, which first has to be set up manually.
!     We specify the bilinear form (grad Psi_j, grad Phi_i) for the
!     scalar system matrix in 2D.
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y

!     In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = rproblem%dmu
    rform%Dcoefficients(2)  = 2*rproblem%dmu + rproblem%dlambda
              
!     Now we can build the matrix entries.
!     We specify the callback function coeff_Laplace for the coefficients.
!     As long as we use constant coefficients, this routine is not used.
!     By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!     the framework will call the callback routine to get analytical
!     data.
   call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(2,2),&
                                 coeff_Laplace_2D)
              
!     Although we could manually create the solution/RHS vector,
!     the easiest way to set up the vector structure is
!     to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (rmatrix,rrhs, .false.)
    call lsysbl_createVecBlockIndMat (rmatrix,rvector, .false.)

!     The vector structure is ready but the entries are missing. 
!     So the next thing is to calculate the content of that vector.
!     
!     At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
!     Clear the solution vector on the finest level.
    call lsysbl_clearVector(rrhs)

!     ... and then discretise the RHS to the first two subvectors of
!     the block vector using the discretisation structure of the 
!     corresponding blocks.
!     
!     Note that the vector is unsorted after calling this routine!
!       For U1
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                  rlinform,.true.,rrhs%RvectorBlock(1),coeff_RHS_Vol_u1_2D)

!       For U2
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(2),&
                  rlinform,.true.,rrhs%RvectorBlock(2),coeff_RHS_Vol_u2_2D)

!      print # of DOF
    call lsysbl_getbase_double (rVector,p_Ddata)

    call output_line ('Number of DOF: ' // sys_siL(size(p_Ddata),12) )

!    Clear the solution vector on the finest level.
    call lsysbl_clearVector(rvector)

!     For implementing boundary conditions, we use a 'filter technique with
!     discretised boundary conditions'. This means, we first have to calculate
!     a discrete version of the analytic BC, which we can implement into the
!     solution/RHS vectors using the corresponsolidmechding filter.
!     
!     Create a t_discreteBC structure where we store all discretised boundary
!     conditions.
    call bcasm_initDiscreteBC(rdiscreteBC)
    
!     We first set up the boundary conditions for the U1-velocity, then those
!     of the U2-velocity.
!     We 'know' already (from the problem definition) that we have four boundary
!     segments in the domain. Each of these, we want to use for enforcing
!     some kind of boundary condition.
!     
!     We ask the bondary routines to create a 'boundary region' - which is
!     simply a part of the boundary corresponding to a boundary segment.
!     A boundary region roughly contains the type, the min/max parameter value
!     and whether the endpoints are inside the region or not.
    do i = 1,size(rproblem%Sbc)
	call boundary_createRegion(rboundary,1,i,rboundaryRegion)
	if (rproblem%Sbc(i) .eq. 'D') then
!     The endpoint of this segment should also be Dirichlet. We set this by
!     changing the region properties in rboundaryRegion.
		rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND

!     We use this boundary region and specify that we want to have Dirichlet
!     boundary there. The following call does the following:
!     - Create Dirichlet boundary conditions on the region rboundaryRegion.
!       We specify icomponent='1' to indicate that we set up the
!       Dirichlet BC's for the first (here: one and only) component in the 
!       solution vector.
!     - Discretise the boundary condition so that the BC's can be applied
!       to matrices and vectorssolidmech
!     - Add the calculated discrete BC's to rdiscreteBC for later use.

		call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
						rboundaryRegion,rdiscreteBC,&
						getBoundaryValues_2D)

!     Now continue with defining the boundary conditions of the U2-velocity:	
		call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
						rboundaryRegion,rdiscreteBC,&
						getBoundaryValues_2D)
	else if (rproblem%Sbc(i) .eq. 'N') then
		rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND

!       For Element Q1
		if (rproblem%selement .EQ. 'Q1') then

!       We use this boundary region and specify that we want to have Neumann boundary there
!       and add this to volumetric part(the neumann boundary part on rhs) for U1
			call linf_buildVectorScalarBdr2d(rlinform,CUB_G2_1D,.false.,&
					rrhs%RvectorBlock(1),coeff_RHS_neumBdr_u1_2D,rboundaryRegion)

!       Now its done for U2
			call linf_buildVectorScalarBdr2d(rlinform,CUB_G2_1D,.false.,&
					rrhs%RvectorBlock(2),coeff_RHS_neumBdr_u2_2D,rboundaryRegion)

!       For Element Q2
		else if (rproblem%selement .EQ. 'Q2') then

!       For U1
			call linf_buildVectorScalarBdr2d(rlinform,CUB_G3_1D,.false.,&
					rrhs%RvectorBlock(1),coeff_RHS_neumBdr_u1_2D,rboundaryRegion)

!       For U2
			call linf_buildVectorScalarBdr2d(rlinform,CUB_G3_1D,.false.,&
					rrhs%RvectorBlock(2),coeff_RHS_neumBdr_u2_2D,rboundaryRegion)
		end if
	else
	print *, ' Invalid Input for Boundary Condition'
	stop
	end if
   end do

    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.!
    rmatrix%p_rdiscreteBC => rdiscreteBC
    rrhs%p_rdiscreteBC => rdiscreteBC
    rvector%p_rdiscreteBC => rdiscreteBC
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhs)
    call vecfil_discreteBCsol (rvector)
    call matfil_discreteBC (rmatrix)

    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with 
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions 
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
!     RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver.
    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
!     p_RfilterChain => RfilterChain
!     nullify(p_rpreconditioner)
!     call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)

!    For Direct solver we use following 
    call linsol_initUMFPACK4 (p_rsolverNode)
    p_rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput = 0

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! We will allow the solver to perform 200 iterations
    p_rsolverNode%nmaxIterations = 10000

    p_rsolverNode%depsRel = 1E-10_DP
   ! p_rsolverNode%depsAbs = 1E-10_DP

    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrix/)
    call linsol_setMatrices(p_rsolverNode,Rmatrices)
    
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

!    For Bending in the gmv
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_DvertexCoords)

     ! Write velocity field
    call lsyssc_getbase_double (rVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rVector%RvectorBlock(2),p_Ddata2)
  
!    we add sol. vector to coordinates to see the bending
    if (rproblem%DEFORMATION .eq. ON) then
	do i = 1,rtriangulation%NVT
		p_Dvertexcoords(1,i) = p_Dvertexcoords(1,i) + p_Ddata(i)
		p_Dvertexcoords(2,i) = p_dvertexCoords(2,i) + p_Ddata2(i)
	end do
    end if

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing. 
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/u2d_0_simple.gmv')

 

    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

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
    call lsysbl_releaseMatrix (rmatrix)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)

 end subroutine 

end module

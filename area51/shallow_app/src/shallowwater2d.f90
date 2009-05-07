!##############################################################################
!# ****************************************************************************
!# <name> Shallowwater2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a 2D Shallowwater
!# problem.
!# </purpose>
!##############################################################################

MODULE shallowwater2d

  USE fsystem
  USE genoutput
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE ucd
  USE pprocerror
  USE genoutput
  USE paramlist
  USE stdoperators

  USE shallowwater2d_routines
  
  IMPLICIT NONE

    ! This is defined in shallowwater2d_routines
	!TYPE t_array
	!	! Pointer to the double-valued matrix or vector data
    !	REAL(DP), DIMENSION(:), POINTER :: Da
	!END TYPE t_array
    
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE shallowwater2d_0
  
!<description>
  ! This is a 2d shallowwater solver.
  ! The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up basic matrices and vectors
  ! 4.) Set initial conditions
  ! 5.) Create solver structure
  ! 6.) Start a timestepping loop with defect correction
  ! 7.) 	Build the Galerkin Operator
  ! 8.) 	Apply TVD/FCT corrections
  ! 9.) 	Apply boundary-conditions
  !10.) 	Solve the linear system
  !11.) Write GMV-file
  !12.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !

	! Number of Variables (h, hu, hv)
	! (h=Waterheights, u/v=speed in x/y-direction)
	INTEGER, PARAMETER :: nvar2d = 3

    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! Scalar matrices
	! (consistent and lumped mass matrix,
    TYPE(t_matrixScalar) :: rmatrixMC, rmatrixML, rmatrixCX, rmatrixCY
    
	! Array with pointers to the datas of the diagonal blocks of P,
	! the preconditioner of the outer defect correction loop
	TYPE(t_array), DIMENSION(nvar2d)    :: rarrayP

	! Array with pointers to the datas of the different components of the
	! rhs, sol, temp, soldot, and defect
	TYPE(t_array), DIMENSION(nvar2d)    :: rarrayRhs, rarraySol, rarrayRstemp, rarraySolDot, rarrayDef

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrixBlockP
    TYPE(t_vectorBlock) :: rrhsBlock, rsolBlock, rdefBlock, rstempBlock, rsolDotBlock

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode, p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! NLMAX receives the level of the grid where we want to solve.
    INTEGER :: NLMAX
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Error of FE function to reference function
    REAL(DP) :: derror
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport)               :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata, p_Ddata1, p_Ddata2
   
    ! Current time, final time and timestepsize
    REAL(dp):: ttime, ttfinal, dt
    
	! Pointer to the entries of the matrices K, L and D
	! They are needed to apply upwind diffusion for example
	! Kedge is an array, which saves the corresponding grid-points
	! for each edge
    INTEGER, DIMENSION(:), POINTER		    :: p_Kld, p_Kcol, p_Kdiagonal
    INTEGER(I32), DIMENSION(:), POINTER     :: p_Ksep
    INTEGER(I32), DIMENSION(:,:), POINTER   :: p_Kedge
	REAL(dp), DIMENSION(:), POINTER	        :: p_CXdata, p_CYdata, p_MLdata, p_MCdata
    
    ! Size of the 2D-array Kedge
    INTEGER(I32), DIMENSION(2)              :: iSize
    
    ! Number of edges of the grid
    INTEGER                                 :: nedge
    
    ! some variables needed to apply upwind diffusion
    INTEGER							:: ieq, ild, icol, ios, icount, icount2, iedge
	INTEGER							:: ii, ij, ji, jj, i ,j, k, l, d
    REAL(DP)						:: dtmp, dtmp1, dtmp2, dtmp3, dtmp4
    
    ! To measure the time the computation took
    REAL(DP)						:: dtime1, dtime2
    
    ! Norm of defect vector
    REAL(DP)						:: dinitialDefect, dcurrentDefect
    
    ! Implicitness Parameter
	! (0=exp. euler, 1/2=crank nicolson, 1=imp. euler)
    REAL(DP)						:: theta
    
    ! Handles
    INTEGER                         :: h_Sep, h_Kedge

    ! For the defect correction loop
	! number of iterations and maximum number of iterations
    INTEGER                         :: ite, itemax

	! Loop parameter
	INTEGER							:: ivar

	! The Jacobi matrices for x and y direction
	! evaluated at the Roe-meanvalues for edge ij
	REAL(DP), DIMENSION(nvar2d,nvar2d)	:: JacobixRoeij, JacobiyRoeij

	! The antisymmetric and symmetric part of the Galerkin operator
	! and the atrifical viscosity (at edge ij)
	REAL(DP), DIMENSION(nvar2d,nvar2d)	:: Aij, Bij, Dij

	! Temporary solution vectors
	REAL(DP), DIMENSION(nvar2d) :: Qi, Qj, Qroe, Qroeij, deltaQij, deltaKi, deltaKj
	REAL(DP), DIMENSION(nvar2d) :: Udoti, Udotj
	REAL(DP), DIMENSION(nvar2d) :: deltaDi, deltaDj

	! Coeffitients for the jacobi matrices
	REAL(DP)		:: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
	
	! Pointer to vertex coordinates
	REAL(DP), DIMENSION(:,:), POINTER   :: p_dVertexCoords

	! temporary variables for the low order operator
	REAL(DP)    :: cRoe, uRoe, vRoe, scalarDissipation, lambda
	REAL(DP)    :: scalefactor, scalefactor1, scalefactor2

	! gravity constant g (usually 9.81...)
	REAL(DP)    :: gravconst

	! Pointers to the flux limiter data arrays
    REAL(DP), DIMENSION(:,:), POINTER   :: p_fld1, p_fld2

    ! Handles
    INTEGER                         :: h_fld1, h_fld2

    ! Parameter list
    TYPE(t_parlist) :: rparlist

    ! For the flux limiter
    REAL(DP), DIMENSION(nvar2d) :: deltaWij
	REAL(DP), DIMENSION(nvar2d) :: deltaFij, deltaGij
	REAL(DP), DIMENSION(nvar2d,nvar2d)	:: Rij, invRij, Eye
    REAL(DP), DIMENSION(nvar2d) :: eigenvalues
    REAL(DP) :: deltaaplus, deltaaminus, deltabplus, deltabminus
    INTEGER :: inode, upwindnode
    REAL(DP) :: deltaak , deltaWijHat

    ! Shall the Preconditioner always be updated?
    ! 0 = no, 1 = yes
    INTEGER :: alwaysupdatepreconditioner

    ! String to read info from ini file
    CHARACTER (LEN=100) :: sstring

    ! The kind of used FE
    INTEGER :: FEkind

    ! What kind of limiter to use for TVD?
    INTEGER :: limiter

    ! Number of equations
    INTEGER :: NEQ

    ! Choosing startvalues
    INTEGER :: Startvalues

    ! What low order Method to use?
    INTEGER :: Method

    ! Apply prelimiting for FCT?
    INTEGER :: prelimiting

    ! For saving the video files
    CHARACTER (LEN=10) :: sfilenumber
    real(DP) :: videotimestep, videotime
    INTEGER :: ifilenumber, makevideo

	! Stopping criteria for the nonlinear and linear solver
	REAL(DP) :: nonlinabsdef, nonlinreldef, nonlinsolup
	REAL(DP) :: linabsdef, linreldef
	
	! Shall boundary conditions be applied in a corner?
	INTEGER :: boundarycorner


! OK, LET'S START

    ! Time measurement
    CALL CPU_TIME(dtime1)


    
! Initialise some values
	! Read parameter file
    CALL parlst_init(rparlist)
    CALL parlst_readFromFile(rparlist,'./dat/1.dat')

	! We want to solve our problem on level... Default=1
    CALL parlst_getvalue_int(rparlist, 'TRIANGULATION', 'NLMAX', nlmax, 1)

    ! And with timestepsize
    CALL parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'dt', dt)
    
    ! To the final time
    CALL parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'ttfinal', ttfinal)

    ! Set implicitness parameter. Default=Crank-Nicolson
    CALL parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'theta', theta, 0.5_DP)
    
    ! Always update Preconditioner? Default=Yes
    CALL parlst_getvalue_int(rparlist, 'SOLVER', &
                            'alwaysupdatepreconditioner', alwaysupdatepreconditioner, 1)
    
    ! What kind of limiter to use? Default=Van Leer
    CALL parlst_getvalue_int(rparlist, 'METHOD', 'limiter', limiter, 2)
    
    ! Gravitational constant. Default=9.81
    CALL parlst_getvalue_double(rparlist, 'PROBLEM', 'gravconst', gravconst, 9.81_DP)
    
    ! Choosing the Startvalues
    CALL parlst_getvalue_int(rparlist, 'PROBLEM', 'Startvalues', Startvalues)
    
    ! Choosing the Method to use. Default=TVD
    CALL parlst_getvalue_int(rparlist, 'METHOD', 'Method', Method, 3)

    ! Apply prelimiting for FCT? Default=Yes.
    CALL parlst_getvalue_int(rparlist, 'METHOD', 'prelimiting', prelimiting, 1)
    
    ! Make gmv snapshots for video? Default=No.
    CALL parlst_getvalue_int(rparlist, 'VIDEO', 'makevideo', makevideo, 0)
    
    ! Make gmv snapshot every ... seconds (must be n*dt)
    CALL parlst_getvalue_double(rparlist, 'VIDEO', 'videotimestep', videotimestep, 50000.0_DP)

	! Maximum number of iterations for the nonlinear solver
	CALL parlst_getvalue_int(rparlist, 'SOLVER', 'itemax', itemax, 20)
	
	! Shall boundary conditions be applied in a corner?
	CALL parlst_getvalue_int(rparlist, 'SOLVER', 'boundarycorner', boundarycorner, 0)

	! Stopping criterium for the nonlinear solver
	! (this is where you determine how accurate your solution is - this should be small)
	! Absolute value of the norm of the defect
	CALL parlst_getvalue_double(rparlist, 'SOLVER', 'nonlinabsdef', nonlinabsdef, 1e-11_DP)
	! Relative value of the norm of the defect
	CALL parlst_getvalue_double(rparlist, 'SOLVER', 'nonlinreldef', nonlinreldef, 1e-10_DP)
	! Norm of solutionupdate
	CALL parlst_getvalue_double(rparlist, 'SOLVER', 'nonlinsolup', nonlinsolup, 1e-12_DP)

	! Stopping criterium for linear solver
	! (this is only for the inner defect correction loop - doesn't have to be too small)
	! Absolute value of the norm of the defect
	CALL parlst_getvalue_double(rparlist, 'SOLVER', 'linabsdef', linabsdef, 1e-7_DP)
	! Relative value of the norm of the defect
	CALL parlst_getvalue_double(rparlist, 'SOLVER', 'linreldef', linreldef, 1e-4_DP)

   
        
! Create the discretisation of the domain.
!	1. Read in the PARametrisation of the boundary
!	2. Read in the TRIangulation
!	3. Refine it
!	4. Create the actual mesh from this raw data
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL parlst_getvalue_string (rparlist, 'TRIANGULATION', &
                                            'prmname', sstring)
    CALL boundary_read_prm(rboundary, sstring)

    ! Now read in the basic triangulation.
    CALL parlst_getvalue_string (rparlist, 'TRIANGULATION', &
                                            'triname', sstring)
    CALL tria_readTriFile2D (rtriangulation, sstring, rboundary, .TRUE.)

    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (NLMAX-1, rtriangulation, rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation, rboundary)
    
    
    
! Initialise the discretisation

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. For the shallow water problem we need three blocks
	CALL spdiscr_initBlockDiscr (rdiscretisation, nvar2d, rtriangulation)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:

	! Initalise the first block of the discretisation
	! Get the type of FE from the parameter file
	CALL parlst_getvalue_int(rparlist, 'TRIANGULATION', 'FEkind', FEkind)
	IF (FEkind == 0) THEN
	    FEkind = EL_E001
	ELSE IF (FEkind == 1) THEN
	    FEkind = EL_E011
	END IF
	CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   FEkind,SPDISC_CUB_AUTOMATIC, rtriangulation, &
								   rboundary)

	! Now copy this initialised block into the other ones
	! But only create a derived copy, which shares the handles of the original one
	DO ivar = 2, nvar2d
		CALL spdiscr_duplicateDiscrSc (rdiscretisation%Rspatialdiscr(1), &
									  rdiscretisation%Rspatialdiscr(ivar), &
									  .TRUE.)
	END DO
    
    
    
! Now create the system matrices

    ! First create matrix MC - the consistent mass matrix

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of some matrices we will later need to build
	! the preconditioner of the outer defect correction loop.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our solution components
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixMC)
        
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (Psi_j, Phi_i) for the
    ! scalar matrix in 2D.
    
	CALL stdop_assembleSimpleMatrix(rmatrixMC, DER_FUNC, DER_FUNC)


	! Next we create the Matrices CX and CY 
	! They represent the diskretisation of the nabla operator
    
    ! We could do this as we did with the mass matrix MC, but as CX and MC
    ! are of the same structure we can as well create a shared matrix,
    ! that shares the handles for the structure (Kld, Kcol, Kdiagonal)
    ! with the Matrix MC and only creates a new haldle for the entries (Data)
    ! So instead of using
    ! CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
    !                                    LSYSSC_MATRIX9,rmatrixCX)
    ! we use
    CALL lsyssc_duplicateMatrix(rmatrixMC, rmatrixCX,&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

	CALL stdop_assembleSimpleMatrix(rmatrixCX, DER_DERIV_X, DER_FUNC)

	! Now we do the same for CY
    CALL lsyssc_duplicateMatrix(rmatrixMC, rmatrixCY,&
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

	CALL stdop_assembleSimpleMatrix(rmatrixCY, DER_DERIV_Y, DER_FUNC)


   
 	! Now get all the needed pointers to the datas of the matrices
	CALL lsyssc_getbase_double(rmatrixCX,p_CXdata)
 	CALL lsyssc_getbase_double(rmatrixCY,p_CYdata)
 	CALL lsyssc_getbase_Kcol(rmatrixMC,p_Kcol)
 	CALL lsyssc_getbase_Kld(rmatrixMC,p_Kld)
 	CALL lsyssc_getbase_Kdiagonal(rmatrixMC,p_Kdiagonal)


    ! Number of equations
    NEQ = rmatrixMC%NEQ


    ! Create Vector Ksep
	! It is used to to create the array Kedge
    CALL storage_new1D ('tvd', 'Ksep', rmatrixMC%NEQ, ST_INT, &
                     h_Sep, ST_NEWBLOCK_ZERO)
    
    ! and get Pointer to Ksep
	CALL storage_getbase_int (h_Sep, p_Ksep)

	! Calculate number of edges
    nedge = (rMatrixMC%NA-rMatrixMC%neq)/2
	
	! Create Array Kedge
	iSize(1)=4
	iSize(2)=nedge
	 CALL storage_new2D ('tvd', 'Kedge', iSize, ST_INT, h_Kedge, &
                            ST_NEWBLOCK_ZERO)
    
    ! and get Pointer to Kedge
    CALL storage_getbase_int2D (h_Kedge, p_Kedge)


	! Now calculate the entries of Kedge
	iedge = 0
    p_Ksep(1:rmatrixMC%NEQ) = p_Kld(1:rmatrixMC%NEQ)
	! Walk throug all lines
	DO ieq = 1, rmatrixMC%neq-1
		DO ij = p_Kdiagonal(ieq)+1, p_Kld(ieq+1)-1
			! ij is the current position
			icol = p_Kcol(ij)
			! It is located in line ieq and column icol
			! ji is supposed to point to the transposed entrie
			ji = p_Ksep(icol)
			p_Ksep(icol)=p_Ksep(icol)+1
			! save the edges and nodes
			iedge=iedge+1
			p_Kedge(1,iedge)=ieq   ! node i for this edge
			p_Kedge(2,iedge)=icol  ! node j for this edge
			p_Kedge(3,iedge)=ij    ! ij
			p_Kedge(4,iedge)=ji    ! ji
		END DO
	END DO
	
	
	! For the application of the TVD correction we will need the
	! temporary arrays fld1+2

	! Create Array fld1
	iSize(1)=nvar2d*6
	iSize(2)=rMatrixMC%neq
	 CALL storage_new2D ('Limiter', 'fld1', iSize, ST_DOUBLE, h_fld1, &
                            ST_NEWBLOCK_ZERO)
    
    ! and get Pointer to fld1
    CALL storage_getbase_double2D (h_fld1, p_fld1)
	
	! Create Array fld2
	iSize(1)=nvar2d*4
	iSize(2)=nedge
	 CALL storage_new2D ('Limiter', 'fld2', iSize, ST_DOUBLE, h_fld2, &
                            ST_NEWBLOCK_ZERO)
    
    ! and get Pointer to fld2
    CALL storage_getbase_double2D (h_fld2, p_fld2)
	


	! Now create the preconditioner block matrix P
	! First create an empty block matrix structur with nvar2d x nvar2d blocks
	CALL lsysbl_createEmptyMatrix (rmatrixBlockP, nvar2d)

	! Next create the diagonal blocks of P as empty matrices, using the 
	! matrix structur of the matrix MC
	! We will only need the diagonal blocks, as we employ a block jacobi
	! method here
	DO ivar = 1, nvar2d
		CALL lsyssc_duplicateMatrix (rmatrixMC, &
									 rmatrixBlockP%Rmatrixblock(ivar,ivar), &
									 LSYSSC_DUP_SHARE, &
									 LSYSSC_DUP_EMPTY)
	END DO
	
	
    ! Create ML (lumped Mass Matrix)
	! Copy MC into ML
	CALL lsyssc_duplicateMatrix (rMatrixMC,rMatrixML,&
                                     LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    
    ! Apply Mass Lumping
    CALL lsyssc_lumpMatrixScalar (rmatrixML,LSYSSC_LUMP_DIAG,.TRUE.)
    
    ! Get Pointer to MC and ML Datas
    CALL lsyssc_getbase_double(rmatrixMC,p_MCdata)
    CALL lsyssc_getbase_double(rmatrixML,p_MLdata)
  

    ! Initialise Dij, the artifical viscosity for the artificial
	! diffusion operator D
    Dij = 0.0_DP
    
    ! Create the identity matrix
    Eye = 0.0_DP
    DO ivar = 1, nvar2d
        Eye(ivar,ivar) = 1.0_DP
    END DO

                 
! Now create all needed scalar and block vectors

	! First create the block vectors
	! rtempBlock, rrhsBlock, rsolBlock, rdefBlock and rstempBlock, rsolDotBlock
	CALL lsysbl_createVecBlockByDiscr (rDiscretisation,rrhsBlock,.TRUE.,&
                                           ST_DOUBLE)
	CALL lsysbl_createVecBlockByDiscr (rDiscretisation,rsolBlock,.TRUE.,&
                                           ST_DOUBLE)
	CALL lsysbl_createVecBlockByDiscr (rDiscretisation,rdefBlock,.TRUE.,&
                                           ST_DOUBLE)
	CALL lsysbl_createVecBlockByDiscr (rDiscretisation,rstempBlock,.TRUE.,&
                                           ST_DOUBLE)
    CALL lsysbl_createVecBlockByDiscr (rDiscretisation,rsolDotBlock,.TRUE.,&
                                           ST_DOUBLE)
                                                  
                                           

	! Get pointers to the components of P, rhs, sol, rstemp, rsolDotBlock
	DO ivar = 1, nvar2d
    	CALL lsyssc_getbase_double(rmatrixBlockP%RmatrixBlock(ivar,ivar), rarrayP(ivar)%Da)	! P
		CALL lsyssc_getbase_double(rsolBlock%RvectorBlock(ivar), rarraySol(ivar)%Da)	! Sol
		CALL lsyssc_getbase_double(rrhsBlock%RvectorBlock(ivar), rarrayrhs(ivar)%Da)	! Rhs
		CALL lsyssc_getbase_double(rstempBlock%RvectorBlock(ivar), rarrayrstemp(ivar)%Da)	! Rstemp
		CALL lsyssc_getbase_double(rsolDotBlock%RvectorBlock(ivar), rarraysolDot(ivar)%Da)	! rsolDotBlock
		CALL lsyssc_getbase_double(rdefBlock%RvectorBlock(ivar), rarraydef(ivar)%Da)  ! Defect
	END DO

    

    
! Initialise the linear solver

       
    ! Create a BiCGStab-solver without preconditioner.and without a filterchain
	! as we've got our own routine to implement the boundary conditions
	NULLIFY(p_rpreconditioner)
    NULLIFY(p_RfilterChain)
    CALL linsol_initBiCGStab (p_rsolverNode)

    ! The linear solver stops, when this relative or absolut norm of
    ! the residual is reached.
    p_rsolverNode%depsRel = linreldef
    p_rsolverNode%depsAbs = linabsdef
	! we don't have to solve the system to precise, as it is only used for
	! the inner loop of the defect correction

    ! Set the output level of the solver 
    ! 0: No output, just warnings
    ! 2: Some output
    p_rsolverNode%ioutputLevel = 0
       
    ! Attach the preconditioner of the outer defect correction loop
	! to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrixBlockP/)
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices)
       
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initStructure (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! *************************************** !
!      CALL linsol_initBiCGStab (p_rsolverNode1)
!      p_rsolverNode1%depsRel = 1E-4_DP
!      p_rsolverNode1%depsAbs = 1E-7_DP
!      p_rsolverNode1%ioutputLevel = 0
!      Rmatrices1 = (/rmatrixBlockP1/)
!      CALL linsol_setMatrices(p_RsolverNode1,Rmatrices1)
!      CALL linsol_initStructure (p_rsolverNode1, ierror)
!      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
!      CALL linsol_initData (p_rsolverNode1, ierror)
!      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
!      
!      CALL linsol_initBiCGStab (p_rsolverNode2)
!      p_rsolverNode2%depsRel = 1E-4_DP
!      p_rsolverNode2%depsAbs = 1E-7_DP
!      p_rsolverNode2%ioutputLevel = 0
!      Rmatrices2 = (/rmatrixBlockP2/)
!      CALL linsol_setMatrices(p_RsolverNode2,Rmatrices2)
!      CALL linsol_initStructure (p_rsolverNode2, ierror)
!      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
!      CALL linsol_initData (p_rsolverNode2, ierror)
!      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
!      
!      CALL linsol_initBiCGStab (p_rsolverNode3)
!      p_rsolverNode3%depsRel = 1E-4_DP
!      p_rsolverNode3%depsAbs = 1E-7_DP
!      p_rsolverNode3%ioutputLevel = 0
!      Rmatrices1 = (/rmatrixBlockP3/)
!      CALL linsol_setMatrices(p_RsolverNode3,Rmatrices1)
!      CALL linsol_initStructure (p_rsolverNode3, ierror)
!      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
!      CALL linsol_initData (p_rsolverNode3, ierror)
!      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! *************************************** !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Bevor we start the timestepping we should set the initial values
    ! get Pointer to the vertex coordinates
	CALL storage_getbase_double2d(rtriangulation%h_DvertexCoords, p_dVertexCoords)

    
    IF (Startvalues==1) THEN

    rarraySol(1)%Da(1:rmatrixMC%neq) = 0.1_DP
    rarraySol(2)%Da(1:rmatrixMC%neq) = 0_DP
    rarraySol(3)%Da(1:rmatrixMC%neq) = 0_DP
     rarraySol(1)%Da(1:rmatrixMC%neq) = 1_DP
    rarraySol(2)%Da(1:rmatrixMC%neq) = -5_DP
    rarraySol(3)%Da(1:rmatrixMC%neq) = 0_DP
    DO i = 1, rmatrixMC%neq
        ! calculate squared distance to (0.0,0.0)
        dtmp1 = ((p_DvertexCoords(1,i)-0.0_DP)**2.0_DP+(p_DvertexCoords(2,i)-0.0_DP)**2.0_DP)
        IF (dtmp1 < 2.5_DP**2.0_DP) THEN
            !rarraySol(1)%Da(i) = 1.0_DP + 0.4_DP*2.71_DP**(-5.0_DP*dtmp1)
    !        rarraySol(1)%Da(i) = 10.5_DP
        END IF
        IF (p_DvertexCoords(1,i)>0) then
         rarraySol(1)%Da(i) = 1_DP
    rarraySol(2)%Da(i) = 5.0_DP
    
        END IF
      
    END DO

    ELSEIF (Startvalues==2) THEN

    rarraySol(1)%Da(1:rmatrixMC%neq) = 1.0_DP
    DO i = 1, rmatrixMC%neq
        ! calculate squared distance to (0,0)
        dtmp1 = ((p_DvertexCoords(1,i))**2.0_DP+(p_DvertexCoords(2,i))**2.0_DP)
        IF (dtmp1 < 25) THEN
            
            rarraySol(1)%Da(i) = 10.0_DP
        END IF
    END DO
    
    ELSEIF (Startvalues==3) THEN

    rarraySol(1)%Da(1:rmatrixMC%neq) = 1.0_DP
    DO i = 1, rmatrixMC%neq
        ! calculate squared distance to (0,0)
        dtmp1 = ((p_DvertexCoords(1,i))**2.0_DP+(p_DvertexCoords(2,i))**2.0_DP)
        IF (dtmp1 < 81) THEN
            
            rarraySol(1)%Da(i) = 3.0_DP
        END IF
    END DO
    
    ELSEIF (Startvalues==4) THEN

    rarraySol(1)%Da(1:rmatrixMC%neq) = 1.0_DP
    DO i = 1, rmatrixMC%neq
        ! calculate squared distance to (0,0)
        dtmp1 = ((p_DvertexCoords(1,i)-8)**2.0_DP+(p_DvertexCoords(2,i))**2.0_DP)
        IF (dtmp1 < 1) THEN
            
            rarraySol(1)%Da(i) = 3.0_DP
            rarraySol(2)%Da(i) = 3.0_DP*3.0_DP
        END IF
    END DO

ELSEIF (Startvalues==5) THEN

    rarraySol(1)%Da(1:rmatrixMC%neq) = 1.0_DP
    WHERE(p_DvertexCoords(1,:)<0)
        
            rarraySol(1)%Da(:) = 3.0_DP
           
        END WHERE

ELSEIF (Startvalues==6) THEN

    rarraySol(1)%Da(1:rmatrixMC%neq) = 1.0_DP
    DO i = 1, rmatrixMC%neq
        ! calculate squared distance to (0,0)
        dtmp1 = ((p_DvertexCoords(1,i)-7)**2.0_DP+(p_DvertexCoords(2,i))**2.0_DP)
        IF (dtmp1 < 4) THEN
            
            rarraySol(1)%Da(i) = 10.0_DP
            
        END IF
    END DO

    ELSEIF (Startvalues==7) THEN   ! Dam Break
        ! Set initial values
        rarraySol(1)%Da(1:rmatrixMC%neq) = 5.0_DP
        WHERE(p_DvertexCoords(1,:)<101)
                rarraySol(1)%Da(:) = 10.0_DP
        END WHERE
        ! Set boundary conditions 
        
        
    ELSEIF (Startvalues==8) THEN   ! Balken
        ! Set initial values
        rarraySol(1)%Da(1:rmatrixMC%neq) = 5.0_DP
        WHERE(p_DvertexCoords(1,:)<12)
                rarraySol(1)%Da(:) = 10.0_DP
        END WHERE
        ! Set boundary conditions 
        

    ELSEIF (Startvalues==9) THEN   ! Zylinder
        ! Set initial values
        rarraySol(1)%Da(1:rmatrixMC%neq) = 5.0_DP
	rarraySol(2)%Da(1:rmatrixMC%neq) = 0.0_DP
	rarraySol(3)%Da(1:rmatrixMC%neq) = 0.0_DP
        WHERE(p_DvertexCoords(1,:)>40)
                rarraySol(1)%Da(:) = 10.0_DP
        END WHERE
        ! Set boundary conditions 
        
        
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    ! Write first video file (the initial conditions)
	!  if we want to make a video
    IF (makevideo == 1) THEN
        ! Start UCD export to GMV file:
        CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/video1.gmv')
    
        CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
        CALL ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
	    CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
        CALL ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
	    CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
        CALL ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
	    CALL ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
	
        ! Write the file to disc
        CALL ucd_write (rexport)
        CALL ucd_release (rexport)

        videotime = videotimestep
        ! Number of first video file
        ifilenumber = 1
    END IF




    
! Begin timestepping

    ttime = 0.0_DP

    timestepping: DO
        
       ! Compute solution from time step t^n to time step t^{n+1}
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*) 'TIME STEP:', ttime
       WRITE(*,*)
       
       ! Print what method is used
       IF (Method==0) THEN
        WRITE(*,*), 'Calculating the high order (group FEM Galerkin) solution'
       ELSE IF (Method==1) THEN
        WRITE(*,*), 'Calculating the low order (scalar dissipation) solution'
       ELSE IF (Method==2) THEN
        WRITE(*,*), 'Calculating the low order (tensorial dissipation) solution'
       ELSE IF (Method==3) THEN
        WRITE(*,*), 'Calculating the TVD solution'
        IF (limiter==1) THEN
         WRITE(*,*), 'Limiter: Minmod'
        ELSE IF (limiter==2) THEN
         WRITE(*,*), 'Limiter: VanLeer'
        ELSE IF (limiter==3) THEN
         WRITE(*,*), 'Limiter: MC'
        ELSE IF (limiter==4) THEN
         WRITE(*,*), 'Limiter: Superbee'
        END IF
       ELSE IF (Method==4) THEN
        WRITE(*,*), 'FCT: Calculating the low order (scalar dissipation) predictor'
       ELSE IF (Method==5) THEN
        WRITE(*,*), 'FCT: Calculating the low order (tensorial dissipation) predictor'
       END IF

      
      ! Assemble the preconditioner of the outer defect correction loop
      CALL BuildShallowWaterPreconditioner (rmatrixBlockP, &
	                rarrayP, rarraySol, p_CXdata, p_CYdata, &
	                p_MLdata, p_Kdiagonal, p_kedge, &
	                NEQ, nedge, theta, dt, gravconst)
	  
	  ! Now assemble RHS
      CALL BuildShallowWaterRHS (&
	                rarrayRhs, rarraySol, rrhsBlock, rsolBlock, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
	                h_fld1, p_fld1, p_fld2, &
	                p_Kdiagonal, p_Kedge, NEQ, nedge, &
	                theta, dt, gravconst, Method, limiter)
      

  ! Here starts the defect correction loop
       def_corr: DO ite = 1, itemax
       
            IF ((alwaysupdatepreconditioner==1).AND.(ite.NE.1)) THEN
            ! Reassemble the Preconditioner
                CALL BuildShallowWaterPreconditioner (rmatrixBlockP, &
	                rarrayP, rarraySol, p_CXdata, p_CYdata, &
	                p_MLdata, p_Kdiagonal, p_kedge, &
	                NEQ, nedge, theta, dt, gravconst)
            END IF
       
            ! Assemble defect
            CALL BuildShallowWaterDefect (&
	                rdefBlock, rstempBlock, rrhsBlock, rsolBlock, &
	                rarrayRhs, rarraySol, rarrayRstemp, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
	                h_fld1, p_fld1, p_fld2, &
	                p_Kdiagonal, p_Kedge, NEQ, nedge, &
	                theta, dt, gravconst, Method, limiter)
	                
	        ! Take care of Boundary Conditions
	        CALL ImplementShallowWaterBCs (&
	                rboundary, rtriangulation, &
	                rarrayP, rarraySol, rarrayDef, &
	                p_Kdiagonal, p_Kld, &
	                gravconst, boundarycorner)

       
            ! Compute norm of defect
			dcurrentDefect = lsysbl_vectorNorm (rdefBlock,LINALG_NORML2)

			! Save the initial defect
            IF (ite == 1) THEN
                dinitialDefect = dcurrentDefect
			END IF

            ! Print norm of the current defect
            PRINT *, "DEFCOR", ite, dcurrentDefect

            ! If norm of defect is small enough then leave the
            ! defect correction loop
            IF (dcurrentDefect < nonlinabsdef .OR.&
                dcurrentDefect < nonlinreldef*dinitialDefect) EXIT def_corr
       
            ! Finally solve the system.
            ! We have calculated the defect and want to calculate a
            ! solution update that we can add to the solution.
            ! So we use
            ! CALL linsol_precondDefect 
            ! If we wanted to solve Au=b instead with
            ! b being the real RHS and u being the real solution vector,
            ! we would use 
            ! CALL linsol_solveAdaptively (p_rsolverNode,rsolBlock,rrhsBlock,rtempBlock).
            CALL linsol_precondDefect (p_rsolverNode,rdefBlock)		
            ! Now we have the solution update in rdef.
       
            ! Add the solution update to the solution: rsol=rsol+rdef
            CALL lsysbl_vectorLinearComb (rdefBlock,rsolBlock,1.0_dp,1.0_dp)

			! If norm of solution update is small enough then leave the
            ! defect correction loop
			dcurrentDefect = lsysbl_vectorNorm (rdefBlock,LINALG_NORML2)
            IF (dcurrentDefect < nonlinsolup) EXIT def_corr
       
       END DO def_corr ! This is the end of the defect correction loop

	! Test if maximum number of iterations was reached and output warning
	IF (ite == itemax+1) THEN
		WRITE(*,*) ''
		WRITE(*,*) '***********************************************************************'
		WRITE(*,*) '* WARNING: Stopping criteria of the nonlinear solver were not reached *'
		WRITE(*,*) '***********************************************************************'
		WRITE(*,*) ''
	END IF
       
    
    ! If we are using TVD, then we are finished here (at least for this timestep)
	! But if we are using the linearised FCT, we've got only the low order predictor
	! and need to add explicit linearised limited antidiffusion
    IF ((Method==4).OR.(Method==5)) THEN
        WRITE(*,*) 'FCT: Adding limited antidiffusion'
    
        ! This routine adds limited antidiffusion according to the linearised
        ! FCT limiter
        CALL FctShallowWaterAddLimitedAntidiffusion(&
	                rarraySol, rarraySolDot, rarrayRhs,&
	                rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
                    h_fld1, p_fld1, p_fld2, &
                    p_Kdiagonal, p_Kedge, NEQ, nedge, &
                    gravconst, dt, Method, prelimiting)
    END IF

	! That's it. RvectorBlock now contains our solution at the current time

    ! write gmvfiles for video (if needed)
    IF ((makevideo==1).AND.(ABS(videotime-ttime)<(dt*0.2_DP))) THEN
    
        WRITE(*,*)
        WRITE(*,*) 'Writing Videofile'
        WRITE(*,*)
    
        ifilenumber = ifilenumber + 1
    
        WRITE(sfilenumber,'(i0)') ifilenumber
        videotime = videotime + videotimestep
    
        CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/video' // TRIM(sfilenumber) // '.gmv')
    
        CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
        CALL ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
	    CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
        CALL ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
	    CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
        CALL ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
	    CALL ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
    
        ! Write the file to disc, that's it.
        CALL ucd_write (rexport)
        CALL ucd_release (rexport)
    END IF

	! Go on to the next time step
    ttime = ttime + dt

    ! Leave the time stepping loop if final time is reached
    IF (ttime .GE. ttfinal-0.001_DP*dt) EXIT timestepping

END DO timestepping


! Release the memory of the solver now

	! Release solver data and structure
	CALL linsol_doneData (p_rsolverNode)
	CALL linsol_doneStructure (p_rsolverNode)
       
	! Release the solver node and all subnodes attached to it (if at all):
	CALL linsol_releaseSolver (p_rsolverNode)

    
    
! Write the calculated solution to a gmv file
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    WRITE(*,*)
    WRITE(*,*) 'Writing Solution at final time',ttime,'to File'
    WRITE(*,*)
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/u2d.gmv')
    
    ! We could write Soldot, the approximation of the time derivative
	! of the solution, which was used while applying the linearised
	! FCT-Limiter, to the file, too.
    !CALL ucd_addVariableVertexBased (rexport,'solDot_uh',UCD_VAR_STANDARD, rarraySolDot(1)%Da)
    
	! Get pointers to the solution data
    CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
	CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
    CALL ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
	CALL lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
    CALL ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
	CALL ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    
! Now release all memory
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    
    ! Release the block matrix/vectors
	CALL lsysbl_releaseVector (rstempBlock)
    CALL lsysbl_releaseVector (rsolBlock)
    CALL lsysbl_releaseVector (rsolDotBlock)
    CALL lsysbl_releaseVector (rdefBlock)
    CALL lsysbl_releaseVector (rrhsBlock)
	CALL lsysbl_releaseMatrix (rmatrixBlockP)
    
    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    CALL storage_free (h_Sep)
    CALL storage_free (h_Kedge)
    CALL storage_free (h_fld1)
    CALL storage_free (h_fld2)
	CALL lsyssc_releaseMatrix (rmatrixCX)
    CALL lsyssc_releaseMatrix (rmatrixCY)
    CALL lsyssc_releaseMatrix (rmatrixMC)
    CALL lsyssc_releaseMatrix (rmatrixML)

	! Release the memory associated to the parameter file
	CALL parlst_done (rparlist)
    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
    ! Time measurement
    CALL CPU_TIME(dtime2)
    WRITE(*,*) "Computational time used:", dtime2-dtime1
    WRITE(*,*)
    
  END SUBROUTINE

END MODULE

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

module shallowwater2d

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use linearalgebra
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use pprocerror
  use genoutput
  use paramlist
  use stdoperators
  use fparser
  use element
  use derivatives
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use filtersupport
  use linearsolver
  use feevaluation

  use shallowwater2d_routines

  implicit none

  ! This is defined in shallowwater2d_routines
  !TYPE t_array
  !	! Pointer to the double-valued matrix or vector data
  !	REAL(DP), DIMENSION(:), POINTER :: Da
  !END TYPE t_array

contains

  ! ***************************************************************************

  !<subroutine>

  subroutine shallowwater2d_0

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
    integer, parameter :: nvar2d = 3

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

    ! Scalar matrices
    ! (consistent and lumped mass matrix,
    type(t_matrixScalar) :: rmatrixMC, rmatrixML, rmatrixCX, rmatrixCY, rmatrixBX, rmatrixBY

    ! Array with pointers to the datas of the diagonal blocks of P,
    ! the preconditioner of the outer defect correction loop
    type(t_array), dimension(nvar2d)    :: rarrayP

    ! Array with pointers to the datas of the different components of the
    ! rhs, sol, temp, soldot, and defect
    type(t_array), dimension(nvar2d)    :: rarrayRhs, rarraySol, rarrayRstemp, rarraySolDot, rarrayDef, rarraySource

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlockP
    type(t_vectorBlock), target :: rrhsBlock, rsolBlock, rdefBlock, rstempBlock, rsolDotBlock, rsourceBlock, rBottomBlock
    
    ! Scalar vector to describe the bottom profile
    type(t_vectorScalar), target :: rvectorbottom, rvecBSX, rvecBSY

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! The collection to give data to the callback routine
    type(t_collection) :: Rcollection

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! NLMAX receives the level of the grid where we want to solve.
    integer :: NLMAX

    ! Error indicator during initialisation of the solver
    integer :: ierror

    ! Error of FE function to reference function
    real(DP) :: derror

    ! Output block for UCD output to GMV file
    type(t_ucdExport)               :: rexport
    real(DP), dimension(:), pointer :: p_Ddata, p_Ddata1, p_Ddata2

    ! Current time, final time and timestepsize
    real(dp):: ttime, ttfinal, dt

    ! Pointer to the entries of the matrices K, L and D
    ! They are needed to apply upwind diffusion for example
    ! Kedge is an array, which saves the corresponding grid-points
    ! for each edge
    integer, dimension(:), pointer     :: p_Kld, p_Kcol, p_Kdiagonal
    integer, dimension(:), pointer     :: p_Ksep
    integer, dimension(:,:), pointer   :: p_Kedge
    real(dp), dimension(:), pointer	   :: p_CXdata, p_CYdata, p_MLdata, p_MCdata
    real(dp), dimension(:), pointer    :: p_BXdata, p_BYdata, p_Bottom, p_BSXdata, p_BSYdata
    real(dp), dimension(:), pointer    :: p_bottomvector

    ! Size of the 2D-array Kedge
    integer, dimension(2)              :: iSize

    ! Number of edges of the grid
    integer                            :: nedge

    ! some variables needed to apply upwind diffusion
    integer							:: ieq, ild, icol, ios, icount, icount2, iedge
    integer							:: ii, ij, ji, jj, i ,j, k, l, d
    real(DP)						:: dtmp, dtmp1, dtmp2, dtmp3, dtmp4

    ! To measure the time the computation took
    real(DP)						:: dtime1, dtime2

    ! Norm of defect vector
    real(DP)						:: dinitialDefect, dcurrentDefect

    ! Implicitness Parameter
    ! (0=exp. euler, 1/2=crank nicolson, 1=imp. euler)
    real(DP)						:: theta

    ! Handles
    integer                         :: h_Sep, h_Kedge, h_Bottom

    ! For the defect correction loop
    ! number of iterations and maximum number of iterations
    integer                         :: ite, itemax

    ! Loop parameter
    integer							:: ivar

    ! The Jacobi matrices for x and y direction
    ! evaluated at the Roe-meanvalues for edge ij
    real(DP), dimension(nvar2d,nvar2d)	:: JacobixRoeij, JacobiyRoeij

    ! The antisymmetric and symmetric part of the Galerkin operator
    ! and the atrifical viscosity (at edge ij)
    real(DP), dimension(nvar2d,nvar2d)	:: Aij, Bij, Dij

    ! Temporary solution vectors
    real(DP), dimension(nvar2d) :: Qi, Qj, Qroe, Qroeij, deltaQij, deltaKi, deltaKj
    real(DP), dimension(nvar2d) :: Udoti, Udotj
    real(DP), dimension(nvar2d) :: deltaDi, deltaDj

    ! Coeffitients for the jacobi matrices
    real(DP)		:: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB

    ! Pointer to vertex coordinates
    real(DP), dimension(:,:), pointer   :: p_dVertexCoords

    ! temporary variables for the low order operator
    real(DP)    :: cRoe, uRoe, vRoe, scalarDissipation, lambda
    real(DP)    :: scalefactor, scalefactor1, scalefactor2

    ! gravity constant g (usually 9.81...)
    real(DP)    :: gravconst

    ! Pointers to the flux limiter data arrays
    real(DP), dimension(:,:), pointer   :: p_fld1, p_fld2

    ! Handles
    integer                         :: h_fld1, h_fld2

    ! Parameter list
    type(t_parlist) :: rparlist

    ! For the flux limiter
    real(DP), dimension(nvar2d) :: deltaWij
    real(DP), dimension(nvar2d) :: deltaFij, deltaGij
    real(DP), dimension(nvar2d,nvar2d)	:: Rij, invRij, Eye
    real(DP), dimension(nvar2d) :: eigenvalues
    real(DP) :: deltaaplus, deltaaminus, deltabplus, deltabminus
    integer :: inode, upwindnode
    real(DP) :: deltaak , deltaWijHat

    ! Shall the Preconditioner always be updated?
    ! 0 = no, 1 = yes
    integer :: alwaysupdatepreconditioner

    ! String to read info from ini file and the string of the function of the bottom profile
    character (LEN=200) :: sstring, sbottomstring
    
    ! Name of output file
    character (LEN=200) :: ofile

    ! The kind of used FE
    integer :: FEkind
    integer(I32) :: celement

    ! What kind of limiter to use for TVD?
    integer :: limiter

    ! Number of equations
    integer :: NEQ

    ! What low order Method to use?
    integer :: Method

    ! Apply prelimiting for FCT?
    integer :: prelimiting

    ! For saving the video files
    character (LEN=10) :: sfilenumber
    real(DP) :: videotimestep, videotime
    integer :: ifilenumber, makevideo

    ! Stopping criteria for the nonlinear and linear solver
    real(DP) :: nonlinabsdef, nonlinreldef, nonlinsolup
    real(DP) :: linabsdef, linreldef

    ! Shall boundary conditions be applied in a corner?
    integer :: boundarycorner
    
    ! Shall reflecting boundary conditions be applied?
    integer :: reflectingbcs

    ! Syncronisation method for FCT limiting
    integer :: syncromethod
    
    ! Is the initial value given by absolute or relative value
    integer :: absrel
    
    ! Add bottom profile to output of h?
    integer :: addbottomtoout
    
    ! Use the bottom-sourceterm?
    integer :: bottomterm
    
    ! Command line and name of the paramater file
    character(LEN=SYS_STRLEN) :: cbuffer, sparameterfileName

    ! Function parser
    type(t_fparser) :: rfparser
    character(LEN=*), dimension(2), parameter ::&
         cvariables = (/ (/'x'/), (/'y'/) /)
         
    ! Some auxiliary values for the bottom source term
    real(DP), dimension(2):: b_x, b_y
    real(DP), dimension(2,2):: Dpoints
    
    ! Counter for number of warnings
    integer :: numwarnings


    ! OK, LET'S START

    ! Time measurement
    call cpu_time(dtime1)
    
    ! Reset the numer of warnings
    numwarnings = 0
    
    
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
  

    ! Initialise some values
    ! Read parameter file
    call parlst_init(rparlist)
    call parlst_readFromFile(rparlist,sparameterfileName)

    ! We want to solve our problem on level... Default=1
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'NLMAX', nlmax, 1)
    
    ! Is the water height given by absolute or relative value
    call parlst_getvalue_int(rparlist, 'PROBLEM', 'ABSREL', absrel, 0)
    
    ! Use the bottom-sourceterm?
    call parlst_getvalue_int(rparlist, 'PROBLEM', 'bottomterm', bottomterm, 0)
    
    ! Add bottom profile to output to get position of the water surface, Default = Yes
    call parlst_getvalue_int(rparlist, 'OUTPUT', 'ADDBOTTOMTOOUT', addbottomtoout, 1)
    
    ! Apply reflecting boundary conditions?
    call parlst_getvalue_int(rparlist, 'METHOD', 'REFLECTINGBCS', reflectingbcs, 1)

    ! And with timestepsize
    call parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'dt', dt)

    ! To the final time
    call parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'ttfinal', ttfinal)

    ! Set implicitness parameter. Default=Crank-Nicolson
    call parlst_getvalue_double(rparlist, 'TIMESTEPPING', 'theta', theta, 0.5_DP)

    ! Always update Preconditioner? Default=Yes
    call parlst_getvalue_int(rparlist, 'SOLVER', &
         'alwaysupdatepreconditioner', alwaysupdatepreconditioner, 1)

    ! What kind of limiter to use? Default=Van Leer
    call parlst_getvalue_int(rparlist, 'METHOD', 'limiter', limiter, 2)

    ! What kind of syncromethod to use? Default=indicator variable
    call parlst_getvalue_int(rparlist, 'METHOD', 'syncromethod', syncromethod, 1)

    ! Gravitational constant. Default=9.81
    call parlst_getvalue_double(rparlist, 'PROBLEM', 'gravconst', gravconst, 9.81_DP)

    ! Choosing the Method to use. Default=TVD
    call parlst_getvalue_int(rparlist, 'METHOD', 'Method', Method, 3)

    ! Apply prelimiting for FCT? Default=Yes.
    call parlst_getvalue_int(rparlist, 'METHOD', 'prelimiting', prelimiting, 1)

    ! Make gmv snapshots for video? Default=No.
    call parlst_getvalue_int(rparlist, 'OUTPUT', 'makevideo', makevideo, 0)

    ! Make gmv snapshot every ... seconds (must be n*dt)
    call parlst_getvalue_double(rparlist, 'OUTPUT', 'videotimestep', videotimestep, 50000.0_DP)

    ! Maximum number of iterations for the nonlinear solver
    call parlst_getvalue_int(rparlist, 'SOLVER', 'itemax', itemax, 20)

    ! Shall boundary conditions be applied in a corner?
    call parlst_getvalue_int(rparlist, 'SOLVER', 'boundarycorner', boundarycorner, 0)

    ! Stopping criterium for the nonlinear solver
    ! (this is where you determine how accurate your solution is - this should be small)
    ! Absolute value of the norm of the defect
    call parlst_getvalue_double(rparlist, 'SOLVER', 'nonlinabsdef', nonlinabsdef, 1e-11_DP)
    ! Relative value of the norm of the defect
    call parlst_getvalue_double(rparlist, 'SOLVER', 'nonlinreldef', nonlinreldef, 1e-10_DP)
    ! Norm of solutionupdate
    call parlst_getvalue_double(rparlist, 'SOLVER', 'nonlinsolup', nonlinsolup, 1e-12_DP)

    ! Stopping criterium for linear solver
    ! (this is only for the inner defect correction loop - doesn't have to be too small)
    ! Absolute value of the norm of the defect
    call parlst_getvalue_double(rparlist, 'SOLVER', 'linabsdef', linabsdef, 1e-7_DP)
    ! Relative value of the norm of the defect
    call parlst_getvalue_double(rparlist, 'SOLVER', 'linreldef', linreldef, 1e-4_DP)
    ! The output file
    call parlst_getvalue_string (rparlist, 'OUTPUT', 'ofile', ofile, 'gmv/u2d.gmv')



    ! Create the discretisation of the domain.
    !	1. Read in the PARametrisation of the boundary
    !	2. Read in the TRIangulation
    !	3. Refine it
    !	4. Create the actual mesh from this raw data

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call parlst_getvalue_string (rparlist, 'TRIANGULATION', &
         'prmname', sstring)
    call boundary_read_prm(rboundary, sstring)

    ! Now read in the basic triangulation.
    call parlst_getvalue_string (rparlist, 'TRIANGULATION', &
         'triname', sstring)
    call tria_readTriFile2D (rtriangulation, sstring, rboundary, .true.)
    
    !call tria_rawGridToTri (rtriangulation)

    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1, rtriangulation, rboundary)

    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation, rboundary)



    ! Initialise the discretisation

    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. For the shallow water problem we need three blocks
    call spdiscr_initBlockDiscr (rdiscretisation, nvar2d, rtriangulation)

    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:

    ! Initalise the first block of the discretisation
    ! Get the type of FE from the parameter file
    call parlst_getvalue_int(rparlist, 'TRIANGULATION', 'FEkind', FEkind)
    if (FEkind == 0) then
       celement = EL_E001
    else if (FEkind == 1) then
       celement = EL_E011
    end if
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
         celement,SPDISC_CUB_AUTOMATIC, rtriangulation, &
         rboundary)

    ! Create bottom vector
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1), &
                                  rvectorbottom,.true.,st_double)

    ! and get the pointer to this vector
    call lsyssc_getbase_double (rvectorbottom,p_bottomvector)
    
    ! and make it a blockvector
    call lsysbl_createVecFromScalar (rvectorbottom,rbottomblock)
    
    ! Create source vector with the non-conservative bottom-part
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1), &
                                  rvecBSX,.true.,st_double)

    ! and get the pointer to this vector
    call lsyssc_getbase_double (rvecBSX,p_BSXdata)
    
    ! Create source vector with the non-conservative bottom-part
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1), &
                                  rvecBSY,.true.,st_double)

    ! and get the pointer to this vector
    call lsyssc_getbase_double (rvecBSY,p_BSYdata)
    

    ! Now copy this initialised block into the other ones
    ! But only create a derived copy, which shares the handles of the original one
    do ivar = 2, nvar2d
       call spdiscr_duplicateDiscrSc (rdiscretisation%Rspatialdiscr(1), &
            rdiscretisation%Rspatialdiscr(ivar), &
            .true.)
    end do



    ! Now create the system matrices

    ! First create matrix MC - the consistent mass matrix

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of some matrices we will later need to build
    ! the preconditioner of the outer defect correction loop.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our solution components
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
         LSYSSC_MATRIX9,rmatrixMC)

    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (Psi_j, Phi_i) for the
    ! scalar matrix in 2D.

    call stdop_assembleSimpleMatrix(rmatrixMC, DER_FUNC, DER_FUNC)


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
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixCX,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    call stdop_assembleSimpleMatrix(rmatrixCX, DER_DERIV_X, DER_FUNC)

    ! Now we do the same for CY
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixCY,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    call stdop_assembleSimpleMatrix(rmatrixCY, DER_DERIV_Y, DER_FUNC)
    
    ! Now get all the needed pointers to the datas of the matrices
    call lsyssc_getbase_double(rmatrixCX,p_CXdata)
    call lsyssc_getbase_double(rmatrixCY,p_CYdata)
    call lsyssc_getbase_Kcol(rmatrixMC,p_Kcol)
    call lsyssc_getbase_Kld(rmatrixMC,p_Kld)
    call lsyssc_getbase_Kdiagonal(rmatrixMC,p_Kdiagonal)


    ! Number of equations
    NEQ = rmatrixMC%NEQ


    ! Create Vector Ksep
    ! It is used to to create the array Kedge
    call storage_new ('tvd', 'Ksep', rmatrixMC%NEQ, ST_INT, &
         h_Sep, ST_NEWBLOCK_ZERO)

    ! and get Pointer to Ksep
    call storage_getbase_int (h_Sep, p_Ksep)

    ! Calculate number of edges
    nedge = (rMatrixMC%NA-rMatrixMC%neq)/2

    ! Create Array Kedge
    iSize(1)=4
    iSize(2)=nedge
    call storage_new ('tvd', 'Kedge', iSize, ST_INT, h_Kedge, &
         ST_NEWBLOCK_ZERO)

    ! and get Pointer to Kedge
    call storage_getbase_int2D (h_Kedge, p_Kedge)


    ! Now calculate the entries of Kedge
    iedge = 0
    p_Ksep(1:rmatrixMC%NEQ) = p_Kld(1:rmatrixMC%NEQ)
    ! Walk throug all lines
    do ieq = 1, rmatrixMC%neq-1
       do ij = p_Kdiagonal(ieq)+1, p_Kld(ieq+1)-1
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
       end do
    end do
    

    ! For the application of the TVD correction we will need the
    ! temporary arrays fld1+2

    ! Create Array fld1
    iSize(1)=nvar2d*6
    iSize(2)=rMatrixMC%neq
    call storage_new ('Limiter', 'fld1', iSize, ST_DOUBLE, h_fld1, &
         ST_NEWBLOCK_ZERO)

    ! and get Pointer to fld1
    call storage_getbase_double2D (h_fld1, p_fld1)

    ! Create Array fld2
    iSize(1)=nvar2d*4
    iSize(2)=nedge
    call storage_new ('Limiter', 'fld2', iSize, ST_DOUBLE, h_fld2, &
         ST_NEWBLOCK_ZERO)

    ! and get Pointer to fld2
    call storage_getbase_double2D (h_fld2, p_fld2)



    ! Now create the preconditioner block matrix P
    ! First create an empty block matrix structur with nvar2d x nvar2d blocks
    call lsysbl_createEmptyMatrix (rmatrixBlockP, nvar2d)

    ! Next create the diagonal blocks of P as empty matrices, using the
    ! matrix structur of the matrix MC
    ! We will only need the diagonal blocks, as we employ a block jacobi
    ! method here
    do ivar = 1, nvar2d
       call lsyssc_duplicateMatrix (rmatrixMC, &
            rmatrixBlockP%Rmatrixblock(ivar,ivar), &
            LSYSSC_DUP_SHARE, &
            LSYSSC_DUP_EMPTY)
    end do


    ! Create ML (lumped Mass Matrix)
    ! Copy MC into ML
    call lsyssc_duplicateMatrix (rMatrixMC,rMatrixML,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)

    ! Apply Mass Lumping
    call lsyssc_lumpMatrixScalar (rmatrixML,LSYSSC_LUMP_DIAG,.true.)

    ! Get Pointer to MC and ML Datas
    call lsyssc_getbase_double(rmatrixMC,p_MCdata)
    call lsyssc_getbase_double(rmatrixML,p_MLdata)


    ! Initialise Dij, the artifical viscosity for the artificial
    ! diffusion operator D
    Dij = 0.0_DP

    ! Create the identity matrix
    Eye = 0.0_DP
    do ivar = 1, nvar2d
       Eye(ivar,ivar) = 1.0_DP
    end do
    
    
    ! Make vector with bottom heigth values
    call storage_new ('Bottom Profile', 'Bottom Vector', neq, ST_DOUBLE, h_bottom, &
         ST_NEWBLOCK_ZERO)

    ! and get Pointer to bottom vector
    call storage_getbase_double (h_bottom, p_bottom)


    ! Now create all needed scalar and block vectors

    ! First create the block vectors
    ! rtempBlock, rrhsBlock, rsolBlock, rdefBlock and rstempBlock, rsolDotBlock
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rrhsBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rdefBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rstempBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsolDotBlock,.true.,&
         ST_DOUBLE)
    call lsysbl_createVecBlockByDiscr (rDiscretisation,rsourceBlock,.true.,&
         ST_DOUBLE)


    ! Initialise the collection to later give the solution vector to the callback routine
    ! that builds the source term and needs to evaluate the heigth field of the solution
    call collct_init(Rcollection)

    ! and hang in the pointer to the solution vector
    Rcollection%p_rvectorQuickAccess1 => rsolBlock


    ! Get pointers to the components of P, rhs, sol, rstemp, rsolDotBlock
    do ivar = 1, nvar2d
       call lsyssc_getbase_double(rmatrixBlockP%RmatrixBlock(ivar,ivar), rarrayP(ivar)%Da)	! P
       call lsyssc_getbase_double(rsolBlock%RvectorBlock(ivar), rarraySol(ivar)%Da)	! Sol
       call lsyssc_getbase_double(rrhsBlock%RvectorBlock(ivar), rarrayrhs(ivar)%Da)	! Rhs
       call lsyssc_getbase_double(rstempBlock%RvectorBlock(ivar), rarrayrstemp(ivar)%Da)	! Rstemp
       call lsyssc_getbase_double(rsolDotBlock%RvectorBlock(ivar), rarraysolDot(ivar)%Da)	! rsolDotBlock
       call lsyssc_getbase_double(rdefBlock%RvectorBlock(ivar), rarraydef(ivar)%Da)  ! Defect
       call lsyssc_getbase_double(rSourceBlock%RvectorBlock(ivar), rarraySource(ivar)%Da)  ! Source
    end do




    ! Initialise the linear solver


    ! Create a BiCGStab-solver without preconditioner.and without a filterchain
    ! as we've got our own routine to implement the boundary conditions
    nullify(p_rpreconditioner)
    nullify(p_RfilterChain)
    call linsol_initBiCGStab (p_rsolverNode)

    ! The linear solver stops, when this relative or absolut norm of
    ! the residual is reached.
    p_rsolverNode%depsRel = linreldef
    p_rsolverNode%depsAbs = linabsdef
    ! we don't have to solve the system too precise, as it is only used for
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
    call linsol_setMatrices(p_RsolverNode,Rmatrices)

    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop



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


    ! Bevor we start the timestepping we should set the initial values and read in the bottom profile
    
    ! get Pointer to the vertex coordinates
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords, p_dVertexCoords)

    ! Initialise function parser
    call fparser_init()
    call fparser_create(rfparser, 4)

    ! Get function, that describes the startvalues of h
    call parlst_getvalue_string (rparlist, 'PROBLEM', 'hstart', sstring);
    call fparser_parseFunction(rfparser, 1, trim(adjustl(sstring)), cvariables)

    ! Get function, that describes the startvalues of hu
    call parlst_getvalue_string (rparlist, 'PROBLEM', 'hustart', sstring);
    call fparser_parseFunction(rfparser, 2, trim(adjustl(sstring)), cvariables)

    ! Get function, that describes the startvalues of h
    call parlst_getvalue_string (rparlist, 'PROBLEM', 'hvstart', sstring);
    call fparser_parseFunction(rfparser, 3, trim(adjustl(sstring)), cvariables)
    
    ! Get function, that describes the bottom profile
    call parlst_getvalue_string (rparlist, 'PROBLEM', 'bottom', sbottomstring);
    call fparser_parseFunction(rfparser, 4, trim(adjustl(sbottomstring)), cvariables)

    ! Fill in initialconditions for h, hu, hv and bottom
    do ieq = 1, rmatrixMC%neq
       call fparser_evalFunction(rfparser, 1, p_DvertexCoords(:,ieq), rarraySol(1)%Da(ieq))
       call fparser_evalFunction(rfparser, 2, p_DvertexCoords(:,ieq), rarraySol(2)%Da(ieq))
       call fparser_evalFunction(rfparser, 3, p_DvertexCoords(:,ieq), rarraySol(3)%Da(ieq))
       call fparser_evalFunction(rfparser, 4, p_DvertexCoords(:,ieq), p_bottom(ieq))
       p_bottomvector(ieq) = p_bottom(ieq)
    end do

    ! Release the function parser
    call fparser_release(rfparser)
    
    
    ! Now create bottom matrices
    
    ! First create the structures
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixBX,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixBY,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
         
    call lsyssc_getbase_double(rmatrixBX,p_BXdata)
    call lsyssc_getbase_double(rmatrixBY,p_BYdata)
    
    ! Now fill the entries
!     do iedge = 1, nedge
!        i  = p_Kedge(1, iedge)
!        j  = p_Kedge(2, iedge)
!        ij = p_Kedge(3, iedge)
!        ji = p_Kedge(4, iedge)
!        ii = p_Kdiagonal(i)
!        jj = p_Kdiagonal(j)
!
! !        p_BXdata(ij) = -gravconst*p_CXdata(ij)*(p_bottom(j)-p_bottom(i))
! !        p_BYdata(ij) = -gravconst*p_CYdata(ij)*(p_bottom(j)-p_bottom(i))
! !
! !        p_BXdata(ji) = -gravconst*p_CXdata(ji)*(p_bottom(i)-p_bottom(j))
! !        p_BYdata(ji) = -gravconst*p_CYdata(ji)*(p_bottom(i)-p_bottom(j))
!
! ! Now evaluate the derivative of the bottom profile
!        Dpoints(:,1) = p_DvertexCoords(:,j)
!        Dpoints(:,2) = p_DvertexCoords(:,i)
!        call fevl_evaluate (der_deriv_x, b_x, rvectorbottom, Dpoints)
!        call fevl_evaluate (der_deriv_y, b_y, rvectorbottom, Dpoints)
!
!        p_BXdata(ij) = -gravconst*p_MCdata(ij)*b_x(1)
!        p_BYdata(ij) = -gravconst*p_MCdata(ij)*b_y(1)
!
!        p_BXdata(ji) = -gravconst*p_MCdata(ji)*b_x(2)
!        p_BYdata(ji) = -gravconst*p_MCdata(ji)*b_y(2)
!
!        p_BSXdata(i) = p_BSXdata(i)-gravconst*p_MCdata(ij)*b_x(1)
!        p_BSXdata(j) = p_BSXdata(j)-gravconst*p_MCdata(ji)*b_x(2)
!
!        p_BSYdata(i) = p_BSYdata(i)-gravconst*p_MCdata(ij)*b_y(1)
!        p_BSYdata(j) = p_BSYdata(j)-gravconst*p_MCdata(ji)*b_y(2)
!
!     end do
!
!     do i = 1, neq
!       Dpoints(:,1) = p_DvertexCoords(:,j)
!       Dpoints(:,2) = p_DvertexCoords(:,i)
!       call fevl_evaluate (der_deriv_x, b_x, rvectorbottom, Dpoints)
!       call fevl_evaluate (der_deriv_y, b_y, rvectorbottom, Dpoints)
!       ii = p_Kdiagonal(i)
!       p_BSXdata(i) = p_BSXdata(i)-gravconst*p_MCdata(ii)*b_x(2)
!       p_BSYdata(i) = p_BSYdata(i)-gravconst*p_MCdata(ii)*b_y(2)
!
!     end do
    
    
    ! If the heigth values are given as relative values, then substract the bottom profile
    ! to get the relative values, that the algorithm can work with
    if (absrel==0) then
      call SubstractBottomAfterWrite(rarraySol,neq,h_bottom)
    end if

    !  Write first video file (the initial conditions)
    ! if we want to make a video
    if (makevideo == 1) then
       ! Start UCD export to GMV file:
       call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
            'gmv/video1.gmv')
            
       ! Before writing add the bottom profile
       if (addbottomtoout==1) then
        call AddBottomBeforeWrite(rarraySol,neq,h_bottom)
       end if


       call lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
       call ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
       call lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
       call ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
       call lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
       call ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
       call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
       
        call lsyssc_getbase_double (rsourceBlock%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'Source_1',UCD_VAR_STANDARD, p_Ddata)
        call lsyssc_getbase_double (rsourceBlock%RvectorBlock(2),p_Ddata1)
        call ucd_addVariableVertexBased (rexport,'Source_2',UCD_VAR_STANDARD, p_Ddata1)
        call lsyssc_getbase_double (rsourceBlock%RvectorBlock(3),p_Ddata2)
        call ucd_addVariableVertexBased (rexport,'Source_3',UCD_VAR_STANDARD, p_Ddata2)

       ! Write the file to disc
       call ucd_write (rexport)
       call ucd_release (rexport)
       
       ! After writing substract the bottom profile
       if (addbottomtoout==1) then
        call SubstractBottomAfterWrite(rarraySol,neq,h_bottom)
       end if

       videotime = videotimestep
       ! Number of first video file
       ifilenumber = 1
    end if





    ! Begin timestepping

    ttime = 0.0_DP

    timestepping: do

       ! Compute solution from time step t^n to time step t^{n+1}
       write(*,*)
       write(*,*)
       write(*,*) 'TIME STEP:', ttime
       write(*,*)

       ! Print what method is used
       if (Method==0) then
          write(*,*), 'Calculating the high order (group FEM Galerkin) solution'
       else if (Method==1) then
          write(*,*), 'Calculating the low order (scalar dissipation) solution'
       else if (Method==2) then
          write(*,*), 'Calculating the low order (tensorial dissipation) solution'
       else if (Method==3) then
          write(*,*), 'Calculating the TVD solution'
          if (limiter==1) then
             write(*,*), 'Limiter: Minmod'
          else if (limiter==2) then
             write(*,*), 'Limiter: VanLeer'
          else if (limiter==3) then
             write(*,*), 'Limiter: MC'
          else if (limiter==4) then
             write(*,*), 'Limiter: Superbee'
          end if
       else if (Method==4) then
          write(*,*), 'char. FCT: Calculating the low order (scalar dissipation) predictor'
       else if (Method==5) then
          write(*,*), 'char. FCT: Calculating the low order (tensorial dissipation) predictor'
       else if (Method==6) then
          write(*,*), 'syncron. FCT: Calculating the low order (scalar dissipation) predictor'
       else if (Method==7) then
          write(*,*), 'syncron. FCT: Calculating the low order (tensorial dissipation) predictor'
       end if


       ! Assemble the preconditioner of the outer defect correction loop
       call BuildShallowWaterPreconditioner (rmatrixBlockP, &
            rarrayP, rarraySol, p_CXdata, p_CYdata, &
            p_MLdata, p_Kdiagonal, p_kedge, &
            NEQ, nedge, theta, dt, gravconst)

       ! Now assemble RHS
       call BuildShallowWaterRHS (&
            rarrayRhs, rarraySol, rrhsBlock, rsolBlock, &
            rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
            p_BXdata, p_BYdata, p_BSXdata, p_BSYdata, &
            h_fld1, p_fld1, p_fld2, &
            p_Kdiagonal, p_Kedge, NEQ, nedge, &
            theta, dt, gravconst, Method, limiter)

       if (bottomterm == 1) then
        call addBottomTermToVec(rSolBlock,rRhsBlock,rsourceBlock,rbottomBlock, &
                                Rcollection,(1-theta)*dt,gravconst)
       end if


       ! Here starts the defect correction loop
       def_corr: do ite = 1, itemax

          if ((alwaysupdatepreconditioner==1).and.(ite.ne.1)) then
             ! Reassemble the Preconditioner
             call BuildShallowWaterPreconditioner (rmatrixBlockP, &
                  rarrayP, rarraySol, p_CXdata, p_CYdata, &
                  p_MLdata, p_Kdiagonal, p_kedge, &
                  NEQ, nedge, theta, dt, gravconst)
          end if

          ! Assemble defect
          call BuildShallowWaterDefect (&
               rdefBlock, rstempBlock, rrhsBlock, rsolBlock, &
               rarrayRhs, rarraySol, rarrayRstemp, &
               p_CXdata, p_CYdata, p_MLdata, rmatrixML, &
               p_BXdata, p_BYdata,  p_BSXdata, p_BSYdata, &
               h_fld1, p_fld1, p_fld2, &
               p_Kdiagonal, p_Kedge, NEQ, nedge, &
               theta, dt, gravconst, Method, limiter)

          if (bottomterm == 1) then
           call addBottomTermToVec(rSolBlock,rDefBlock,rsourceBlock,rbottomBlock, &
                                   Rcollection,-theta*dt,gravconst)
          end if

          ! Take care of Boundary Conditions
          if (reflectingbcs==1) then
          call ImplementShallowWaterBCs (&
               rboundary, rtriangulation, &
               rarrayP, rarraySol, rarrayDef, &
               p_Kdiagonal, p_Kld, p_Kcol, &
               gravconst, boundarycorner)
          end if


          ! Compute norm of defect
          dcurrentDefect = lsysbl_vectorNorm (rdefBlock,LINALG_NORML2)

          ! Save the initial defect
          if (ite == 1) then
             dinitialDefect = dcurrentDefect
          end if

          ! Print norm of the current defect
          print *, "DEFCOR", ite, dcurrentDefect

          ! If norm of defect is small enough then leave the
          ! defect correction loop
          if (dcurrentDefect < nonlinabsdef .or.&
               dcurrentDefect < nonlinreldef*dinitialDefect) exit def_corr

          ! Finally solve the system.
          ! We have calculated the defect and want to calculate a
          ! solution update that we can add to the solution.
          ! So we use
          ! CALL linsol_precondDefect
          ! If we wanted to solve Au=b instead with
          ! b being the real RHS and u being the real solution vector,
          ! we would use
          ! CALL linsol_solveAdaptively (p_rsolverNode,rsolBlock,rrhsBlock,rtempBlock).
          call linsol_precondDefect (p_rsolverNode,rdefBlock)
          ! Now we have the solution update in rdef.

          ! Add the solution update to the solution: rsol=rsol+rdef
          call lsysbl_vectorLinearComb (rdefBlock,rsolBlock,1.0_dp,1.0_dp)
          
          ! As we now have dry bed handling, we clip all height variables
          ! smaller than eps
          call ClipHeight (rarraySol, neq)

          ! If norm of solution update is small enough then leave the
          ! defect correction loop
          dcurrentDefect = lsysbl_vectorNorm (rdefBlock,LINALG_NORML2)
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
!
!           ! write gmvfiles for video (if needed)
!        if (dcurrentDefect > 1000.0_dp*dinitialDefect) then
!
!           write(*,*)
!           write(*,*) 'Writing Videofile'
!           write(*,*)
!
!           write(*,*) dcurrentDefect , dinitialDefect
!           pause
!
!           ifilenumber = ifilenumber + 1
!
!           write(sfilenumber,'(i0)') ifilenumber
!
!           ! Before writing add the bottom profile
!           if (addbottomtoout==1) then
!             call AddBottomBeforeWrite(rarraySol,neq,h_bottom)
!            end if
!
!           call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                'gmv/error' // trim(sfilenumber) // '.gmv')
!
!           call lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
!           call ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
!           call lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
!           call ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
!           call lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
!           call ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
!           call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
!
!           call lsyssc_getbase_double (rsourceBlock%RvectorBlock(1),p_Ddata)
!           call ucd_addVariableVertexBased (rexport,'Source_1',UCD_VAR_STANDARD, p_Ddata)
!           call lsyssc_getbase_double (rsourceBlock%RvectorBlock(2),p_Ddata1)
!           call ucd_addVariableVertexBased (rexport,'Source_2',UCD_VAR_STANDARD, p_Ddata1)
!           call lsyssc_getbase_double (rsourceBlock%RvectorBlock(3),p_Ddata2)
!           call ucd_addVariableVertexBased (rexport,'Source_3',UCD_VAR_STANDARD, p_Ddata2)
!
!           ! Write the file to disc, that's it.
!           call ucd_write (rexport)
!           call ucd_release (rexport)
!
!           ! After writing substract the bottom profile
!           if (addbottomtoout==1) then
!             call SubstractBottomAfterWrite(rarraySol,neq,h_bottom)
!           end if
!        end if
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          if (dcurrentDefect < nonlinsolup) exit def_corr

       end do def_corr ! This is the end of the defect correction loop

       ! Test if maximum number of iterations was reached and output warning
       if (ite == itemax+1) then
          write(*,*) ''
          write(*,*) '***********************************************************************'
          write(*,*) '* WARNING: Stopping criteria of the nonlinear solver were not reached *'
          write(*,*) '***********************************************************************'
          write(*,*) ''
          ! increase the counter of the number of warnings
          numwarnings = numwarnings +1
       end if


       ! If we are using TVD, then we are finished here (at least for this timestep)
       ! But if we are using the linearised FCT, we've got only the low order predictor
       ! and need to add explicit linearised limited antidiffusion
       if ((Method==4).or.(Method==5)) then
          write(*,*) 'linearised characteristic FCT: Adding limited antidiffusion'

          ! This routine adds limited antidiffusion according to the linearised
          ! FCT limiter
          call linFctShallowWaterAddLimitedAntidiffusion_characteristic(&
               rarraySol, rarraySolDot, rarrayRhs,&
               rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
               rmatrixML, p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
               h_fld1, p_fld1, p_fld2, &
               p_Kdiagonal, p_Kedge, NEQ, nedge, &
               gravconst, dt, Method, prelimiting)

          call BuildShallowWaterPreconditioner (rmatrixBlockP, &
               rarrayP, rarraySol, p_CXdata, p_CYdata, &
               p_MLdata, p_Kdiagonal, p_kedge, &
               NEQ, nedge, theta, dt, gravconst)

          ! Take care of Boundary Conditions after this fct correction
          if (reflectingbcs==1) then
          call ImplementShallowWaterBCs (&
               rboundary, rtriangulation, &
               rarrayP, rarraySol, rarrayDef, &
               p_Kdiagonal, p_Kld, p_Kcol, &
               gravconst, boundarycorner)
          end if
       end if

       if ((Method==6).or.(Method==7)) then
          write(*,*) 'linearised syncronized FCT: Adding limited antidiffusion'

          ! This routine adds limited antidiffusion according to the linearised
          ! FCT limiter
!           call linFctShallowWaterAddLimitedAntidiffusion_syncronized(&
!                rarraySol, rarraySolDot, rarrayRhs,&
!                rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
!                rmatrixML, p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
!                h_fld1, p_fld1, p_fld2, &
!                p_Kdiagonal, p_Kedge, NEQ, nedge, &
!                gravconst, dt, Method, prelimiting, syncromethod, &
!                rtriangulation)
!
!           call new_syncronized(&
!                rarraySol, rarraySolDot, rarrayRhs,&
!                rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
!                rmatrixML, p_CXdata, p_CYdata, p_BXdata, p_BYdata, p_MLdata, p_MCdata, &
!                h_fld1, p_fld1, p_fld2, &
!                p_Kdiagonal, p_Kedge, NEQ, nedge, &
!                gravconst, dt, Method, prelimiting, syncromethod, &
!                rtriangulation)
               
          call lfctsync(rarraySol, rarraySolDot, rSolDotBlock, &
                        p_Kedge, p_Kdiagonal, NEQ, nedge,&
                        p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
                        gravconst, dt, syncromethod)
               
          ! As we now have a simple dry bed handling, we clip all height variables
          ! smaller than eps
          call ClipHeight (rarraySol, neq)

          call BuildShallowWaterPreconditioner (rmatrixBlockP, &
               rarrayP, rarraySol, p_CXdata, p_CYdata, &
               p_MLdata, p_Kdiagonal, p_kedge, &
               NEQ, nedge, theta, dt, gravconst)

          ! Take care of Boundary Conditions after this fct correction
          if (reflectingbcs==1) then
          call ImplementShallowWaterBCs (&
               rboundary, rtriangulation, &
               rarrayP, rarraySol, rarrayDef, &
               p_Kdiagonal, p_Kld, p_Kcol, &
               gravconst, boundarycorner)
          end if
       end if


       ! As we now have a simple dry bed handling, we clip all height variables
       ! smaller than eps
       call ClipHeight (rarraySol, neq)
       
!       ! Add Sourceterm explicitely
!       call AddExplicitSourceTerm(rarraySol,dt,neq,rtriangulation%h_DvertexCoords,gravconst)


       ! That's it. RvectorBlock now contains our solution at the current time

       ! write gmvfiles for video (if needed)
       if ((makevideo==1).and.(ttime.ge.videotime-0.001_DP*dt)) then

          write(*,*)
          write(*,*) 'Writing Videofile'
          write(*,*)

          ifilenumber = ifilenumber + 1

          write(sfilenumber,'(i0)') ifilenumber
          videotime = videotime + videotimestep
          
          ! Before writing add the bottom profile
          if (addbottomtoout==1) then
            call AddBottomBeforeWrite(rarraySol,neq,h_bottom)
           end if
          
          call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
               'gmv/video' // trim(sfilenumber) // '.gmv')

          call lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
          call ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
          call lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
          call ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
          call lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
          call ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
          call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
          
          call lsyssc_getbase_double (rsourceBlock%RvectorBlock(1),p_Ddata)
          call ucd_addVariableVertexBased (rexport,'Source_1',UCD_VAR_STANDARD, p_Ddata)
          call lsyssc_getbase_double (rsourceBlock%RvectorBlock(2),p_Ddata1)
          call ucd_addVariableVertexBased (rexport,'Source_2',UCD_VAR_STANDARD, p_Ddata1)
          call lsyssc_getbase_double (rsourceBlock%RvectorBlock(3),p_Ddata2)
          call ucd_addVariableVertexBased (rexport,'Source_3',UCD_VAR_STANDARD, p_Ddata2)

          ! Write the file to disc, that's it.
          call ucd_write (rexport)
          call ucd_release (rexport)
          
          ! After writing substract the bottom profile
          if (addbottomtoout==1) then
            call SubstractBottomAfterWrite(rarraySol,neq,h_bottom)
          end if
       end if

       ! Go on to the next time step
       ttime = ttime + dt

       ! Leave the time stepping loop if final time is reached
       if (ttime .ge. ttfinal-0.001_DP*dt) exit timestepping

    end do timestepping


    ! Release the memory of the solver now

    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)

    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)

    ! Release the collection
    call collct_done(Rcollection)



    ! Write the calculated solution to a gmv file

    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to GMV file:
    write(*,*)
    write(*,*) 'Writing Solution at final time',ttime,'to File'
    write(*,*)
    
    ! Before writing add the bottom profile
    if (addbottomtoout==1) then
      call AddBottomBeforeWrite(rarraySol,neq,h_bottom)
    end if
    
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
         ofile)

    ! We could write Soldot, the approximation of the time derivative
    ! of the solution, which was used while applying the linearised
    ! FCT-Limiter, to the file, too.
    !CALL ucd_addVariableVertexBased (rexport,'solDot_uh',UCD_VAR_STANDARD, rarraySolDot(1)%Da)

    ! Get pointers to the solution data
    call lsyssc_getbase_double (rsolBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol_h',UCD_VAR_STANDARD, p_Ddata)
    call lsyssc_getbase_double (rsolBlock%RvectorBlock(2),p_Ddata1)
    call ucd_addVariableVertexBased (rexport,'sol_uh',UCD_VAR_STANDARD, p_Ddata1)
    call lsyssc_getbase_double (rsolBlock%RvectorBlock(3),p_Ddata2)
    call ucd_addVariableVertexBased (rexport,'sol_vh',UCD_VAR_STANDARD, p_Ddata2)
    call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1, p_Ddata2)
    
    call lsyssc_getbase_double (rsourceBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'Source_1',UCD_VAR_STANDARD, p_Ddata)
    call lsyssc_getbase_double (rsourceBlock%RvectorBlock(2),p_Ddata1)
    call ucd_addVariableVertexBased (rexport,'Source_2',UCD_VAR_STANDARD, p_Ddata1)
    call lsyssc_getbase_double (rsourceBlock%RvectorBlock(3),p_Ddata2)
    call ucd_addVariableVertexBased (rexport,'Source_3',UCD_VAR_STANDARD, p_Ddata2)


    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! After writing substract the bottom profile
    call SubstractBottomAfterWrite(rarraySol,neq,h_bottom)


    ! Now release all memory

    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.

    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rstempBlock)
    call lsysbl_releaseVector (rsolBlock)
    call lsysbl_releaseVector (rsolDotBlock)
    call lsysbl_releaseVector (rdefBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseVector (rsourceblock)
    call lsysbl_releaseVector (rbottomblock)
    call lsysbl_releaseMatrix (rmatrixBlockP)
    

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    call storage_free (h_Sep)
    call storage_free (h_Kedge)
    call storage_free (h_fld1)
    call storage_free (h_fld2)
    call storage_free (h_bottom)
    call lsyssc_releaseMatrix (rmatrixCX)
    call lsyssc_releaseMatrix (rmatrixCY)
    call lsyssc_releaseMatrix (rmatrixBX)
    call lsyssc_releaseMatrix (rmatrixBY)
    call lsyssc_releaseMatrix (rmatrixMC)
    call lsyssc_releaseMatrix (rmatrixML)
    call lsyssc_releaseVector (rvecBSX)
    call lsyssc_releaseVector (rvecBSY)
    call lsyssc_releaseVector (rvectorbottom)

    ! Release the memory associated to the parameter file
    call parlst_done (rparlist)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)

    ! Release the triangulation.
    call tria_done (rtriangulation)

    ! Finally release the domain, that's it.
    call boundary_release (rboundary)

    ! Time measurement
    call cpu_time(dtime2)
    write(*,*) "Computational time used:", dtime2-dtime1
    write(*,*)

    if (numwarnings>0) then
          write(*,*) ''
          write(*,*) '**********************************************************************'
          write(*,*) '* WARNING: In',numwarnings, 'timesteps there where convergence problems *'
          write(*,*) '**********************************************************************'
          write(*,*) ''
      
    end if

  end subroutine shallowwater2d_0

end module shallowwater2d

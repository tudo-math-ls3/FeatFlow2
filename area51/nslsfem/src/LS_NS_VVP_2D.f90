!##############################################################################
!# ****************************************************************************
!# <name> LS_NS_VVP_2D </name>
!# ****************************************************************************
!# <purpose>
!# This module solves the 2D Navier-Stokes (NS) eq. using LSFEM.
!#
!# The second-order elliptic NS equations are reformulated into 
!# first-order equations using the definition of the vorticity:
!#   vorticity:  w = curl(u).
!# The resulting system called Velocity-Vorticity-Pressure (VVP). 
!#
!# The problem is solved in a coupled manner for the solution of:
!#   1- velocity components   u1, u2 (in 2D)
!#   2- vorticity function  w
!#   see the following link for the definition of vorticity:
!#   http://www.student.math.uwaterloo.ca/
!#         ~amat361/Fluid%20Mechanics/topics/vorticity.htm
!#   3- pressure         p
!# variables.
!#
!# The nonlinear term is first linearized using Newton method.
!# The LSFEM formulation then applied which yiedls a symmetric-
!# positive definite linear system of equations.
!# The routine uses the standard iterative\direct linear solvers.
!# The discretisation uses the block assembly method to evaluate the
!# mattrix all-in-one.
!# </purpose>
!##############################################################################

module LS_NS_VVP_2D

  use fsystem
  use storage
  use genoutput
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use linearalgebra
  use discretebc
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use spdiscprojection
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use linearsolver
  use ucd
  
  use LS_NS_callback
  
  use blockmatassemblybase
  use blockmatassembly
  use collection
  use feevaluation2
  use feevaluation
  use pprocnavierstokes
  use vectorio
  use paramlist  
  use pprocerror
  use convection 
  use matrixmodification
  
  implicit none

contains
  
  !****************************************************************************

!<subroutine>
  subroutine ls_vvp_2d
  
!<description>
  ! This is a compact LSFEM navier-stokes solver.
  ! The routine performs the following tasks:
  !
  ! 0) Read the parameters to solve the problem
  ! a) Read in triangulation and prepare grids
  ! b) Set up matrix structure
  ! c) Initialize the BCs, and the solution/RHS vectors
  ! d) Do the nonlinear iteratins, it has 6 main steps!
  ! e) Write solution to GMV/VTK file, postprocessing
  ! f) Release all variables, finish!
!</description>

!</subroutine>

  ! Definitions of variables.
  !
  ! We need a couple of variables for this problem. Let us see...
  !
  ! An object for saving the domain:
  type(t_boundary) :: rboundary

  ! An object for saving the triangulation on the domain
  type(t_triangulation) :: rtriangulation

  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation) :: rdiscretisation
  
  ! A block matrix and a couple of block vectors. These will be filled
  ! with data for the linear solver.
  type(t_matrixBlock) :: rmatrix
  type(t_vectorBlock) :: rvector,rvector_old,rrhs
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo) :: rcubatureInfo
  
  ! A variable describing the discrete boundary conditions.
  type(t_discreteBC), target :: rdiscreteBC
    
  ! Parameters used in nonliner loop
  integer :: inl,NLN_Max
  real(DP) :: dNLEpsi
  ! Convergence parameter, either by error or by NLN_Max
  logical :: converged,diverged
  
  ! Collection of vectors to evaluate in nonlinear terms
  type(t_fev2Vectors) :: revalVectors
    
  ! Viscosity parameter nu = 1/Re
  real(DP) :: dnu
  
  ! Collection structure for callback routines
  type(t_collection) :: rcollection
  
  ! All parameter of the LSFEM solver
  type(t_parlist) :: rparams 
  
  ! Ok, let us start.
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 0)-Read all the parameters from data file and initialize the collection.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call ls_initParams(rparams,NLN_MAX,dNLEpsi,rcollection)
  

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! a)-Read the domain, the mesh and refine it. All about computational GRID.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_grid(rboundary,rtriangulation,rparams,rcollection)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! b)- Set up a discretisation and cubature info structure which tells 
  ! the code which finite element to use.
  ! Also, create a 4*4 block matrix structure.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_structure(rboundary,rtriangulation,rdiscretisation,&
             rcubatureInfo,rmatrix,rparams)

  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! c)- Initialization of Boundary conditions, and the solution/RHS vectors.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize the discrete boundary conditions
  call ls_BCs_onetime(rdiscretisation,rboundary,rdiscreteBC)

  call ls_Init_RhsndSolution(rdiscretisation,rrhs,rvector_old,&
                          rvector,rparams)
                              
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! d)- Here the non-linear loop starts.
  ! The loop performs a maximum of NLN_Max iterations. The stopping criteria
  ! is controled with dNLEpsi.
  ! 
  ! Every loop performs the following series of operations:
  !   1- System matrix assembly (requires the evaluation of 
  !                nonlinear deferred velocities)
  !   1-1 And, jump stabilization is set up if required
  !   2- RHS assembly (nonlinear deferred velocities must be released 
  !          right after this step!!)
  !   3- Boundary conditions implementation, excluding the one time
  !     calculation of descretised boundary conditions which is
  !     done earlier in 'ls_BCs_onetime' subroutine
  !   4- Solver setup, solution of the final system, solver release
  !   5- Check for the non-linear loop convergence/divergence
  !   6- Update initial guess 'if (.not. converged) .and. (.not. diverged)'
  !    (all matrix and RHS vectors must be cleaned, 
  !     zero-valued, after this!!)
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do inl=1,NLN_Max

    ! Linearization Scheme
    ! If it is a mixed type, check when to shift to
    ! Newton's
    if (rcollection%IquickAccess(2) .eq. 3) then
      if (inl .gt. rcollection%IquickAccess(3)) then
        ! It's time to shift to Newton's method 
        rcollection%DquickAccess(4) = 1.0_DP
      end if
    end if
  
    ! +++++++++++++++++++++++++
    ! 1- System Matrix Assembly
    ! +++++++++++++++++++++++++
    ! Preparing a collection of vectors to evaluate nonlinear terms
    call ls_vec_collection(revalVectors,rvector_old)
    
    ! Assemble the whole system matrix 
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,ls_ns2D_Matrix,&
       rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
       revalVectors=revalVectors)
 
    ! 1-1 Set up jump stabilization
    call ls_jump(rmatrix,rparams,rcollection)
    
    
    ! +++++++++++++++
    ! 2- RHS Assembly
    ! +++++++++++++++
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,ls_ns2D_rhs,&
       rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
       revalVectors=revalVectors)  
    
    ! Release the vector structure used in linearization
    call fev2_releaseVectorList(revalVectors)
  
    
    ! ++++++++++++++++++++++++++++++++
    ! 3- Implement Boundary Conditions
    ! ++++++++++++++++++++++++++++++++
    call ls_BCs_iterative(rmatrix,rrhs,rvector,rvector_old,rdiscreteBC)
    
    
    ! +++++++++++++++++++
    ! 4- Solve the System
    ! +++++++++++++++++++
    call ls_Solver_iterative(rmatrix,rrhs,rvector)
    
    
    ! ++++++++++++++++++++++++++++++++++++++
    ! 5- Check for Convergence or Divergence
    ! ++++++++++++++++++++++++++++++++++++++
    call ls_con_di_verge(converged,diverged,rvector,&
                rvector_old,NLN_Max,inl,dNLEpsi)
    
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 6- Update the Initial Solution, clear matrix, vectors
    !  or Terminate the Loop
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ((.not. converged) .and. (.not. diverged)) then
    
      ! Copy the current solution to the old solution
      call lsysbl_copyVector (rvector,rvector_old)
      !*** Clear all data in matrix and RHS ***!
      call lsysbl_clearVector (rrhs)
      call lsysbl_clearMatrix (rmatrix)
    
    else
    
      if (diverged) then
        call output_lbrk()
        call output_line ('Nonlinear Loop is diverged :(')
        call output_lbrk()
      end if
      
      EXIT
    
    end if
        
  end do  ! inl


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! e)- Post-processing starts here. Writing the solution to file ...
  ! export GMV/VTK files and calculate drag/lift forces ...
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_postprocess(rboundary,rmatrix,rvector,rtriangulation,&
                    rcubatureInfo,rdiscretisation,rparams)
  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! f)- Clean up, free all the memory which is used for variables.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_cleanup(rvector,rvector_old,rrhs,rmatrix,rcubatureInfo,&
             rdiscreteBC,rdiscretisation,rtriangulation,rboundary)

  end subroutine


  !****************************************************************************
  
!<subroutine>
  subroutine ls_initParams(rparams,NLN_MAX,dNLEpsi,rcollection)

 !<description>  
  ! In this subroutine, the collection structrure is initialized.
 !</description>                

 !<output>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(out) :: rparams 

  ! Nonlinear loop stopping criteria
  real(DP), intent(out) :: dNLEpsi

  ! Nonlinear loop Max. number of iterations
  integer, intent(out) :: NLN_MAX
  
  ! Collection structure for callback routines  
  type(t_collection), intent(out) :: rcollection
 !</output>

 !</subroutine>

  ! Local variables
  ! Kinematic Viscosity 
  real(DP) :: dnu
  
  ! determine whether scaling factors (physical and mesh dependent)
  ! should be multiplied by our least-squares functionals. 
  integer :: scPhysic, scADN
  
  ! Continuity equation scaling factor
  !  default = 1.0_DP  no scaling
  real(DP) :: alpha

  ! Linearization method parameters
  integer :: LinTyp, FPIter
  
  ! reading data from the *.dat file 
  call parlst_init(rparams)
  call parlst_readfromfile (rparams, "./data/lsvvp.dat") 
     
  ! Viscosity parameter: noo = 1/Re
  call parlst_getvalue_double (rparams, 'GFPROPER', 'dnu', dnu, 1.0_DP)   
  
  ! Nonlinear loop Max. number of iterations
  call parlst_getvalue_int (rparams, 'NLLOOP', 'NLN_MAX', NLN_MAX, 3)
  
  ! Nonlinear loop stopping criteria
  ! Serves as a measure of relative residual.
  call parlst_getvalue_double (rparams, 'NLLOOP', 'dNLEpsi', dNLEpsi, 0.001_DP)  
  
  ! Initializing the collections
  ! Put kinematic viscosity there, to be used in nonlinear assembly
  rcollection%DquickAccess(1) = dnu

  ! Scale least-squares functionals or not.
  call parlst_getvalue_int (rparams, 'GFPROPER', 'scPhysic', scPhysic, 0)
  call parlst_getvalue_int (rparams, 'GFPROPER', 'scADN', scADN, 0)
  call parlst_getvalue_double (rparams, 'GFPROPER', 'alpha', alpha, 1.0_DP)
  
  ! Linearization method
  call parlst_getvalue_int (rparams, 'NLLOOP', 'LinTyp', LinTyp, 0)
  call parlst_getvalue_int (rparams, 'NLLOOP', 'FPIter', FPIter, 0)
  
  if (scPhysic .eq. 1) then
    rcollection%DquickAccess(2) = 1.0_DP/(dnu*dnu)
  else
    rcollection%DquickAccess(2) = 1.0_DP
  end if
  
  rcollection%IquickAccess(1) = scADN

  ! Continuity equation scaling
  !  \alpha * || \nabla \cdot \bu ||_0
  rcollection%DquickAccess(3) = alpha
    
  ! Linearization Scheme
  rcollection%IquickAccess(2) = LinTyp
  select case (LinTyp)
  case (1)

   ! Fixed-point
   rcollection%DquickAccess(4) = 0.0_DP

  case (2)

   ! Newton's method
   rcollection%DquickAccess(4) = 1.0_DP

  case (3)

   ! Mixed method
   rcollection%DquickAccess(4) = 0.0_DP
   
   if (FPIter .le. NLN_MAX) then
    rcollection%IquickAccess(3) = FPIter
   else
    rcollection%IquickAccess(3) = NLN_MAX
   end if   

  end select
    
  end subroutine
  

  !****************************************************************************
  
!<subroutine>
  subroutine ls_grid(rboundary,rtriangulation,rparams,rcollection)

 !<description>  
  ! In this subroutine the initial mesh is read and all the postprocessing
  ! steps are done to prepare the final mesh.
 !</description>                

 !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams 
 !</input> 

 !<output>
  ! An object for saving the domain:
  type(t_boundary), intent(out) :: rboundary

  ! An object for saving the triangulation on the domain
  type(t_triangulation), intent(out) :: rtriangulation
 !</output>

 !<inputoutput>
  ! Collection structure for callback routines  
  type(t_collection), intent(inout) :: rcollection
 !</inputoutput>


 !</subroutine>

  ! Local variables
  ! Path to the mesh  
  character(LEN=SYS_STRLEN) :: sfile

  ! NLMAX, the fine level to solve the problem
  integer :: NLMAX, ihandle
    
  ! We want to solve our NS problem on level...
  call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 3)
  
  call parlst_getvalue_string (rparams, 'MESH', &
                 'sFilenamePathMesh',sfile, """""", bdequote=.true.)  


  ! At first, read in the parametrisation of the boundary and save
  ! it to rboundary.
  call boundary_read_prm(rboundary, trim(sfile)//'.prm')
    ! bench1
  ! Now read in the basic triangulation.
  call tria_readTriFile2D (rtriangulation, trim(sfile)//'.tri', rboundary)
  
  ! Refine the mesh up to the maximum level
  call tria_quickRefine2LevelOrdering(NLMAX-1,rtriangulation,rboundary)
  
  ! Create information about adjacencies and everything one needs from
  ! a triangulation. Afterwards, we have the coarse mesh.
  call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

!  ihandle = rtriangulation%h_Delementvolume
!  call storage_getbase_double (ihandle, rcollection%p_Rvectorquickaccess1)
   
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_structure(rboundary,rtriangulation,rdiscretisation,&
             rcubatureInfo,rmatrix,rparams)

 !<description>  
  ! Set up a discretisation and cubature info structure. This tells
  ! the code which finite element to use.
  ! Also, create a block matrix structure.
 !</description>                

 !<input>
  ! An object for saving the domain:
  type(t_boundary), intent(in) :: rboundary

  ! An object for saving the triangulation on the domain
  type(t_triangulation), intent(in) :: rtriangulation

  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams 
 !</input>

 !<output>
  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation), intent(out) :: rdiscretisation
  
  ! A block matrix which will be filled
  ! with data for the linear solver.
  type(t_matrixBlock), intent(out) :: rmatrix
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(out) :: rcubatureInfo  
 !</output>
 
!</subroutine>
 
 
  ! local variables
  ! String variable
  character(len=SYS_STRLEN) :: sstring
    
  ! Type of finite elements
  integer(I32) :: Velm, Pelm, Welm  
  
  ! Jump stabilization parameters
  integer :: detVJump, detWJump

  ! Type of cubature rule for numerical integration
  integer(I32) :: ccubType
  
  ! Now we can start to initialise the discretisation. At first, set up
  ! a block discretisation structure that specifies 4 blocks in the
  ! solution vector.
  call spdiscr_initBlockDiscr (rdiscretisation,4,&
                 rtriangulation, rboundary)

  ! rdiscretisation%RspatialDiscr is a list of scalar
  ! discretisation structures for every component of the solution vector.
  ! We have a solution vector with four components:
  !  Component 1 = X-velocity
  !  Component 2 = Y-velocity
  !  Component 3 = Pressure
  !  Component 4 = Vorticity   
  ! For simplicity, we set up one discretisation structure for the
  ! velocity...
  
  ! Read the finite element for velocities
  call parlst_getvalue_string (rparams, 'MESH', 'Velm', sstring)
  Velm = elem_igetID(sstring)
!  call parlst_getvalue_int (rparams, 'GENERAL', 'Vtild', Vtild) 
       
  call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
        Velm, rtriangulation, rboundary)
        
  ! ...and copy this structure also to the discretisation structure
  ! of the 2nd component (Y-velocity). This needs no additional memory,
  ! as both structures will share the same dynamic information afterwards.
  call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
    rdiscretisation%RspatialDiscr(2))

  ! For the pressure (3rd component), we set up a separate discretisation
  ! structure, as this 'MAY' use different finite elements for trial and test
  ! functions.
  ! Read the finite element for Pressure
  call parlst_getvalue_string (rparams, 'MESH', 'Pelm', sstring)
  Pelm = elem_igetID(sstring)
     
  call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
    Pelm, rdiscretisation%RspatialDiscr(3))
  
  ! And for Vorticity, a separate doscretisation structure as well.
  ! Read the finite element for Vorticity
  call parlst_getvalue_string (rparams, 'MESH', 'Welm', sstring)
  Welm = elem_igetID(sstring)  
    
  call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
    Welm, rdiscretisation%RspatialDiscr(4))

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Set up an cubature info structure to tell the code which cubature
  ! formula to use
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! Create an assembly information structure which tells the code
  ! the cubature formula to use.
  ! Get the integration rule from the data file
  call parlst_getvalue_string (rparams,'MESH','ccubType', sstring)
  ccubType = cub_igetID(sstring)  
  call spdiscr_createDefCubStructure(&  
    rdiscretisation%RspatialDiscr(1),rcubatureInfo, ccubType)
  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Create a block matrix structure for the system matrix
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Initialise the block matrix with default values based on
  ! the discretisation.
  call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)
  
  ! Now as the discretisation is set up, we can start to generate
  ! the structure of the system matrix which is to solve.
  !
  ! The global system looks like this, a full 4*4 block matrix
  ! which is dymmetric and positive definite.
  !
  !  ( A11 A12 A13 A14 )
  !  ( A21 A22 A23 A24 )
  !  ( A31 A32 A33 A34 )
  !  ( A41 A42 A43 A44 )
  !  
  ! Create the matrix structure of the X-velocity. Block A11,
  !
  ! Let's check if we have to set up jump stabilization
  ! If so, we need to define the matrix structure accordingly
  call parlst_getvalue_int (rparams, 'JUMP', 'detVJump', detVJump, 0)  

  ! Velocity jump stabilization
  if (detVJump .eq. 1) then  
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                     LSYSSC_MATRIX9, &
        rmatrix%RmatrixBlock(1,1),cconstrType=BILF_MATC_EDGEBASED)
  else 
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,1))  
  end if
                   
  ! Block A12,
  call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                   LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,2))  
  ! Block A13
  call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(3),&
                   LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,3),&
                   rdiscretisation%RspatialDiscr(1))     
  ! Block A14
  call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(4),&
                   LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,4),&
                   rdiscretisation%RspatialDiscr(1))  
                                  
                                  
  ! Use X-velocity structure for the Y-velocity. Block A22,
  call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
  ! Block A23,
  call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,3),&
    rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
  ! Block A24,
  call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,4),&
    rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)


  ! Create the matrix structure of the Pressure. Block A33,
  call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(3),&
                   LSYSSC_MATRIX9, rmatrix%RmatrixBlock(3,3))
  ! Block A34,
  call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(4),&
                   LSYSSC_MATRIX9, rmatrix%RmatrixBlock(3,4),&
                   rdiscretisation%RspatialDiscr(3))


  ! Create the matrix structure of the Vorticity. Block A44,
  ! 
  ! Let's check if we have to set up jump stabilization
  ! If so, we need to define the matrix structure accordingly
  call parlst_getvalue_int (rparams, 'JUMP', 'detWJump', detWJump, 0)  

  ! Velocity jump stabilization
  if (detWJump .eq. 1) then   
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(4),&
                     LSYSSC_MATRIX9, &
        rmatrix%RmatrixBlock(4,4),cconstrType=BILF_MATC_EDGEBASED)
  else 
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(4),&
                LSYSSC_MATRIX9, rmatrix%RmatrixBlock(4,4))  
  end if
    
    
  ! Create the structure for the transpoed matrices by transposing the structre
  ! of the relevant matrices. Block A21,
  call lsyssc_transposeMatrix (&
    rmatrix%RmatrixBlock(1,2),rmatrix%RmatrixBlock(2,1),LSYSSC_TR_STRUCTURE)
  
  
  ! Block A31,
  call lsyssc_transposeMatrix (&
    rmatrix%RmatrixBlock(1,3),rmatrix%RmatrixBlock(3,1),LSYSSC_TR_STRUCTURE)
  ! Block A32,
  call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(3,1),&
    rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    
    
  ! Block A41,
  call lsyssc_transposeMatrix (&
    rmatrix%RmatrixBlock(1,4),rmatrix%RmatrixBlock(4,1),LSYSSC_TR_STRUCTURE)
  ! Block A42,
  call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(4,1),&
    rmatrix%RmatrixBlock(4,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
  ! Block A43  
  call lsyssc_transposeMatrix (&
    rmatrix%RmatrixBlock(3,4),rmatrix%RmatrixBlock(4,3),LSYSSC_TR_STRUCTURE)
  
  ! Now re-assign the block discretisation structure to all matrices
  call lsysbl_assignDiscrDirectMat (rmatrix,rdiscretisation)
  
  ! Allocate memory for the matrix
  call lsysbl_allocEmptyMatrix (rmatrix,LSYSSC_SETM_ZERO)

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_onetime(rdiscretisation,rboundary,rdiscreteBC)
                
 !<description>  
  ! In this subroutine we discretise the boundary conditions and
  ! prepare them to be applied to the matrix/RHS/sol in the 
  ! nonlinear iteration loop.
 !</description>                

 !<output>
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC),intent(out), target :: rdiscreteBC
 !</output>

 !<input>
  ! An object for saving the domain:
  type(t_boundary), intent(in) :: rboundary

  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation), intent(in) :: rdiscretisation  
 !</input>

!</subroutine>

  ! Local variables
  type(t_boundaryRegion) :: rboundaryRegion

  ! For implementing boundary conditions, we use a `filter technique with
  ! discretised boundary conditions`. This means, we first have to calculate
  ! a discrete version of the analytic BC, which we can implement into the
  ! solution/RHS vectors using the corresponding filter.
  !
  ! Create a t_discreteBC structure where we store all discretised boundary
  ! conditions.
  call bcasm_initDiscreteBC(rdiscreteBC)
  
  ! We first set up the boundary conditions for the X-velocity, then those
  ! of the Y-velocity.
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
                     getBoundaryValues_2D)
               
!  ! edge 2 of boundary component 1.
!  call boundary_createregion(rboundary,1,2,rboundaryregion)
!  call bcasm_newdirichletbconrealbd (rdiscretisation,1,&
!                     rboundaryregion,rdiscretebc,&
!                     getboundaryvalues_2d)
               
  ! Edge 3 of boundary component 1.
  call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!  rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)
  
  ! Edge 4 of boundary component 1. That is it.
  call boundary_createRegion(rboundary,1,4,rboundaryRegion)
  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)

!  ! Edge 5 of boundary component 1. That is it.
!  call boundary_createRegion(rboundary,1,5,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)

  ! Now continue with defining the boundary conditions of the Y-velocity:
  !
  ! Define edge 1.
  call boundary_createRegion(rboundary,1,1,rboundaryRegion)
  ! As we define the Y-velocity, we now set icomponent=2 in the following call.
  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)
               
!  ! Edge 2 of boundary component 1.
!  call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
               
  ! Edge 3 of boundary component 1.
  call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!  rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)
  
  ! Edge 4 of boundary component 1. That is it.
  call boundary_createRegion(rboundary,1,4,rboundaryRegion)
  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)


!  ! Edge 5 of boundary component 1. That is it.
!  call boundary_createRegion(rboundary,1,5,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)

  ! Pressure BCs.
  ! The pressure and vorticity do not need boundary conditions.
  ! Edge 2 of boundary component 1. That is it.
  call boundary_createRegion(rboundary,1,2,rboundaryRegion)
  call bcasm_newDirichletBConRealBD (rdiscretisation,3,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)


  ! Flow around cylinder
  ! The pressure and vorticity do not need boundary conditions.
  ! Edge 2 of boundary component 2. That is it.
  call boundary_createRegion(rboundary,2,1,rboundaryRegion)
  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)


  ! The pressure and vorticity do not need boundary conditions.
  ! Edge 2 of boundary component 2. That is it.
  call boundary_createRegion(rboundary,2,1,rboundaryRegion)
  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
                     getBoundaryValues_2D)

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Init_RhsndSolution(rdiscretisation,rrhs,rvector_old,&
      rvector,rparams)
                
 !<description>  
  ! Initializing the RHS and the solution vector
 !</description>                

 !<output>
  ! Block vectors
  type(t_vectorBlock), intent(out) :: rvector_old,rvector,rrhs
 !</output>
 
 !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams 

  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
 !</input>

!</subroutine>

  ! local variables
  integer :: nlinit

  ! Initial value for the 1st step of nonliner loop
  real(DP) :: dinit_vect(4)

  ! Path to the data file which has the initial solution
  character(LEN=SYS_STRLEN) :: sfile, sstring, sarray
  
  ! Create a RHS and a solution vector(s) based on the discretisation.
  ! Fill with zero.
  call lsysbl_createVectorBlock (rdiscretisation,rrhs,.true.)
  call lsysbl_createVectorBlock (rdiscretisation,rvector_old,.true.)
  call lsysbl_createVectorBlock (rdiscretisation,rvector,.true.)   
  
  ! Determine how to setup initial nonlinear solution
  call parlst_getvalue_int (rparams, 'ISOLUTION', 'nlinit', nlinit, 0)  
  
  if (nlinit .eq. 0) then
  
    ! Initialize the RHS and solution vector(s)
    call parlst_getvalue_string (rparams, 'ISOLUTION', 'initValues',&
                   sstring, '0.0_DP 0.0_DP 0.0_DP 0.0_DP')
    read (sstring,*) dinit_vect(1), dinit_vect(2), dinit_vect(3), dinit_vect(4)
    ! Scale the sub-vectors to initialize the nonlineaer iteration loop
    call lsyssc_clearVector (rvector_old%RvectorBlock(1),dinit_vect(1))
    call lsyssc_clearVector (rvector_old%RvectorBlock(2),dinit_vect(2))
    call lsyssc_clearVector (rvector_old%RvectorBlock(3),dinit_vect(3))
    call lsyssc_clearVector (rvector_old%RvectorBlock(4),dinit_vect(4))    
    
  else      
  
    ! Ignor the initial values, read from file
    call parlst_getvalue_string (rparams, 'ISOLUTION', &
           'sFilenamePathResult',sfile, """""", bdequote=.true.)  
    call vecio_readBlockVectorHR (rvector_old, sarray, .true.,&
                      0, sfile, .true.)
                      
  end if
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_iterative(rmatrix,rrhs,rvector,rvector_old,rdiscreteBC)
                
 !<description>  
  ! Implementing BCs to the matrix and solution/RHS vectors.
 !</description>                

 !<inputoutput>
  ! Block matrix and vectors
  type(t_matrixBlock), intent(inout) :: rmatrix  
  type(t_vectorBlock), intent(inout) :: rvector,rvector_old,rrhs
 !</inputoutput>
 
 !<input>
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC),intent(inout), target :: rdiscreteBC
 !</input>

!</subroutine>

  ! Assign the boundary conditions to the matrix and the vectors.
  call lsysbl_assignDiscreteBC(rmatrix,rdiscreteBC)
  call lsysbl_assignDiscreteBC(rrhs,rdiscreteBC)
  call lsysbl_assignDiscreteBC(rvector,rdiscreteBC)
  call lsysbl_assignDiscreteBC(rvector_old,rdiscreteBC)
  
  ! Next step is to implement boundary conditions into the RHS,
  ! solution and matrix. This is done using a vector/matrix filter
  ! for discrete boundary conditions.
  ! The discrete boundary conditions are already attached to the
  ! vectors/matrix. Call the appropriate vector/matrix filter that
  ! modifies the vectors/matrix according to the boundary conditions.
  call vecfil_discreteBCrhs (rrhs)
  call vecfil_discreteBCsol (rvector)
  call vecfil_discreteBCsol (rvector_old)
  call matfil_discreteBC (rmatrix)
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Solver_iterative(rmatrix,rrhs,rvector)
                
 !<description>  
  ! Set up a linear solver, solve the problem, release the solver.
 !</description>                

 !<inputoutput>
  ! Solution Vector  
  type(t_vectorBlock), intent(inout) :: rvector  
 !</inputoutput>
 
  !<input>
   ! Block matrix and vectors
  type(t_matrixBlock), intent(inout) :: rmatrix  
  type(t_vectorBlock), intent(inout) :: rrhs  
  !</input>
 
!</subroutine>

  ! Local variables
  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  type(t_matrixBlock), dimension(1) :: Rmatrices
  
  ! A temporary vector
  type(t_vectorBlock) :: rtempBlock
  
  ! Error indicator during initialisation of the solver
  integer :: ierror  
  
  ! A solver node that accepts parameters for the linear solver 
  type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner  
  
  ! A filter chain that describes how to filter the matrix/vector
  ! before/during the solution process. The filters usually implement
  ! boundary conditions.
  type(t_filterChain), dimension(2), target :: RfilterChain

  integer, dimension(1) :: Irows = (/1/)
 
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! During the linear solver, the boundary conditions must
  ! frequently be imposed to the vectors. This is done using
  ! a filter chain. As the linear solver does not work with
  ! the actual solution vectors but with defect vectors instead,
  ! a filter for implementing the real boundary conditions
  ! would be wrong.
  ! Therefore, create a filter chain with one filter only,
  ! which implements Dirichlet-conditions into a defect vector.
!  RfilterChain(1)%ifilterType = FIlTER_DISCBCDEFREAL
  
!  RfilterChain(2)%ifilterType = FILTER_TOL20
!  RfilterChain(2)%itoL20component = 3
  
  ! Create a BiCGStab-solver with VANCA preconditioner.
  ! Attach the above filter chain to the solver, so that the solver
  ! automatically filters the vector during the solution process.
!  nullify(p_rpreconditioner)
!  call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
!  call linsol_initSSOR (p_rpreconditioner, 1.0_DP, .true.)
!  call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)

  ! We will allow the solver to perform 1000 iterations
!  p_rsolverNode%nmaxIterations = 200
!  p_rsolverNode%depsRel = 1E-5

  call linsol_initUMFPACK4 (p_rsolverNode)

  ! Set the output level of the solver to 2 for some output
  p_rsolverNode%ioutputLevel = -1
  
  ! Attach the system matrix to the solver.
  ! First create an array with the matrix data (on all levels, but we
  ! only have one level here), then call the initialisation
  ! routine to attach all these matrices.
  ! Remark: Do not make a call like
  !  CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
  ! This does not work on all compilers, since the compiler would have
  ! to create a temp array on the stack - which does not always work!


!  ! Making pressure unique
!  call vecfil_subvectorToL20 (rrhs,3)
!  call mmod_replaceLinesByUnitBlk (rmatrix,3,Irows)  
  
  
  Rmatrices = (/rmatrix/)
  call linsol_setMatrices(p_rsolverNode,Rmatrices)

  call linsol_initStructure (p_rsolverNode, ierror)

  if (ierror .ne. LINSOL_ERR_NOERROR) then
    call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
    call sys_halt()
  end if

  ! Initialise structure/data of the solver. This allows the
  ! solver to allocate memory / perform some precalculation
  ! to the problem.
  call linsol_initData (p_rsolverNode, ierror)

  if (ierror .ne. LINSOL_ERR_NOERROR) then
    call output_line("Matrix singular!",OU_CLASS_ERROR)
    call sys_halt()
  end if

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Solve the system
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! Finally solve the system. As we want to solve Ax=b with
  ! b being the real RHS and x being the real solution vector,
  ! we use linsol_solveAdaptively. If b is a defect
  ! RHS and x a defect update to be added to a solution vector,
  ! we would have to use linsol_precondDefect instead.
  call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)  

  ! Release solver data and structure
  call linsol_doneData (p_rsolverNode)
  call linsol_doneStructure (p_rsolverNode)
  
  ! Release the solver node and all subnodes attached to it (if at all):
  call linsol_releaseSolver (p_rsolverNode)
  
  ! Release temporary vector
  call lsysbl_releaseVector (rtempBlock)
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_vec_collection(revalVectors,rvector_old)
                
 !<description>  
  ! Calculate the vector collection to be used in nonlinear terms
  ! Nonlinear deferred velocities/derivatives are calculated in all
  ! cubature points of the given element set.
 !</description>                

 !<output>
  ! Values of FEM functions automatically evaluated in the
  ! cubature points.
  type(t_fev2Vectors), intent(inout) :: revalVectors 
 !</output>
 
  !<input>
   ! Solution Vector in the current nonliner iteration  
  type(t_vectorBlock), intent(inout) :: rvector_old
  !</input>
 
!</subroutine>
 
  ! The routine <verb>fev2_addVectorToEvalList</verb> allows to define
  ! the evaluation of derivatives.

  call fev2_addVectorToEvalList(revalVectors,&
     rvector_old%RvectorBlock(1),1)   ! u1
  call fev2_addVectorToEvalList(revalVectors,&
     rvector_old%RvectorBlock(2),1)   ! u2
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_jump(rmatrix,rparams,rcollection)
                
 !<description>  
  ! Set up the jump stabilization, if required. Based on the parameters 
  ! in the section [JUMP], this subroutine may be active or not.
 !</description>                

 !<inputoutput>
  ! Block matrix
  type(t_matrixBlock), intent(inout) :: rmatrix  
 !</inputoutput>
 
  !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams

  ! Collection structure for callback routines
  type(t_collection), intent(in) :: rcollection  
  !</input>
 
!</subroutine>
 
  ! Local variables

  ! Jump stabilization structure
  type(t_jumpStabilisation) :: rjumpStabil
  
  ! Jump stabiliztion parameters
  integer :: detVJump
  integer :: detWJump
  real(DP) :: dJumpV, dJumpStarV, deojEdgeExpV
  real(DP) :: dJumpW, dJumpStarW, deojEdgeExpW
  
  ! Let's check if we realy have to set up jump stabilization
  call parlst_getvalue_int (rparams, 'JUMP', 'detVJump', detVJump, 0)
  call parlst_getvalue_int (rparams, 'JUMP', 'detWJump', detWJump, 0)  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Velocity jump stabilization
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (detVJump .eq. 1) then
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpV', &
                          dJumpV, 0.01_DP)  
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpStarV',&
                        dJumpStarV, 0.0_DP)  
    call parlst_getvalue_double (rparams, 'JUMP', 'deojEdgeExpV',&
                        deojEdgeExpV, 2.0_DP)                            

    ! Set up the jump stabilisation structure.
    ! The kinematic viscosity 1/Re
    rjumpStabil%dnu = rcollection%DquickAccess(1)

    ! Set stabilisation parameter
    rjumpStabil%dgamma = dJumpV
    rjumpStabil%dgammastar = dJumpStarV
    rjumpStabil%deojEdgeExp = deojEdgeExpV

    ! Matrix weight, =0 no jump stabilization will be added
    rjumpStabil%dtheta = 1.0_DP

    ! Cubature formula to be used in jump term calculations
    ! over the edges
    rjumpStabil%ccubType = CUB_G2_1D

    ! Call the jump stabilisation technique for the 1st velocity.
    call conv_jumpStabilisation2d (rjumpStabil, CONV_MODMATRIX, &
                      rmatrix%RmatrixBlock(1,1)) 
    
    ! Call the jump stabilisation technique for the 2nd velocity.
    call conv_jumpStabilisation2d (rjumpStabil, CONV_MODMATRIX, &
                      rmatrix%RmatrixBlock(2,2)) 

  end if
  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Vorticity jump stabilization
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (detWJump .eq. 1) then
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpW', &
                          dJumpW, 0.01_DP)  
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpStarW',&
                        dJumpStarW, 0.0_DP)  
    call parlst_getvalue_double (rparams, 'JUMP', 'deojEdgeExpW',&
                        deojEdgeExpW, 2.0_DP)                            

    ! Set up the jump stabilisation structure.
    ! The kinematic viscosity 1/Re
    rjumpStabil%dnu = rcollection%DquickAccess(1)

    ! Set stabilisation parameter
    rjumpStabil%dgamma = dJumpW
    rjumpStabil%dgammastar = dJumpStarW
    rjumpStabil%deojEdgeExp = deojEdgeExpW

    ! Matrix weight, =0 no jump stabilization will be added
    rjumpStabil%dtheta = 1.0_DP

    ! Cubature formula to be used in jump term calculations
    ! over the edges
    rjumpStabil%ccubType = CUB_G2_1D

    ! Call the jump stabilisation technique for the 1st velocity.
    call conv_jumpStabilisation2d (rjumpStabil, CONV_MODMATRIX, &
                      rmatrix%RmatrixBlock(4,4)) 
    
  end if  
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_con_di_verge(converged,diverged,rvector,rvector_old,&
               NLN_Max,inl,dNLEpsi)
                
 !<description>  
  ! Check the nonlinear loop convergence/divergence.
  ! Print the residuals.
 !</description>                

 !<output>
  ! Convergence parameter, either by error or by NLN_Max
  logical, intent(out) :: converged,diverged
 !</output>
 
  !<input>
   ! Solution vectors in the current/previous nonliner iterations  
  type(t_vectorBlock) :: rvector,rvector_old
  
  ! Nonlinear loop's maximum/current number of iterations
  integer, intent(in) :: NLN_Max,inl
  
  ! Nonlinear loop's stopping criteria for the norms
  real(DP), intent(in) :: dNLEpsi 
  !</input>
 
!</subroutine>
 
   ! Local variables
  real(DP) :: Dres(2),Dresv(2),Dres_rel(2)
  integer, dimension(2) :: Cnorms
  integer :: i,isum
  
  ! Scaling factors
  real(DP) :: cx,cy
  
  ! Difference of the vectors in the current nonlinear iteration
  type(t_vectorBlock) :: rdiff
  
  
  ! Euclidian vector norm: (vector,vector) 0
  ! $l_2$-norm: 1/sqrt(NEQ) * (vector,vector) 2
  ! max-norm 3   
  ! Normal L^2 Norm
  Cnorms(:) = 0
  
  ! Initialize the 'rdiff' structure and set the values to zero
  call lsysbl_createVecBlockIndirect (rvector,rdiff,.true.)
  
  ! Perform a linear combination: rdiff = cx * rvector  +  cy * rvector_old
  cx = 1.0_DP
  cy = -1.0_DP
  call lsysbl_vectorLinearComb (rvector,rvector_old,cx,cy,rdiff)
  
  ! Calculate the norm of the difference of the velocity sub-vectors
  call lsysbl_vectorNormBlock (rdiff,Cnorms,Dres)

  ! Calculate the norm of all current iteration velocity sub-vectors
  call lsysbl_vectorNormBlock (rvector,Cnorms,Dresv)
  
  ! Calculate the relative error of velocity sub-vectors 
  Dres_rel(1) = Dres(1)/Dresv(1)
  Dres_rel(2) = Dres(2)/Dresv(2)
  
  ! Convergence check
  converged = .false.
  isum = 0   
  ! Iteration number control
  if (inl .eq. NLN_Max) then
    converged = .true.
  else
    ! Norm control
    do i=1,2
      if (Dres_rel(i) .lt. dNLEpsi) then
         isum = isum + 1
      end if
    end do
    if (isum .eq. 2) then
      converged = .true.
    end if
  end if  

  ! Divergence check
  diverged = .false.
  diverged = .not. (Dres_rel(1) .lt. 1E5 .and. Dres_rel(2) .lt. 1E5)

  
  ! Release the block vector
  call lsysbl_releaseVector (rdiff)  


  ! Some output data
  if (inl .eq. 1) then
    call output_line ('Iter. ' //' U1 Rel. Err. ' //' U2 Rel. Err. ')
    call output_line ('----------------------------------')
    call output_line (sys_siL(inl, 5) //'  '&
    //trim(sys_sdEL(Dres_rel(1),6))//'  '&
    //trim(sys_sdEL(Dres_rel(2),6)))   
  else
    call output_line (sys_siL(inl, 5) //'  '&
    //trim(sys_sdEL(Dres_rel(1),6))//'  '&
    //trim(sys_sdEL(Dres_rel(2),6)))   
    if ( (mod(inl,10) .eq. 0) .and. (inl .ne. NLN_Max) &
      .and. (.not. converged) .and. (.not. diverged)) then
      call output_lbrk()
      call output_line ('Iter. ' &
      //' U1 Rel. Err. ' //' U2 Rel. Err. ')
      call output_line ('----------------------------------')
    end if
  end if

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_postprocess(rboundary,rmatrix,rvector,rtriangulation,&
                    rcubatureInfo,rdiscretisation,rparams)
                
 !<description>  
  ! Postprocessing of the LSFEM solution is done here.
  ! Writing the GMV/VTK files, calculating the darg/lift forces, ...
 !</description>                

  !<input>
   ! Final solution vector
  type(t_vectorBlock), intent(in) :: rvector
  
  ! An object for saving the triangulation on the domain
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! An object for saving the domain:
  type(t_boundary), intent(in) :: rboundary
    
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo

  ! An object specifying the discretisation.
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
  
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams   
  !</input>
 
  !<inputoutput>  
  ! The system block matrix.
  type(t_matrixBlock), intent(inout) :: rmatrix
  !</inputoutput>

!</subroutine>
 
   ! Local variables
  ! Determine whether
  !   to write the final solution in a data file
  !   to calculate the flow around cylinder parameters
  !   to export GMV/VTK outputs
  !   to calculate the Kinetic energy
  !   to write the real/projected data
  integer :: detWriteResult, LiftDragASO, ExporType
  integer :: KEnergy, detKEnergy, Ensto, detEnsto
  integer :: Vtild, Ptild, Wtild
  
  ! Kinematic viscosity noo = 1/Re  
  real(DP) :: dnu
     
   ! Path to the data file which has the initial solution
  character(LEN=SYS_STRLEN) :: sfile
 
  ! Output block for UCD output to GMV/VTK file
  type(t_ucdExport) :: rexport
  character(len=SYS_STRLEN) :: sucddir
  real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
  
  ! An object specifying the discretisation.
  type(t_blockDiscretisation) :: rprjDiscretisation,rprjDiscretisationn
  
  ! A block vector which contains projected data
  type(t_vectorBlock) :: rprjVector,rprjVectorr

  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_boundaryRegion) :: rboundaryRegion
  type(t_discreteBC), target :: rprjDiscreteBC
  
  ! Forces on the objects
  real(DP), dimension(2) :: Dforces, Dvalues   

  ! A list of points where to evaluate FEM data.
  ! DIMENSION(1..ndim,1..npoints)
  real(DP), dimension(2,2) :: Dpoints

  ! The 1x1 block mass matrix, created for the kinetic energy calculations.
  type(t_matrixBlock) :: rmass_matrix
  
  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation) :: rblockDiscr  
  
  ! The block vector, created for the kinetic energy calculations.
  type(t_vectorBlock) :: rxvel_vector, ryvel_vector, rvort_vector, ru1, ru2
  
  ! Kinetic energy and its sub values
  real(DP) :: dE, dU1, dU2

  ! Global Mass Conservation parameters
  integer :: detGMC
  real(DP) :: Dfluxi, Dfluxo, Dgmc
  real(DP), dimension(2,2) :: Dcoords
  character(len=SYS_STRLEN) :: sstring
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Write the final result in a data file. This can be later read as an
  !  initial solution for the non-linear loop.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! Determine whether to write the final solution in a data file
  call parlst_getvalue_int (rparams, 'POST', 'detWriteResult', detWriteResult, 0)   
  if (detWriteResult .eq. 1) then
  
    ! Write the final solution to a data file
    call parlst_getvalue_string (rparams, 'ISOLUTION', &
         'sFilenamePathResult',sfile, "", bdequote=.true.)
      
    call vecio_writeBlockVectorHR (rvector, 'SOLUTION', .true.,&
                       0, sfile, '(E22.15)')
  end if

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate drag-/lift coefficients on the 2nd boundary component.
  ! This is for the benchmark problem: flow around cylinder!
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

  ! Determine whether to calculate flow around cylinder parameters or not
  call parlst_getvalue_int (rparams, 'POST', 'LiftDragASO', LiftDragASO, 0)
  
  if (LiftDragASO .eq. 1) then
    
    call boundary_createRegion (rboundary,2,0, rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! Kinematic viscosity, noo = 1/Re
    call parlst_getvalue_double (rparams, 'GFPROPER', 'dnu', dnu, 1.0_DP)  
         
    call ppns2D_bdforces_uniform (rvector,rboundaryRegion,Dforces,CUB_G3_1D,&
        dnu,2.0_DP)

    call output_lbrk()
    call output_line ('Body forces')
    call output_line ('-----------')
    call output_line ('Drag/Lift')
    call output_line (trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
              //trim(sys_sdEP(Dforces(2),15,6)))

    call output_lbrk()
    call output_line ('Coefficients (Line Integration)')
    call output_line ('------------------------------')
    call output_line ('Drag/Lift')
    call output_line (trim(sys_sdEP(Dforces(1)*500.0_DP,15,6)) // ' / '&
              //trim(sys_sdEP(Dforces(2)*500.0_DP,15,6)))


    call ppns2D_bdforces_uniform (rvector,rboundaryRegion,Dforces,CUB_G3X3,&
        dnu,2.0_DP)
    !CUB_GEN_AUTO_G3
    call output_lbrk()
    call output_line ('Coefficients (Volume Integration)')
    call output_line ('--------------------------------')
    call output_line ('Drag/Lift')
    call output_line (trim(sys_sdEP(Dforces(1)*500.0_DP,15,6)) // ' / '&
      //trim(sys_sdEP(Dforces(2)*500.0_DP,15,6)))

    ! Calculate the pressure drop accross the cylinder
    Dpoints(1,1) = 0.15_DP
    Dpoints(2,1) = 0.2_DP
    Dpoints(1,2) = 0.25_DP
    Dpoints(2,2) = 0.2_DP
    call fevl_evaluate (DER_FUNC, Dvalues, rvector%RvectorBlock(3), Dpoints)
    
    call output_lbrk()
    call output_line ('Pressure Drop')
    call output_line ('------------')
    call output_line (trim(sys_sdEP(Dvalues(1)-Dvalues(2),15,6)))

  end if


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Writing the solution to GMV/VTK files.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! We can now start the exporting the results.
  ! Get the path for writing postprocessing files from the environment variable
  ! $UCDDIR. If that does not exist, write to the directory "./gmv".
  if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

  ! Determine which typeof export do we use, GMV/VTK
  call parlst_getvalue_int (rparams, 'POST', 'ExporType', ExporType, 0)

  ! Detemine whether we need to project the solution to a GMV/VTK readable
  ! type.
  call parlst_getvalue_int (rparams, 'MESH', 'Vtild', Vtild, 0)
  call parlst_getvalue_int (rparams, 'MESH', 'Ptild', Ptild, 0)
  call parlst_getvalue_int (rparams, 'MESH', 'Wtild', Wtild, 0)

  if ( (Vtild .eq. 1) .or. (Ptild .eq. 1) .or. (Wtild .eq. 1)) then
      
    ! make a new discretization structure for the projected data
    ! and modify its sub-structures if required
    call spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
    
    if (Vtild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
      EL_Q1, CUB_G2_2D, rprjDiscretisation%RspatialDiscr(1))

      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
      EL_Q1, CUB_G2_2D, rprjDiscretisation%RspatialDiscr(2))
    endif

    if (Ptild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(3), &
      EL_Q1, CUB_G2_2D, rprjDiscretisation%RspatialDiscr(3))         
    endif    

    if (Wtild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(4), &
      EL_Q1, CUB_G2_2D, rprjDiscretisation%RspatialDiscr(4))         
    endif 
   
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)

    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)

    !!!! THE SECOND PHASE !!!!

    ! make a new discretization structure for the projected data
    ! and modify its sub-structures if required
    call spdiscr_duplicateBlockDiscr (rprjDiscretisation,rprjDiscretisationn)
    
    if (Vtild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rprjDiscretisation%RspatialDiscr(1), &
      EL_Q2, CUB_G2_2D, rprjDiscretisationn%RspatialDiscr(1))

      call spdiscr_deriveSimpleDiscrSc (rprjDiscretisation%RspatialDiscr(2), &
      EL_Q2, CUB_G2_2D, rprjDiscretisationn%RspatialDiscr(2))
    endif

    if (Ptild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rprjDiscretisation%RspatialDiscr(3), &
      EL_Q1, CUB_G2_2D, rprjDiscretisationn%RspatialDiscr(3))         
    endif    

    if (Wtild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rprjDiscretisation%RspatialDiscr(4), &
      EL_Q1, CUB_G2_2D, rprjDiscretisationn%RspatialDiscr(4))         
    endif 
   
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisationn,rprjVectorr,.false.)

    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rprjVector,rprjVectorr)

    ! Discretise the boundary conditions according to the new discretisation
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rprjDiscreteBC)

    ! Prepare the discrete BCs. data
    call ls_BCs_onetime(rprjDiscretisationn,rboundary,rprjDiscreteBC)

    ! Hang the pointer into the vector.
    rprjVectorr%p_rdiscreteBC => rprjDiscreteBC

    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVectorr) 
    
    if (ExporType .eq. 0) then
    
      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nflsfem.vtk')

      ! Write Pressure
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

      ! Write velocity field
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

      ! Write Vorticity
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)
    
    else
      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nflsfem.gmv')  

      ! Write Pressure
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
      
      ! Write velocity field     
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
      
      ! Write Vorticity
      call lsyssc_getbase_double (rprjVectorr%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)     
    
    end if
    
    ! Release the temporary projected vector
    call lsysbl_releaseVector (rprjVector)
    call lsysbl_releaseVector (rprjVectorr)
    ! Release our discrete version of the projected boundary conditions
    call bcasm_releaseDiscreteBC (rprjDiscreteBC)
    ! Release the projected discretisation structure and 
    ! all spatial discretisation structures in it.
    call spdiscr_releaseBlockDiscr(rprjDiscretisation)    
    call spdiscr_releaseBlockDiscr(rprjDiscretisationn)
  
  else ! real data will be used
  
    if (ExporType .eq. 0) then
    
      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nflsfem.vtk')

      ! Write Pressure
      call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

      ! Write velocity field
      call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

      ! Write Vorticity
      call lsyssc_getbase_double (rvector%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)
    
    else
      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nflsfem.gmv')  

      ! Write Pressure
      call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
      
      ! Write velocity field     
      call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
      
      ! Write Vorticity
      call lsyssc_getbase_double (rvector%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)     
    
    end if
  
  end if ! end of real or projected data condition
    
  ! Write the file to disc, that is it.
  call ucd_write (rexport)
  call ucd_release (rexport)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculating the Kinetic Energy.
  !   E = 1/2 \int{u^2}
  !   using the definition of the velocities based on FEM we end up with:
  !   E = 1/2*[u^T][M][u] = 1/2*([u1^T][M11][u1] + [u2^T][M22][u2])
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  ! Determine whether to calculate Kinetic energy
  call parlst_getvalue_int (rparams, 'POST', 'KEnergy', KEnergy, 0)  
    
  if (KEnergy .eq. 1) then
  
    ! Determine how to calculate Kinetic energy
    call parlst_getvalue_int (rparams, 'POST', 'detKEnergy', detKEnergy, 0) 
    select case (detKEnergy)

    case (0)
      ! detKEnergy =0  use the FEM definition
      ! E = 1/2 [u]^T[M][u]
      
      ! Build up a one block mass matrix.
      ! Get the structure and data from the system matrix A11.
      ! BUT, first clear the system matrix data.
      call lsysbl_clearMatrix (rmatrix)
      
      ! Create the mass matrix discretization structure
      call spdiscr_createBlockDiscrInd (&
               rmatrix%p_rblockDiscrTrial%rspatialDiscr(1),&
                        rblockDiscr)  
      call lsysbl_createMatFromScalar (rmatrix%RmatrixBlock(1,1),&
                  rmass_matrix, rblockDiscr,rblockDiscr)
                  
      ! Bulid the mass matrix
      call bma_buildMatrix (rmass_matrix,BMA_CALC_STANDARD,ls_Mass,&
                      rcubatureInfo=rcubatureInfo)
       
      ! Extract the first block of the solution matrix,
      ! X-velocity block
      call lsysbl_createVecFromScalar(rvector%RvectorBlock(1),&
                       rxvel_vector, rblockDiscr)
      ! Create a temporary vector  
      call lsysbl_createVecBlockIndirect (rxvel_vector,ru1,.true.)
      
      ! Do the matrix-vector multiplication
      ! ru1   =   cx * rmass_matrix * rxvel_vector   +   cy * ru1
      call lsysbl_blockMatVec(rmass_matrix, rxvel_vector,&
                     ru1, cx=1.0_DP, cy=0.0_DP)
      
      ! Do the vector-vector multiplication
      dU1 = lsysbl_scalarProduct(ru1,rxvel_vector)
      
      ! Extract the second block of the solution matrix,
      ! Y-velocity block
      call lsysbl_createVecFromScalar(rvector%RvectorBlock(2),&
                        ryvel_vector, rblockDiscr)
      ! Create a temporary vector
      call lsysbl_createVecBlockIndirect (ryvel_vector,ru2,.true.)
      ! Do the matrix-vector multiplication
      ! ru2   =   cx * rmass_matrix * ryvel_vector   +   cy * ru2
      call lsysbl_blockMatVec(rmass_matrix, ryvel_vector, ru2,&
                           cx=1.0_DP, cy=0.0_DP)
      
      ! Do the vector-vector multiplication
      dU2 = lsysbl_scalarProduct(ru2,ryvel_vector)
      
      ! Kinetic energy
      dE = 0.5_DP*(dU1 + dU2)
      
      ! Print the Kinetic energy value
      call output_lbrk()
      call output_line ('Kinetic energy - based on mass matrix')
      call output_line ('-------------------------------------')
      call output_line (trim(sys_sdEP(dE,15,6)))  

      ! Release the discretisation structure
       !and all spatial discretisation structures in it.
      call spdiscr_releaseBlockDiscr(rblockDiscr)
        
      ! Release the temporary vectors and matrix  
      call lsysbl_releaseVector (ru1)
      call lsysbl_releaseVector (rxvel_vector)
      call lsysbl_releaseVector (ru2)
      call lsysbl_releaseVector (ryvel_vector)
      call lsysbl_releaseMatrix (rmass_matrix)  
        
    case (1)
      ! detKEnergy = 1  simply take the L^2 norm of velocity vectors
      !  E = 1/2||u||^2_{L^2} 
          
      ! Call the error analysis subroutine without an analytical function
      ! to calculate the L^2 norms: ||u1||_{L^2} , ||u2||_{L^2}
      call pperr_scalar (PPERR_L2ERROR,dU1,rvector%RvectorBlock(1),&
        rcubatureInfo=rcubatureInfo)
      call pperr_scalar (PPERR_L2ERROR,dU2,rvector%RvectorBlock(2),&
        rcubatureInfo=rcubatureInfo)
      
      ! Kinetic energy       
      dE = 0.5_DP*(dU1**2 + dU2**2)

      ! Print the Kinetic energy value
      call output_lbrk()
      call output_line ('Kinetic energy - based on L^2 norms')
      call output_line ('-----------------------------------')
      call output_line (trim(sys_sdEP(dE,15,6)))  
       
    end select
    
  end if 


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculating the Enstophy.
  !   Z = 1/2 \int{w^2}  see: http://en.wikipedia.org/wiki/Enstrophy
  !   using the definition of the vorticity based on FEM we end up with:
  !   Z = 1/2*[w^T][M][w]
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  ! Determine whether to calculate flow around cylinder parameters or not
  call parlst_getvalue_int (rparams, 'POST', 'Ensto', Ensto, 0)  
  
  if (Ensto .eq. 1) then

    ! Determine how to calculate Kinetic energy
    call parlst_getvalue_int (rparams, 'POST', 'detEnsto', detEnsto, 0) 
    select case (detEnsto)

    case (0)
      ! detKEnergy =0  use the FEM definition
      ! E = 1/2 [w]^T[M][w]
      
      ! Build up a one block mass matrix. Get the structure and data from
      ! the system matrix A11. BUT, first clear the system matrix data.
      call lsysbl_clearMatrix (rmatrix)
      
      ! Create the mass matrix discretization structure
      call spdiscr_createBlockDiscrInd (&
           rmatrix%p_rblockDiscrTrial%rspatialDiscr(4),rblockDiscr)  
      call lsysbl_createMatFromScalar (rmatrix%RmatrixBlock(4,4),&
                   rmass_matrix, rblockDiscr,rblockDiscr)
                  
      ! Bulid the mass matrix
      call bma_buildMatrix (rmass_matrix,BMA_CALC_STANDARD,ls_Mass,&
                      rcubatureInfo=rcubatureInfo)
       
      ! Extract the first block of the solution matrix, X-velocity block
      call lsysbl_createVecFromScalar(rvector%RvectorBlock(4),&
                       rvort_vector, rblockDiscr)
      ! Create a temporary vector  
      call lsysbl_createVecBlockIndirect (rvort_vector,ru1,.true.)
      
      ! Do the matrix-vector multiplication
      ! ru1   =   cx * rmass_matrix * rvort_vector   +   cy * ru1
      call lsysbl_blockMatVec(rmass_matrix, rvort_vector, &
                       ru1, cx=1.0_DP, cy=0.0_DP)
      
      ! Do the vector-vector multiplication
      dU1 = lsysbl_scalarProduct(ru1,rvort_vector)
      
      ! Enstrophy
      dE = 0.5_DP*dU1
      
      ! Print the Enstrophy value
      call output_lbrk()
      call output_line ('Enstophy - based on mass matrix')
      call output_line ('-------------------------------')
      call output_line (trim(sys_sdEP(dE,15,6)))  

      ! Release the discretisation structure and all 
      ! spatial discretisation structures in it.
      call spdiscr_releaseBlockDiscr(rblockDiscr)
        
      ! Release the temporary vectors and matrix  
      call lsysbl_releaseVector (ru1)
      call lsysbl_releaseVector (rvort_vector)
      call lsysbl_releaseMatrix (rmass_matrix)  
        
    case (1)
      ! detEnsto = 1  simply take the L^2 norm of vorticity
      !  Z = 1/2||w||^2_{L^2} 
          
      ! Call the error analysis subroutine without an analytical function
      ! to calculate the L^2 norms: ||w||_{L^2}
      call pperr_scalar (PPERR_L2ERROR,dU1,rvector%RvectorBlock(4),&
                         rcubatureInfo=rcubatureInfo)
      
      ! Kinetic energy       
      dE = 0.5_DP*(dU1**2)

      ! Print the Enstrophy value
      call output_lbrk()
      call output_line ('Enstophy - based on L^2 norms')
      call output_line ('---------')
      call output_line (trim(sys_sdEP(dE,15,6)))   

    end select

  end if 


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate Global Mass Conservation (GMC)
  !  This applies to the channel flows ONLY. The GMC is the normalised 
  !  difference between the input and output velocity fluxes (mass flow rate)
  !  of the domain. It is defined as:
  !
  !		  \int_{\Gamma_i}{n.v} - \int_{\Gamma_o}{n.v}
  !  GMC = --------------------------------------------- * 100
  !				   \int_{\Gamma_i}{n.v}
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! Check whether to calculate the GMC or not.
  call parlst_getvalue_int (rparams, 'POST', 'detGMC',detGMC,0)
  
  if (detGMC .eq. 1) then
    ! Input line coordinates
    call parlst_getvalue_string (rparams, 'POST', 'inputGMC',sstring)
    read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
             Dcoords(1,2), Dcoords(2,2)
    ! Input line flux
    call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                         Dcoords(1:2,2),Dfluxi)                               

    Dfluxi = abs(Dfluxi)
    
    ! Output line coordinates
    call parlst_getvalue_string (rparams, 'POST', 'outputGMC',sstring)
    read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
             Dcoords(1,2), Dcoords(2,2)
    ! Output line flux
    call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                         Dcoords(1:2,2),Dfluxo)
    
    Dfluxo = abs(Dfluxo)
    
    ! The GMC is then calculated as
    Dgmc = 100.0_DP*(Dfluxi - Dfluxo) / Dfluxi                           
    
    ! Print the GMC value
    call output_lbrk()
    call output_line ('Global Mass Conservation(%)')
    call output_line ('---------------------------')
    call output_line (trim(sys_sdEP(Dgmc,15,6)))     
    
  end if
     
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_cleanup(rvector,rvector_old,rrhs,rmatrix,rcubatureInfo,&
             rdiscreteBC,rdiscretisation,rtriangulation,rboundary)
 
 !<description>  
  ! Release all the memory used in our calculations.
 !</description> 

  !<inputoutput>
  ! An object for saving the domain:
  type(t_boundary) :: rboundary

  ! An object for saving the triangulation on the domain
  type(t_triangulation) :: rtriangulation

  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation) :: rdiscretisation
  
  ! A block matrix and a couple of block vectors. These will be filled
  ! with data for the linear solver.
  type(t_matrixBlock) :: rmatrix
  type(t_vectorBlock) :: rvector,rvector_old,rrhs
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo) :: rcubatureInfo
  
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC), target :: rdiscreteBC
  !</inputoutput>

!</subroutine>

  ! Now, clean up so that all the memory is available again.
    
  ! Release the block matrix/vectors
  call lsysbl_releaseVector (rvector)
  call lsysbl_releaseVector (rvector_old)
  call lsysbl_releaseVector (rrhs)
  call lsysbl_releaseMatrix (rmatrix)
  
  ! Release the cubature formula
  call spdiscr_releaseCubStructure (rcubatureInfo)
  
  ! Release our discrete version of the boundary conditions
  call bcasm_releaseDiscreteBC (rdiscreteBC)

  ! Release the discretisation structure and all spatial discretisation
  ! structures in it.
  call spdiscr_releaseBlockDiscr(rdiscretisation)
    
  ! Release the triangulation.
  call tria_done (rtriangulation)
  
  ! Finally release the domain, that is it.
  call boundary_release (rboundary)

  end subroutine
 
 
  ! ***************************************************************************

!<subroutine>
  subroutine ls_ns2D_Matrix(RmatrixData,rassemblyData,rmatrixAssembly,&
             npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
!</description>

!<inputoutput>
  ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
  ! have to be filled with data.
  type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
  
  ! Number of points per element
  integer, intent(in) :: npointsPerElement
  
  ! Number of elements
  integer, intent(in) :: nelements

  ! Values of FEM functions automatically evaluated in the
  ! cubature points.
  type(t_fev2Vectors), intent(in) :: revalVectors

  ! User defined collection structure
  type(t_collection), intent(inout), target, optional :: rcollection
!</input>
  
!</subroutine>

  ! Local variables
  real(DP) :: dbasI, dbasJ, dval, dbasIx, dbasIy, dbasJx, dbasJy, dnu
  integer :: iel, icubp, idofe, jdofe
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA12
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA13,p_DlocalMatrixA14
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA21,p_DlocalMatrixA22
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA23,p_DlocalMatrixA24
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA33,p_DlocalMatrixA34  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA41,p_DlocalMatrixA42
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA43,p_DlocalMatrixA44  
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA33,p_DbasTestA33
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA44,p_DbasTestA44
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA14,p_DbasTestA14
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA34,p_DbasTestA34

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11,p_rmatrixDataA33,p_rmatrixDataA44
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA13
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA14,p_rmatrixDataA34
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

  ! Element volumes used in ADN scaling
  real(DP), dimension(:), pointer :: p_DelementVolume
  
  ! Scaling factors
  real(DP), dimension(nelements) :: Dadn
  real(DP) :: Dphy, Dfpn
  
  ! Continuity equation scaling
  real(DP) :: alpha
    
  ! A handle to the element volumes
  integer :: ihandle, i, ielreal
  
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dUx, dUy, dVx, dVy
    
  ! Viscosity
  dnu = rcollection%DquickAccess(1)
  
  ! Continuity equation scaling factor
  !  default = 1.0_DP  no scaling
  alpha = rcollection%DquickAccess(3)

  ! Linearization Scheme
  Dfpn = rcollection%DquickAccess(4)
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
  p_rmatrixDataA33 => RmatrixData(3,3)
  p_rmatrixDataA44 => RmatrixData(4,4)
  p_rmatrixDataA13 => RmatrixData(1,3)
  p_rmatrixDataA14 => RmatrixData(1,4)
  p_rmatrixDataA34 => RmatrixData(3,4)
  
  
  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
  p_DbasTrialA33 => RmatrixData(3,3)%p_DbasTrial
  p_DbasTestA33 => RmatrixData(3,3)%p_DbasTest
  p_DbasTrialA44 => RmatrixData(4,4)%p_DbasTrial
  p_DbasTestA44 => RmatrixData(4,4)%p_DbasTest 
  p_DbasTrialA13 => RmatrixData(1,3)%p_DbasTrial
  p_DbasTestA13 => RmatrixData(1,3)%p_DbasTest
  p_DbasTrialA14 => RmatrixData(1,4)%p_DbasTrial
  p_DbasTestA14 => RmatrixData(1,4)%p_DbasTest
  p_DbasTrialA34 => RmatrixData(3,4)%p_DbasTrial
  p_DbasTestA34 => RmatrixData(3,4)%p_DbasTest   
  
  
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
  p_DlocalMatrixA12 => RmatrixData(1,2)%p_Dentry
  p_DlocalMatrixA13 => RmatrixData(1,3)%p_Dentry
  p_DlocalMatrixA14 => RmatrixData(1,4)%p_Dentry
  p_DlocalMatrixA21 => RmatrixData(2,1)%p_Dentry
  p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
  p_DlocalMatrixA23 => RmatrixData(2,3)%p_Dentry
  p_DlocalMatrixA24 => RmatrixData(2,4)%p_Dentry
  p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
  p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
  p_DlocalMatrixA33 => RmatrixData(3,3)%p_Dentry
  p_DlocalMatrixA34 => RmatrixData(3,4)%p_Dentry
  p_DlocalMatrixA41 => RmatrixData(4,1)%p_Dentry
  p_DlocalMatrixA42 => RmatrixData(4,2)%p_Dentry
  p_DlocalMatrixA43 => RmatrixData(4,3)%p_Dentry
  p_DlocalMatrixA44 => RmatrixData(4,4)%p_Dentry   
  
  
  ! Get the velocity field from the parameters
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata  

  ! Get the element volume to be used in ADN scaling
  ! Check if ADN scaling is required
  if (rcollection%IquickAccess(1) .eq. 1) then
    ihandle = rmatrixAssembly%P_RTRIANGULATION%H_DELEMENTVOLUME
    call storage_getbase_double (ihandle, p_DelementVolume)
    do i = 1,nelements
      ielreal = rassemblyData%P_IELEMENTLIST(i)
      Dadn(i) = 1.0_DP/p_DelementVolume(ielreal)
    end do
  else
    Dadn = 1.0_DP  
  end if

  ! Assigne the value of the Physical scaling
  ! This value is set in subroutine 'ls_initParams' 
  ! If no scaling is required, this value is set to 1.0_DP
  Dphy = rcollection%DquickAccess(2)
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! We do the calculation in a block-by-block manner. All the
  ! relevant blocks will be calculated in the same loop over the
  ! elements.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  ! +++++++++++++++++++++++++++++++++++
  ! Calculate blocks A11, A22, A12, A21 
  ! +++++++++++++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Velocity/derivatives field in this cubature point
    dU = p_Du1(icubp,iel,DER_FUNC)
    dV = p_Du2(icubp,iel,DER_FUNC)
    
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA11%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA11(idofe,DER_FUNC,icubp,iel)
      dbasIx = p_DbasTestA11(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA11(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA11%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA11(jdofe,DER_FUNC,icubp,iel)
      dbasJx = p_DbasTrialA11(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA11(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A11
      dval = p_DcubWeight(icubp,iel) * (  Dphy*(dU*dU * dbasJx*dbasIx + &
           dU*dV * dbasJx*dbasIy + dV*dU * dbasJy*dbasIx + &
           dV*dV * dbasJy*dbasIy + Dfpn*(dU*dUx * dbasJx*dbasI + &
           dV*dUx * dbasJy*dbasI + dUx*dU * dbasJ*dbasIx + &
           dUx*dV * dbasJ*dbasIy + dUx*dUx * dbasJ*dbasI + &
           dVx*dVx * dbasJ*dbasI)) + Dadn(iel)*(dbasJy*dbasIy + alpha*dbasJx*dbasIx)  )
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
      ! A22
      dval = p_DcubWeight(icubp,iel) * (  Dphy*(dU*dU * dbasJx*dbasIx + &
           dU*dV * dbasJx*dbasIy + dV*dU * dbasJy*dbasIx + &
           dV*dV * dbasJy*dbasIy + Dfpn*(dU*dVy * dbasJx*dbasI + &
           dV*dVy * dbasJy*dbasI + dVy*dU * dbasJ*dbasIx + &
           dVy*dV * dbasJ*dbasIy + dUy*dUy * dbasJ*dbasI + &
           dVy*dVy * dbasJ*dbasI)) + Dadn(iel)*(alpha*dbasJy*dbasIy + dbasJx*dbasIx)  )      
      p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
      
      ! A12
      dval = p_DcubWeight(icubp,iel) * (  Dfpn*Dphy*(dU*dVx * dbasJx*dbasI + &
           dV*dVx * dbasJy*dbasI + dUy*dU * dbasJ*dbasIx + &
           dUy*dV * dbasJ*dbasIy + dUy*dUx * dbasJ*dbasI + &
           dVy*dVx * dbasJ*dbasI) + Dadn(iel)*(-dbasJx*dbasIy + alpha*dbasJy*dbasIx)  )      
      p_DlocalMatrixA12(jdofe,idofe,iel) = p_DlocalMatrixA12(jdofe,idofe,iel) + dval
      
      ! A21
      p_DlocalMatrixA21(idofe,jdofe,iel) = p_DlocalMatrixA21(idofe,jdofe,iel) + dval      
      
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  

  ! +++++++++++++++++++++++++++++++++++
  ! Calculate blocks A13, A23, A31, A32 
  ! +++++++++++++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Velocity/derivatives field in this cubature point
    dU = p_Du1(icubp,iel,DER_FUNC)
    dV = p_Du2(icubp,iel,DER_FUNC)
    
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA13%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA13(idofe,DER_FUNC,icubp,iel)
      dbasIx = p_DbasTestA13(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA13(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA13%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA13(jdofe,DER_FUNC,icubp,iel)
      dbasJx = p_DbasTrialA13(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA13(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A13
      dval = p_DcubWeight(icubp,iel) * Dphy*(  dU * dbasJx*dbasIx + &
           dV * dbasJx*dbasIy + Dfpn*(dUx * dbasJx*dbasI + dVx * dbasJy*dbasI)  )
      p_DlocalMatrixA13(jdofe,idofe,iel) = p_DlocalMatrixA13(jdofe,idofe,iel) + dval
      
      ! A31
      p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
           
      ! A23
      dval = p_DcubWeight(icubp,iel) * Dphy*(  dU * dbasJy*dbasIx + &
           dV * dbasJy*dbasIy + Dfpn*(dUy * dbasJx*dbasI + dVy * dbasJy*dbasI)  )
      p_DlocalMatrixA23(jdofe,idofe,iel) = p_DlocalMatrixA23(jdofe,idofe,iel) + dval
      
      ! A32
      p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + dval
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  
  ! +++++++++++++++++++++++++++++++++++
  ! Calculate blocks A14, A24, A41, A42 
  ! +++++++++++++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Velocity/derivatives field in this cubature point
    dU = p_Du1(icubp,iel,DER_FUNC)
    dV = p_Du2(icubp,iel,DER_FUNC)
    
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA14%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA14(idofe,DER_FUNC,icubp,iel)
      dbasIx = p_DbasTestA14(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA14(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA14%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA14(jdofe,DER_FUNC,icubp,iel)
      dbasJx = p_DbasTrialA14(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA14(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A14
      dval = p_DcubWeight(icubp,iel) * (  dnu*Dphy*(  dU * dbasJy*dbasIx + &
           dV * dbasJy*dbasIy + Dfpn*(dUx * dbasJy*dbasI - dVx * dbasJx*dbasI)  ) + &
           Dadn(iel)*(dbasJ*dbasIy)  )
      p_DlocalMatrixA14(jdofe,idofe,iel) = p_DlocalMatrixA14(jdofe,idofe,iel) + dval

      ! A41
      p_DlocalMatrixA41(idofe,jdofe,iel) = p_DlocalMatrixA41(idofe,jdofe,iel) + dval

      ! A24
      dval = p_DcubWeight(icubp,iel) * (  dnu*Dphy*(  -dU * dbasJx*dbasIx - &
           dV * dbasJx*dbasIy + Dfpn*(dUy * dbasJy*dbasI - dVy * dbasJx*dbasI)  ) - &
           Dadn(iel)*(dbasJ*dbasIx)  )
      p_DlocalMatrixA24(jdofe,idofe,iel) = p_DlocalMatrixA24(jdofe,idofe,iel) + dval
      
      ! A42
      p_DlocalMatrixA42(idofe,jdofe,iel) = p_DlocalMatrixA42(idofe,jdofe,iel) + dval
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  
  
  ! +++++++++++++++++++
  ! Calculate block A33 
  ! +++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA33%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA33(idofe,DER_FUNC,icubp,iel)
      dbasIx = p_DbasTestA33(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA33(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA33%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA33(jdofe,DER_FUNC,icubp,iel)
      dbasJx = p_DbasTrialA33(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA33(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A33
      dval = p_DcubWeight(icubp,iel) * Dphy*(dbasJx*dbasIx + dbasJy*dbasIy)
      p_DlocalMatrixA33(jdofe,idofe,iel) = p_DlocalMatrixA33(jdofe,idofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  

  ! +++++++++++++++++++++++++
  ! Calculate blocks A34, A43 
  ! +++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA34%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasIx = p_DbasTestA34(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA34(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA34%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJx = p_DbasTrialA34(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA34(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A34
      dval = dnu * Dphy* p_DcubWeight(icubp,iel) * ( dbasJy*dbasIx - dbasJx*dbasIy )
      p_DlocalMatrixA34(jdofe,idofe,iel) = p_DlocalMatrixA34(jdofe,idofe,iel) + dval

      ! A43
      p_DlocalMatrixA43(idofe,jdofe,iel) = p_DlocalMatrixA43(idofe,jdofe,iel) + dval

                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  
  
  ! +++++++++++++++++++
  ! Calculate block A44 
  ! +++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA44%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA44(idofe,DER_FUNC,icubp,iel)      
      dbasIx = p_DbasTestA44(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA44(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA44%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA44(jdofe,DER_FUNC,icubp,iel)      
      dbasJx = p_DbasTrialA44(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA44(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A44
      dval = p_DcubWeight(icubp,iel) * ( dnu*dnu*Dphy*(dbasJy*dbasIy + dbasJx*dbasIx) +&
           Dadn(iel)*(dbasJ*dbasI) )
      p_DlocalMatrixA44(jdofe,idofe,iel) = p_DlocalMatrixA44(jdofe,idofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  
  end subroutine

  !****************************************************************************


!<subroutine>
  subroutine ls_ns2D_rhs(rvectorData,rassemblyData,rvectorAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the RHS vector in a block-by-block procedures.
!</description>

!<inputoutput>
  ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
  ! have to be filled with data.
  type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData
!</inputoutput>

!<input>
  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
  
  ! Number of points per element
  integer, intent(in) :: npointsPerElement
  
  ! Number of elements
  integer, intent(in) :: nelements
  
  ! Values of FEM functions automatically evaluated in the
  ! cubature points.
  type(t_fev2Vectors), intent(in) :: revalVectors

  ! User defined collection structure
  type(t_collection), intent(inout), target, optional :: rcollection
!</input>
  
!<subroutine>

  ! Local variables
  real(DP) :: dbasI,dbasIx,dbasIy, dval1, dval2, dval3, dval4, dnu
  integer :: iel, icubp, idofe
  real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
  real(DP), dimension(:,:), pointer :: p_DlocalVector3, p_DlocalVector4
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTest1,p_DbasTest3,p_DbasTest4
  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData3
  type(t_bmaVectorData), pointer :: p_rvectorData4


  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dUx, dUy, dVx, dVy

  real(DP) :: Dphy, Dfpn
    
  ! Viscosity
  dnu = rcollection%DquickAccess(1)

  ! Assigne the value of the Physical scaling
  ! This value is set in subroutine 'ls_initParams' 
  ! If no scaling is required, this value is set to 1.0_DP
  Dphy = rcollection%DquickAccess(2)
  
  ! Linearization Scheme
  Dfpn = rcollection%DquickAccess(4)
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rvectorData1 => RvectorData(1)
  p_rvectorData3 => RvectorData(3)
  p_rvectorData4 => RvectorData(4)

  p_DlocalVector1 => RvectorData(1)%p_Dentry
  p_DlocalVector2 => RvectorData(2)%p_Dentry
  p_DlocalVector3 => RvectorData(3)%p_Dentry
  p_DlocalVector4 => RvectorData(4)%p_Dentry

  p_DbasTest1 => RvectorData(1)%p_DbasTest
  p_DbasTest3 => RvectorData(3)%p_DbasTest
  p_DbasTest4 => RvectorData(4)%p_DbasTest  
  
  
  ! Get the velocity field from the parameters
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata  
    
  ! Calculate the RHS of the velocities
  
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Velocity/derivatives field in this cubature point
    dU = p_Du1(icubp,iel,DER_FUNC)
    dV = p_Du2(icubp,iel,DER_FUNC)
    
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rvectorData1%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      ! into dbasI
      dbasI = p_DbasTest1(idofe,DER_FUNC,icubp,iel)
      dbasIx = p_DbasTest1(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTest1(idofe,DER_DERIV2D_Y,icubp,iel)
                
      ! Values of the velocity RHS for the X1 and X2 component
      dval1 = 0.0_DP + Dfpn*(   (dU*dUx + dV*dUy) * dU * dbasIx + &
          (dU*dUx + dV*dUy) * dV * dbasIy + &
          (dU*dUx + dV*dUy) * dUx * dbasI + &
          (dU*dVx + dV*dVy) * dVx * dbasI   )
          
      dval2 = 0.0_DP + Dfpn*(   (dU*dVx + dV*dVy) * dU * dbasIx + &
          (dU*dVx + dV*dVy) * dV * dbasIy + &
          (dU*dUx + dV*dUy) * dUy * dbasI + &
          (dU*dVx + dV*dVy) * dVy * dbasI   )
          
      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local vectors.
      p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
        p_DcubWeight(icubp,iel) * dval1 * Dphy
      p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
        p_DcubWeight(icubp,iel) * dval2 * Dphy
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  


  ! Calculate the RHS of the pressure
  
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Velocity/derivatives field in this cubature point
    dU = p_Du1(icubp,iel,DER_FUNC)
    dV = p_Du2(icubp,iel,DER_FUNC)
    
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rvectorData3%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      ! into dbasI
      dbasIx = p_DbasTest3(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTest3(idofe,DER_DERIV2D_Y,icubp,iel)
                
      ! Values of the velocity RHS for pressure
      dval3 = 0.0_DP + Dfpn*(  (dU*dUx + dV*dUy) * dbasIx + &
          (dU*dVx + dV*dVy) * dbasIy  )
          
      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local vectors.
      p_DlocalVector3(idofe,iel) = p_DlocalVector3(idofe,iel) + &
        p_DcubWeight(icubp,iel) * dval3 * Dphy
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel


   ! Calculate the RHS of the vorticity
  
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Velocity/derivatives field in this cubature point
    dU = p_Du1(icubp,iel,DER_FUNC)
    dV = p_Du2(icubp,iel,DER_FUNC)
    
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rvectorData4%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      ! into dbasI
      dbasIx = p_DbasTest4(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTest4(idofe,DER_DERIV2D_Y,icubp,iel)
                
      ! Values of the velocity RHS for vorticity
      dval4 = 0.0_DP + Dfpn*dnu * (  (dU*dUx + dV*dUy) * dbasIy - &
          (dU*dVx + dV*dVy) * dbasIx  )
          
      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local vectors.
      p_DlocalVector4(idofe,iel) = p_DlocalVector4(idofe,iel) + &
        p_DcubWeight(icubp,iel) * dval4 * Dphy
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
             npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble a mass matrix in a block-by-block procedures.
!</description>

!<inputoutput>
  ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
  ! have to be filled with data.
  type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
  
  ! Number of points per element
  integer, intent(in) :: npointsPerElement
  
  ! Number of elements
  integer, intent(in) :: nelements

  ! Values of FEM functions automatically evaluated in the
  ! cubature points.
  type(t_fev2Vectors), intent(in) :: revalVectors

  ! User defined collection structure
  type(t_collection), intent(inout), target, optional :: rcollection
!</input>
  
!<subroutine>

  ! Local variables
  real(DP) :: dbasI, dbasJ
  integer :: iel, icubp, idofe, jdofe
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11  
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11
  

  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
   
  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
  
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry  
    
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! We do the calculation in a block-by-block manner.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  ! ++++++++++++++++++++++++++++++++++
  ! Calculate blocks A11, mass matrix. 
  ! ++++++++++++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA11%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA11(idofe,DER_FUNC,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA11%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA11(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A11
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + &
                         p_DcubWeight(icubp,iel) * ( dbasJ*dbasI)
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  end subroutine

end module

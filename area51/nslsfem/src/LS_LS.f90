!##############################################################################
!# ****************************************************************************
!# <name> LS_LS </name>
!# ****************************************************************************
!# <purpose>
!#  This has been one of my dreams over the past 4-5 years to work on
!#  multiphase flows. Here, I could only implement the level-set equations and 
!#  run very simple examples :(
!#  The module is not mature enough to be coupled with Navier-Stokes yet!
!# </purpose>
!#
!# Author:        Masoud Nickaeen
!# First Version: Aug  30, 2012
!# Last Update:   Sep  06, 2012
!##############################################################################

module LS_LS

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
  use multilevelprojection
  use dofmapping
  use meshmodification
  use jumpstabilisation
  use matrixio
  use analyticprojection
  
  implicit none

contains
  
  !****************************************************************************

!<subroutine>
  subroutine ls_levelset
  
!<description>

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

  ! The || grad(\phi) - 1 ||_{L_2,T} vector
  ! type(t_vectorScalar) :: rindicator

  ! An object specifying the discretisation (structure of the
  ! solution, trial/test functions,...)
  type(t_blockDiscretisation) :: rdiscretisation

  ! The descretized BC object
  type(t_discreteBC) :: rdiscreteBC
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo) :: rcubatureInfo
  
  ! A system matrix.
  type(t_matrixBlock) :: rmatrix
  
  ! A projection matrix.
  type(t_matrixBlock) :: rmatrix_projection
    
  ! A couple of block vectors.
  type(t_vectorBlock) :: rvector,rvector_old,rvector_oldT,rrhs,rvectorU,rvectorV
    
  ! Max. Min. level to solve the problem
  integer :: NLMAX
    
  ! Loop index
  integer :: i

  ! Parameters used in time loop
  integer :: itime,NMaxTime
  real(DP) :: dTimeStep
    
  ! Parameters used in nonliner loop
  integer :: inl,NLN_Max
  
  ! Re-initialization parameters
  integer :: RI_Loop_MAX, RI_MAX, reini_freq, reiType
  
  ! BCs.
  integer :: detBCs
  
  ! Convergence parameter, either by error or by NLN_Max
  logical :: converged,diverged

  ! Nonlinear loop control
  logical :: det
    
  ! Collection structure for callback routines
  type(t_collection) :: rcollection
  
  ! All parameter of the LSFEM solver
  type(t_parlist) :: rparams
  
  ! Error value
  real(DP) :: derror
  
  ! Ok, let's start.
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 0)-Read all the parameters from data file and initialize the collection.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call ls_initParams(rparams,NMaxTime,dTimeStep,RI_MAX,reini_freq,reiType,&
                                                           detBCs,rcollection)
  

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! a)-Read the domain, the mesh and refine it. All about computational GRID.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_grid(rtriangulation,rboundary,rparams,rcollection,NLMAX) !,rindicator)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! b)- Set up a discretisation and cubature info structure which tells 
  ! the code which finite element to use.
  ! Also, create a 4*4 block matrix structure.
  ! Initialize the structure of the deferred velocities to be use in the
  !  nonlinear matrix assembly routine. This is done for all grid levels.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_structure(rtriangulation,rboundary,rparams,rdiscretisation,&
                                  rmatrix, rmatrix_projection,rcubatureInfo)

  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! c)- Initialization of the solution/RHS vectors.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize the discrete boundary conditions
  
  ! Initialize the RHS and fine level solution vectors
  call ls_Init_RhsndSolution(rdiscretisation,rmatrix_projection,rrhs, &
         rvector_old,rvector_oldT,rvector,rvectorU,rvectorV,rparams,&
         rcubatureInfo,rdiscreteBC,rboundary)

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The Nonlinear Loop, the outermost loop, starts here.
  ! The loop performs maximum 'NMaxTime' times.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do itime=1,NMaxTime
     
    call output_line ('-----------------------------------')
    call output_line ('Time step '//trim(sys_siL(itime,6))// &
                      '     Time '//trim(sys_sdL(itime*dTimeStep,5)))
    call output_lbrk ()
    
    det = .false.
    
    ! Real simulation time for time-dependent velocity fields
    rcollection%DquickAccess(3) = (itime-1)*dTimeStep
        
    ! ++++++++++++++++++++++++++++++
    ! Re-initialization loop control
    ! ++++++++++++++++++++++++++++++ 

    !Calculate the H_1 error = ||grad(\phi)||_{L_2}
    call pperr_scalar (PPERR_H1ERROR,derror,rvector_old%RvectorBlock(1),&
         rcubatureInfo=rcubatureInfo) ! ,relementError=rindicator
    call output_line ('||grad(phi)||_0   ' // sys_sdEL(derror,10) )
    call output_lbrk ()
    
    if (reiType == 0) then
    
      ! Adaptive re-initialization
      if ((derror .gt. 1.1_DP) .or. (derror .lt. 0.9_DP)) then
        call output_line ('Adaptive re-initialization invoked!')
        rcollection%IquickAccess(1) = 1
        RI_Loop_MAX = RI_MAX
        
        ! Copy the initial/previous solution to the current solution
        call lsysbl_copyVector (rvector_old,rvector_oldT)
      else
        rcollection%IquickAccess(1) = 0
        RI_Loop_MAX = 1
      end if        
        
    else
    
      ! Manual re-initialization
      if (mod(itime,reini_freq)== 0) then
        call output_line ('Manual re-initialization invoked!')
        rcollection%IquickAccess(1) = 1
        RI_Loop_MAX = RI_MAX
        
        ! Copy the initial/previous solution to the current solution
        call lsysbl_copyVector (rvector_old,rvector_oldT)
      else
        rcollection%IquickAccess(1) = 0
        RI_Loop_MAX = 1
      end if
      
    end if
            
    do i=1,RI_Loop_MAX

      ! +++++++++++++++++++++++++
      ! 1- System Matrix Assembly
      ! +++++++++++++++++++++++++
      call ls_MatAssembly(rmatrix,rvector_old,rvectorU,rvectorV,&
           rcollection,rparams,NLMAX,rcubatureInfo)
    
 
      ! +++++++++++++++
      ! 2- RHS Assembly
      ! +++++++++++++++
      call ls_RHS_Assembly(rrhs,rvector_old,rvector_oldT,rvectorU,rvectorV,&
           rcollection,rcubatureInfo) 
        
      ! ++++++
      ! 3- BCs
      ! ++++++
      if (detBCs == 1) call ls_BCs_Dirichlet(rmatrix,rrhs,rvector,&
                                rvector_old,rvector_oldT,rdiscreteBC)
      
      
      ! +++++++++++++++++++
      ! 4- Solve the System
      ! +++++++++++++++++++
      call ls_Solver(rmatrix,rrhs,rvector,rparams)
          
      
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 6- Update the Initial Solution, clear matrix, vectors
      !  or Terminate the Loop
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
      call ls_update_solution(rmatrix,rvector,rvector_old,rrhs,rcollection,&
                                                                i,RI_Loop_MAX)
  
  
    end do ! i
    
    ! Write the time-dependent data
    call ls_postprocess(rboundary,rmatrix,rvector,rtriangulation,&
     rcubatureInfo,rdiscretisation,rparams, itime)
  
  end do  ! itime
    
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! f)- Clean up, free all the memory which is used for variables.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_cleanup(rvector,rvector_old,rvector_oldT,rvectorU,rvectorV,rrhs,&
  rmatrix,rmatrix_projection,rboundary,rparams,rdiscretisation,&
  rtriangulation,rcubatureInfo,rdiscreteBC)

  end subroutine


  !****************************************************************************
  
!<subroutine>
  subroutine ls_initParams(rparams,NMaxTime,dTimeStep,RI_MAX,reini_freq,&
                                                   reiType,detBCs,rcollection)

 !<description>  
  ! In this subroutine, the collection structrure is initialized.
 !</description>                

 !<output>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(out) :: rparams

  ! Time loop Max. number of iterations
  integer, intent(out) :: NMaxTime
  
  ! Time step size
  real(DP), intent(out) :: dTimeStep
  
  ! Re-initialization parameters
  integer, intent(out) :: RI_MAX,reini_freq,reiType
  
  ! BCs
  integer, intent(out) :: detBCs
  
  ! Collection structure for callback routines  
  type(t_collection), intent(out) :: rcollection
 !</output>

 !</subroutine>

  ! Local variables
  
  ! Time-dependent parameters
  real(DP) :: dTheta
  integer :: detTimeScale, detFEM
  
  ! Scaling factor for re-initialization
  real(DP) :: mbb
  
  ! reading data from the *.dat file 
  call parlst_init(rparams)
  call parlst_readfromfile (rparams, "./data/levelset.dat")
  
  ! Time-dependent parameters
  call parlst_getvalue_double (rparams, 'TIME', 'dTheta', dTheta, 0.5_DP)
  call parlst_getvalue_double (rparams, 'TIME', 'dTimeStep', dTimeStep, 0.1_DP)
  call parlst_getvalue_int (rparams, 'TIME', 'NMaxTime', NMaxTime, 10)
  
  call parlst_getvalue_int (rparams, 'TIME', 'reiType', reiType, 0)
  call parlst_getvalue_int (rparams, 'TIME', 'reini_freq', reini_freq, 10)
  call parlst_getvalue_int (rparams, 'TIME', 'RI_MAX', RI_MAX, 6)
  call parlst_getvalue_double (rparams, 'TIME', 'mbb', mbb, 1.0_DP)  
  
  ! Standard or Galerkin FEM
  call parlst_getvalue_int (rparams, 'MESH', 'detFEM', detFEM, 1)
  
  ! Assigne the theta scheme parameter
  rcollection%DquickAccess(1) = dTheta
  
  ! Also, pass the time step to the collection
  rcollection%DquickAccess(2) = dTimeStep

  ! Passing the current time
  rcollection%DquickAccess(3) = 0.0_DP
  
  ! The re-initialization parameter scale
  rcollection%DquickAccess(4) = mbb
  
  ! Re-initialization parameter, whether it is active or not
  ! in the current time step
  rcollection%IquickAccess(1) = 0
  
  ! Either using Standard or Least-squares
  ! FEM formulation.
  rcollection%IquickAccess(2) = detFEM

  ! Do we have BCs.
  call parlst_getvalue_int (rparams, 'MESH', 'detBCs', detBCs, 0)
  
  end subroutine
  

  !****************************************************************************
  
!<subroutine>
  subroutine ls_grid(rtriangulation,rboundary,rparams,rcollection,NLMAX) !,rindicator)

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
    
  ! Max. level to solve the problem
  integer, intent(out) :: NLMAX
  
  ! The || grad(\phi) - 1 ||_{L_2,T} vector
  !type(t_vectorScalar), intent(out) :: rindicator
  
 !</output>

 !<inputoutput>
  ! Collection structure for callback routines    
  type(t_collection), intent(inout) :: rcollection
 !</inputoutput>
 

 !</subroutine>

  ! Local variables
  ! Path to the mesh  
  character(LEN=SYS_STRLEN) :: sfile
   
  ! Mesh perturbance percentage
  real(DP) :: dPert
    
  ! We want to solve our NS problem on level...
  call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
  
  ! At first, read in the parametrisation of the boundary and save
  ! it to rboundary.
  call parlst_getvalue_string (rparams, 'MESH', &
            'sFilenamePathMesh',sfile, """""", bdequote=.true.)    
  call boundary_read_prm(rboundary, trim(sfile)//'.prm')

  ! Now read in the basic triangulation
  call tria_readTriFile2D (rtriangulation,&
                trim(sfile)//'.tri', rboundary)
  
  ! Refine the mesh up to the coarse grid level
  call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
  
  ! Create information about adjacencies and everything one needs from
  ! a triangulation. Afterwards, we have the coarse mesh.
  call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

  ! Create a vector based on the triangulation to
  ! store the || grad(\phi) - 1 ||_{L_2,T} on each element
  ! call lsyssc_createVector(rindicator,rtriangulation%NEL,.true.)


  ! Perturb the mesh on the coarse level so that all other meshes
  !  in the hierarchy have the same structure
  call parlst_getvalue_double (rparams, 'MESH', 'dPert', dPert, 0.0_DP)
  if (dPert .gt. 0.0_DP) then
    call meshmod_disturbMesh(rtriangulation,dPert)
  end if
   
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_structure(rtriangulation,rboundary,rparams,rdiscretisation,&
                                     rmatrix, rmatrix_projection,rcubatureInfo)

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
  ! An object specifying the discretisation (structure of the
  ! solution, trial/test functions,...)
  type(t_blockDiscretisation), intent(out) :: rdiscretisation
  
  ! A system matrix.
  type(t_matrixBlock), intent(out) :: rmatrix
  
  ! A system matrix.
  type(t_matrixBlock), intent(out) :: rmatrix_projection
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(out) :: rcubatureInfo
 !</output>

 
!</subroutine>
 
  ! Local variables
  ! String variable
  character(len=SYS_STRLEN) :: sstring
  
  ! Type of finite elements
  integer(I32) :: Pelm

  ! Type of cubature rule for numerical integration
  integer(I32) :: ccubType
  
  ! Jump stabilization parameters
  integer :: detPJump
  
  ! Max level
  integer :: NLMAX
  
  ! We want to solve our NS problem on level...
  call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
  
  ! Now we can start to initialise the discretisation. At first, set up
  ! a block discretisation structure that specifies 4 blocks in the
  ! solution vector.
  call spdiscr_initBlockDiscr (rdiscretisation,1,&
                 rtriangulation, rboundary)
  
  ! Read the finite element for velocities
  call parlst_getvalue_string (rparams, 'MESH', 'Pelm', sstring)
  Pelm = elem_igetID(sstring)


  ! Here we set up one discretisation structure for the 
  ! 1st component of vector variable
  call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
    Pelm,rtriangulation, rboundary)
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Set up an cubature info structure to tell the code which cubature
  ! formula to use on all levels
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Create an assembly information structure which tells the code
  ! the cubature formula to use.
  ! Get the integration rule from the data file
  call parlst_getvalue_string (rparams,'MESH','ccubType', sstring)
  ccubType = cub_igetID(sstring)  
   call spdiscr_createDefCubStructure(rdiscretisation%RspatialDiscr(1),&
   rcubatureInfo,ccubType)
  
  ! Initialise the block matrix with default values based on
  ! the discretisation.
  call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)
  
  call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix_projection)
  
  ! Now as the discretisation is set up, we can start to generate
  ! the structure of the system matrix which is to solve.
  !
  ! Let's check if we have to set up jump stabilization
  ! If so, we need to define the matrix structure accordingly
  call parlst_getvalue_int (rparams, 'JUMP', 'detPJump', detPJump, 0)  

  ! Velocity jump stabilization
  if (detPJump .eq. 1) then  
    call bilf_createMatrixStructure (&
    rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
    rmatrix%RmatrixBlock(1,1),cconstrType=BILF_MATC_EDGEBASED)
  else 
    call bilf_createMatrixStructure (&
    rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
    rmatrix%RmatrixBlock(1,1))  
  end if
  ! Allocate memory and zero-valued
  call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(1,1), LSYSSC_SETM_ZERO)  
  
  ! The projection matrix
  call bilf_createMatrixStructure (&
    rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
    rmatrix_projection%RmatrixBlock(1,1))
  
  ! Allocate memory and zero-valued
  call lsyssc_allocEmptyMatrix (rmatrix_projection%RmatrixBlock(1,1),&
                          LSYSSC_SETM_ZERO)
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Init_RhsndSolution(rdiscretisation,rmatrix_projection,rrhs, &
            rvector_old,rvector_oldT,rvector,rvectorU,rvectorV, rparams,&
            rcubatureInfo,rdiscreteBC,rboundary)
                
 !<description>  
  ! Initializing the RHS and the solution vector
 !</description>                

 !<output>
  ! Block vectors
  type(t_vectorBlock), intent(out) :: rvector_old,rvector_oldT,rvector
  type(t_vectorBlock), intent(out) :: rvectorU,rvectorV,rrhs
  
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC),intent(out), target :: rdiscreteBC  
 !</output>
 
 !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams 
  
  ! An object specifying the discretisation (structure of the
  ! solution, trial/test functions,...)
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
  
  ! A projection matrix.
  type(t_matrixBlock), intent(inout) :: rmatrix_projection
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
  
  ! An object for saving the domain:
  type(t_boundary), intent(in) :: rboundary  
 !</input>

!</subroutine>

  ! local variables
  type(t_boundaryRegion) :: rboundaryRegion  
  
  ! Let's start
  ! Create a RHS and solution vectors based on the discretisation.
  ! Fill with zero.
  call lsysbl_createVectorBlock (rdiscretisation, rrhs,.true.)
  call lsysbl_createVectorBlock (rdiscretisation, rvector,.true.)
  call lsysbl_createVectorBlock (rdiscretisation, rvector_old,.true.)
  call lsysbl_createVectorBlock (rdiscretisation, rvector_oldT,.true.)
  call lsysbl_createVectorBlock (rdiscretisation, rvectorU,.true.)
  call lsysbl_createVectorBlock (rdiscretisation, rvectorV,.true.)
    
  ! Bulid the mass matrix
  call bma_buildMatrix (rmatrix_projection,BMA_CALC_STANDARD,ls_Mass,&
      rcubatureInfo=rcubatureInfo)
  
  ! do the projection
  call anprj_analytL2projectionByMass (rvector_old%RvectorBlock(1),& 
   rmatrix_projection%RmatrixBlock(1,1), analyt_project)

  call anprj_analytL2projectionByMass (rvectorU%RvectorBlock(1),& 
   rmatrix_projection%RmatrixBlock(1,1), analyt_project_U)
  
  call anprj_analytL2projectionByMass (rvectorV%RvectorBlock(1),& 
   rmatrix_projection%RmatrixBlock(1,1), analyt_project_V)

  ! ====================
  ! Boundary conditions
  ! ====================
  call bcasm_initDiscreteBC(rdiscreteBC)

!  call boundary_createRegion(rboundary,1,1,rboundaryRegion)
!  rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
!  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                     rboundaryRegion,rdiscreteBC,&
!                                     getBoundaryValues_LS)
!                           
!  ! edge 2 of boundary component 1.
!  call boundary_createregion(rboundary,1,4,rboundaryRegion)
!  rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
!  call bcasm_newdirichletbconrealbd (rdiscretisation,1,&
!                                     rboundaryRegion,rdiscreteBC,&
!                                     getBoundaryValues_LS)
                             

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_MatAssembly(rmatrix,rvector_old,rvectorU,rvectorV,&
                            rcollection,rparams,NLMAX,rcubatureInfo)
                    
 !<description>  
  ! Initializing the solution vectors on all levels by calculating the 
  ! memory required to the interlevel projections.
 !</description>                
 
 !<input>
  ! Level info.
  integer, intent(in) :: NLMAX
  
  ! RHS block vector
  type(t_vectorBlock), intent(inout), target :: rvector_old,rvectorU,rvectorV
  
  ! Collection structure for callback routines  
  type(t_collection) :: rcollection
  
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams  
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo  
 !</input>
  
 !<inputoutput>
  ! System matrix.
  type(t_matrixBlock), intent(inout) :: rmatrix
 !</inputoutput>

!</subroutine>

  ! Local variables
  ! Collection of vectors to evaluate in nonlinear terms
  type(t_fev2Vectors) :: revalVectors
  
  ! Preparing a collection of vectors to evaluate nonlinear terms
  call ls_vec_collection(revalVectors,rvector_old,rvectorU,rvectorV)

  if (rcollection%IquickAccess(2) == 1) then
    ! LSFEM
    ! Assemble the whole system matrix on each level  
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,ls_levelset_mat,&
       rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
       revalVectors=revalVectors)
  else
    ! Galerkin
    ! Assemble the whole system matrix on each level  
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,std_levelset_mat,&
       rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
       revalVectors=revalVectors)
  end if     

  
  ! Set up jump stabilization if there is
  call ls_jump(rmatrix,rparams,rcollection)
   
  ! Release the vector structure used in linearization
  call fev2_releaseVectorList(revalVectors)
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_RHS_Assembly(rrhs,rvector_old,rvector_oldT,rvectorU,rvectorV,&
  rcollection,rcubatureInfo)
                                
 !<description>  
    ! Initializing the solution vectors on all levels by calculating the 
    ! memory required to the interlevel projections.
 !</description>                                
 
 !<input>
    ! Cubature information structure which defines the cubature formula.
    type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
 !</input> 
 
 !<inputoutput>
    ! RHS/solution block vectors
    type(t_vectorBlock), intent(inout) :: rrhs, rvector_old,rvector_oldT
    type(t_vectorBlock), intent(inout) :: rvectorU, rvectorV
    
    ! Collection structure for callback routines    
    type(t_collection) :: rcollection  
 !</inputoutput>
    
!</subroutine>


    ! Local variables    
    ! Collection of vectors to evaluate in nonlinear terms
    type(t_fev2Vectors) :: revalVectors
    
    ! Preparing a collection of vectors to evaluate nonlinear terms
    if (rcollection%IquickAccess(1) == 0) then
      ! No re-initialization considered :-(
      call ls_vec_collection(revalVectors,rvector_old,rvectorU,rvectorV)
    else
      ! Re-initialization considered :-)
      call ls_vec_collection(revalVectors,rvector_old,rvectorU,rvectorV,rvector_oldT)
    end if

  if (rcollection%IquickAccess(2) == 1) then
    ! LSFEM
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,ls_levelset_rhs,&
         rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
         revalVectors=revalVectors)    
  else
    ! Galerkin
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,std_levelset_rhs,&
         rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
         revalVectors=revalVectors)
  end if             
    
    ! Release the vector structure used in linearization
    call fev2_releaseVectorList(revalVectors)
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_Dirichlet(rmatrix,rrhs,rvector,rvector_old,&
                                              rvector_oldT,rdiscreteBC)
                                
 !<description>  
 ! Implementing BCs to the matrix and solution/RHS vectors.
 !</description>                                
 
 !<input> 
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC),intent(in) :: rdiscreteBC 
 !</input>
 
 !<inputoutput>
  ! Block vectors 
  type(t_vectorBlock), intent(inout) :: rvector,rvector_old,rrhs,rvector_oldT
  
  ! Block matrix
  type(t_matrixBlock), intent(inout) :: rmatrix
 !</inputoutput>

 
!</subroutine>
    
    ! Local variables
    
    ! Assign the boundary conditions to the matrix
    call lsysbl_assignDiscreteBC(rmatrix,rdiscreteBC)
    ! Implement the filter
    call matfil_discreteBC (rmatrix)
           
    ! Assign the boundary conditions to the vectors.
    call lsysbl_assignDiscreteBC(rrhs,rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector,rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector_old,rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector_oldT,rdiscreteBC)
    
    ! Implement the filter  
    call vecfil_discreteBCrhs (rrhs)
    call vecfil_discreteBCsol (rvector)
    call vecfil_discreteBCsol (rvector_old)
    call vecfil_discreteBCsol (rvector_oldT)

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Solver(rmatrix,rrhs,rvector,rparams)
                
 !<description>
  ! Set up a linear solver, solve the problem, release the solver.
 !</description>

 !<inputoutput>
  ! Solution Vector  
  type(t_vectorBlock), intent(inout) :: rvector
 !</inputoutput>
 
  !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams
  
   ! RHS vector
  type(t_vectorBlock), intent(inout) :: rrhs  
  
  ! System matrix.
  type(t_matrixBlock), intent(in) :: rmatrix
  !</input>
 
!</subroutine>

  ! Local variables
  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  
  ! A temporary vector
  type(t_vectorBlock) :: rtempBlock
  
  ! Error indicator during initialisation of the solver
  integer :: ierror  

  type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner

  ! Level info.
  integer :: NLMAX
  
  ! Read in some parameters
  ! Level info.
  call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
  call linsol_initUMFPACK4 (p_rsolverNode)
  
  call linsol_setMatrix(p_RsolverNode,rmatrix)
  
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
  subroutine ls_vec_collection(revalVectors,rvector_old,rvectorU,rvectorV,&
                                rvector_oldT)
                
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
 
  !<inputoutput>
   ! Solution vector in the current nonliner/time iteration  
  type(t_vectorBlock), intent(inout) :: rvector_old,rvectorU,rvectorV
  
  type(t_vectorBlock), intent(inout), optional :: rvector_oldT
  !</inputoutput>
 
!</subroutine>
 
  ! The routine <verb>fev2_addVectorToEvalList</verb> allows to define
  ! the evaluation of functions and 1st/2nd derivatives.
  
  ! First we add the nonlinear vectors
  call fev2_addVectorToEvalList(revalVectors,&
     rvectorU%RvectorBlock(1),0)       ! u    (1)
  call fev2_addVectorToEvalList(revalVectors,&
     rvectorV%RvectorBlock(1),0)       ! v    (2)
  call fev2_addVectorToEvalList(revalVectors,&
     rvector_old%RvectorBlock(1),1)    ! phi  (3)   
     
  ! Then, we add the time discretization vector if it is passed
  if (present(rvector_oldT)) then
    call fev2_addVectorToEvalList(revalVectors,&
      rvector_oldT%RvectorBlock(1),1)  ! phi  (4)
  end if
     
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
    integer :: detPJump
    real(DP) :: dJumpP, dJumpStarP, deojEdgeExpP
    
    ! Let's check if we realy have to set up jump stabilization   
    call parlst_getvalue_int (rparams, 'JUMP', 'detPJump', detPJump, 0) 
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Velocity jump stabilization
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    if (detPJump .eq. 1) then
        call parlst_getvalue_double (rparams, 'JUMP', 'dJumpP', &
                                                    dJumpP, 0.01_DP)  
        call parlst_getvalue_double (rparams, 'JUMP', 'dJumpStarP',&
                                                dJumpStarP, 0.0_DP)  
        call parlst_getvalue_double (rparams, 'JUMP', 'deojEdgeExpP',&
                                                deojEdgeExpP, 2.0_DP)                                                        

        ! Set up the jump stabilisation structure.
        ! The kinematic viscosity 1/Re
        rjumpStabil%dnu = 1.0_DP

        ! Set stabilisation parameter
        rjumpStabil%dgamma = dJumpP
        rjumpStabil%dgammastar = dJumpStarP
        rjumpStabil%deojEdgeExp = deojEdgeExpP

        ! Matrix weight, =0 no jump stabilization will be added
        rjumpStabil%dtheta = 1.0_DP

        ! Cubature formula to be used in jump term calculations
        ! over the edges
        rjumpStabil%ccubType = CUB_G3_1D

        ! Call the jump stabilisation technique for the 1st velocity.
        call conv_jumpStabilisation2d (rjumpStabil, CONV_MODMATRIX, &
                                            rmatrix%RmatrixBlock(1,1))
    end if
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_update_solution(rmatrix,rvector,rvector_old,rrhs,rcollection,&
                                                                i,RI_Loop_MAX)
                
 !<description>  
  ! Update Solution
 !</description>                

 !<inputoutput>  
   ! Solution vectors in the current/previous nonliner iterations
   ! and the RHS vector
  type(t_vectorBlock), intent(inout) :: rvector,rvector_old,rrhs
  
  ! Block matrix
  type(t_matrixBlock), intent(inout) :: rmatrix
 !<\inputoutput>
  

 !<input>  
  ! Re-initialization loop control parameters
  integer, intent(in) :: i, RI_Loop_MAX
  
  ! Collection structure for callback routines  
  type(t_collection) :: rcollection  
  
 !<\input>

!</subroutine>
    
  ! Local variables
  real(DP) :: Dres(1),Dresv(1),Dres_rel(1)
  integer, dimension(1) :: Cnorms
  
  ! Scaling factors
  real(DP) :: cx,cy
  
  ! Difference of the vectors in the current nonlinear iteration
  type(t_vectorBlock) :: rdiff  
  
  ! Print 
  if (rcollection%IquickAccess(1) == 1) then   
    
    ! Euclidian vector norm: (vector,vector) 0
    ! $l_2$-norm: 1/sqrt(NEQ) * (vector,vector) 2
    ! max-norm 3   
    ! Normal L^2 Norm
    Cnorms(1) = 0
    
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
    
    ! Some output data
    if (i .eq. 1) then
        call output_line ('Iter. ' //' Phi Rel. Err. ')
        call output_line ('--------------------------')
        call output_line (sys_siL(i, 5) //'  '&
        //trim(sys_sdEL(Dres_rel(1),6)))   
    else
        call output_line (sys_siL(i, 5) //'  '&
        //trim(sys_sdEL(Dres_rel(1),6)))   
        if ( (mod(i,10) .eq. 0) .and. (i .ne. RI_Loop_MAX) ) then
            call output_lbrk()
            call output_line ('Iter. ' &
            //' U1 Rel. Err. ')
            call output_line ('--------------------------')
        end if
    end if
    
    ! Release the block vector
    call lsysbl_releaseVector (rdiff) 
    
 end if    
    
    
  
  ! Copy the current solution to the old solution
  call lsysbl_copyVector (rvector,rvector_old)
  !*** Clear all data in matrix and RHS ***!
  call lsysbl_clearVector (rrhs)
  call lsysbl_clearMatrix (rmatrix)
 
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_postprocess(rboundary,rmatrix,rvector,rtriangulation,&
                                  rcubatureInfo,rdiscretisation,rparams,itime)
                                
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
    
    ! An optional parameter, if it is set the routine is used 
    !  to write the time-dependent data (GMV/VTK)
    integer, intent(in), optional :: itime
    
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
    integer :: ExporType
    integer :: Ptild
           
     ! Path to the data file which has the initial solution
    character(LEN=SYS_STRLEN) :: sfile
 
    ! Output block for UCD output to GMV/VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! An object specifying the discretisation.
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! A block vector which contains projected data
    type(t_vectorBlock) :: rprjVector

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    
    character(len=SYS_STRLEN) :: sstring
    
    if (present(itime)) then
      
      ! ############################################
      ! Use this routine to ONLY write the gmv data
      !  in a time-dependent solution
      ! ############################################
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Writing the solution to GMV/VTK files.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! We can now start the exporting the results.
      ! Get the path for writing postprocessing files from the environment variable
      ! $UCDDIR. If that does not exist, write to the directory "./gmv".
      if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

      ! Determine which type of export do we use, GMV/VTK
      call parlst_getvalue_int (rparams, 'POST', 'ExporType', ExporType, 0)

      ! Detemine whether we need to project the solution to a GMV/VTK readable
      ! type.
      call parlst_getvalue_int (rparams, 'MESH', 'Ptild', Ptild, 0)

      if ( Ptild .eq. 1) then
              
          ! make a new discretization structure for the projected data
          ! and modify its sub-structures if required
          call spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
          
          if (Ptild .eq. 1) then
              call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
              EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(1))                 
          endif        
       
          ! Now set up a new solution vector based on this discretisation,
          ! allocate memory.
          call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)

          ! Then take our original solution vector and convert it according to the
          ! new discretisation:
          call spdp_projectSolution (rvector,rprjVector)
          
          select case (ExporType)
            case (0)
          
              ! Start UCD export to VTK file:
              call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.vtk')

              ! Write Pressure
              call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'phi',UCD_VAR_STANDARD,p_Ddata)
          
            case (1)
              
              ! Start UCD export to GMV file:
              call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

              ! Write Pressure
              call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'phi',UCD_VAR_STANDARD,p_Ddata)
              
      
            case (2)
              
              ! Start UCD export to GMV file, Binary GMV
              call ucd_startBGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

              ! Write Pressure
              call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'phi',UCD_VAR_STANDARD,p_Ddata)
              
            case default
              call output_line ('Invalid visualisation output type.', &
                                OU_CLASS_ERROR,OU_MODE_STD,'ls_postprocess')
              call sys_halt()
          end select
                    
          ! Release the temporary projected vector
          call lsysbl_releaseVector (rprjVector)
          
          ! Release the projected discretisation structure and 
          ! all spatial discretisation structures in it.
          call spdiscr_releaseBlockDiscr(rprjDiscretisation)
      
      else ! real data will be used
      
          select case (ExporType)
            case (0)
          
              ! Start UCD export to VTK file:
              call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.vtk')

              ! Write Pressure
              call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'phi',UCD_VAR_STANDARD,p_Ddata)
          
          case (1)
              ! Start UCD export to GMV file:
              call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

              ! Write Pressure
              call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'phi',UCD_VAR_STANDARD,p_Ddata)
                      
          
          case (2)
              ! Start UCD export to GMV file, Binary GMV
              call ucd_startBGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

              ! Write Pressure
              call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
                          
            case default
              call output_line ('Invalid visualisation output type.', &
                                OU_CLASS_ERROR,OU_MODE_STD,'ls_postprocess')
              call sys_halt()
          
          end select  
      
      end if ! end of real or projected data condition
        
      ! Write the file to disc, that is it.
      call ucd_write (rexport)
      call ucd_release (rexport)
    
      
    else
      ! ##########################################
      ! This is a full post-processing subroutine 
      ! ##########################################
    end if
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_cleanup(rvector,rvector_old,rvector_oldT,rvectorU,rvectorV,&
  rrhs,rmatrix,rmatrix_projection,rboundary,rparams,rdiscretisation,&
  rtriangulation,rcubatureInfo,rdiscreteBC)
 
 !<description>  
    ! Release all the memory used in our calculations.
 !</description> 

  !<input>
    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams
  !</input>  

  !<inputoutput>
  ! An object for saving the domain:
  type(t_boundary) :: rboundary
    
  ! An object for saving the triangulation on the domain
  type(t_triangulation), intent(inout) :: rtriangulation    
    
  ! A couple of block vectors.
  type(t_vectorBlock) :: rvector,rvector_old,rvector_oldT,rvectorU,rvectorV,rrhs
    
  ! An object specifying the discretisation (structure of the
  ! solution, trial/test functions,...)
  type(t_blockDiscretisation), intent(inout) :: rdiscretisation
  
  ! A system matrix.
  type(t_matrixBlock), intent(inout) :: rmatrix, rmatrix_projection
  
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(inout) :: rcubatureInfo    
  
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC),intent(inout), target :: rdiscreteBC
  !</inputoutput>

!</subroutine>
    
    ! Now, clean up so that all the memory is available again.    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rvector_old)
    call lsysbl_releaseVector (rvector_oldT)
    call lsysbl_releaseVector (rvectorU)
    call lsysbl_releaseVector (rvectorV)
    call lsysbl_releaseVector (rrhs)

    call lsysbl_releaseMatrix (rmatrix)
    call lsysbl_releaseMatrix (rmatrix_projection)

    call bcasm_releaseDiscreteBC (rdiscreteBC)


    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo)

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
  subroutine ls_levelset_mat(RmatrixData,rassemblyData,rmatrixAssembly,&
                         npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
  ! This is where the BIG-BANG starts :D
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
  real(DP) :: dbasI, dbasJ, dval, dbasIx, dbasIy, dbasJx, dbasJy
  integer :: iel, icubp, idofe, jdofe
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11 
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2,P_Dphi
  
  ! Real time data
  real(DP) :: Dt, dPI
     
  ! Time-dependent parameters
  real(DP) :: Dtstp, Dtheta
    
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dPx, dPy
    
  ! The re-initialization scale and,
  ! Re-initialization parameter, whether it is active or not
  ! in the current time step
  real(DP) :: mbb, dNx, dNy
  integer :: det
  
  mbb = rcollection%DquickAccess(4)
  det = rcollection%IquickAccess(1)
  
  ! The value of \theta in our temporal discretization
  Dtheta = rcollection%DquickAccess(1)
  
  ! The time step size = \nabla t
  Dtstp = rcollection%DquickAccess(2)
  
  ! Real simulation time
  Dt = rcollection%DquickAccess(3)
  
  ! PI number
  dPI = 3.1415926535897932_DP
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
  
  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
    
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
    
  ! Get the velocity field from the parameters 
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! We do the calculation in a block-by-block manner. All the
  ! relevant blocks will be calculated in the same loop over the
  ! elements.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  ! ++++++++++++++++++++
  ! Calculate blocks A11
  ! ++++++++++++++++++++
  ! Loop over the elements in the current set.
  
  if (det == 0) then
    ! +++++++++++++++++++++++++++++++
    ! NO RE-INITIALIZATION CONSIDERED 
    ! +++++++++++++++++++++++++++++++
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
            
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
        dval = p_DcubWeight(icubp,iel) * (   dbasJ*dbasI  + &
        Dtheta*Dtstp*(dbasJ*dU*dbasIx + dbasJ*dV*dbasIy)  + &
        Dtheta*Dtstp*(dbasJx*dU*dbasI + dbasJy*dV*dbasI)  + &
         (Dtheta*Dtstp)*(Dtheta*Dtstp)*( dU*dU*dbasJx*dbasIx+ &
        dU*dV*dbasJx*dbasIy + dV*dU*dbasJy*dbasIx + dV*dV*dbasJy*dbasIy )  )
             
        p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
                        
        end do ! idofe
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
  
  else

    ! +++++++++++++++++++++++++++++++
    ! RE-INITIALIZATION IS CONSIDERED 
    ! +++++++++++++++++++++++++++++++
    ! Previous re-initialization loop
    P_Dphi=> revalVectors%p_RvectorData(3)%p_Ddata 
    
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
      
      dPx = P_Dphi(icubp,iel,DER_DERIV2D_X)
      dPy = P_Dphi(icubp,iel,DER_DERIV2D_Y)
      
      dNx = dPx/(DSQRT(dPx**2 + dPy**2))
      dNy = dPy/(DSQRT(dPx**2 + dPy**2))
      
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
        dval = p_DcubWeight(icubp,iel) * (   dbasJ*dbasI  + &
        Dtheta*Dtstp*(dbasJ*dU*dbasIx + dbasJ*dV*dbasIy)  + &
        Dtheta*Dtstp*(dbasJx*dU*dbasI + dbasJy*dV*dbasI)  + &
         (Dtheta*Dtstp)*(Dtheta*Dtstp)*( dU*dU*dbasJx*dbasIx+ &
        dU*dV*dbasJx*dbasIy + dV*dU*dbasJy*dbasIx + dV*dV*dbasJy*dbasIy )  )

        dval = dval + mbb*p_DcubWeight(icubp,iel) * (   dNx**2*dbasJx*dbasIx +&
        dNx*dNy*dbasJx*dbasIy + dNy*dNx*dbasJy*dbasIx + dNy**2*dbasJy*dbasIy   )
        
        p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
                        
        end do ! idofe
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
  end if  
  
  end subroutine

  !****************************************************************************


!<subroutine>
  subroutine ls_levelset_rhs(rvectorData,rassemblyData,rvectorAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the RHS vector in a block-by-block procedures.
  ! The rest of the BIG-BANG happens to occure here :D
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
  real(DP) :: dbasI,dbasIx,dbasIy, dval
  integer :: iel, icubp, idofe
  real(DP), dimension(:,:), pointer :: p_DlocalVector
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaVectorData), pointer :: p_rvectorData

  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2,P_Dphi,P_DphiN
   
  ! Time-dependent parameters
  real(DP) :: Dtstp, Dtheta
    
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dPx, dPy, dphi, dPxN, dPyN
    
  ! Real time data
  real(DP) :: Dt, dPI    
    
  ! The re-initialization scale and,
  ! Re-initialization parameter, whether it is active or not
  ! in the current time step
  real(DP) :: mbb, dNx, dNy
  integer :: det
  
  mbb = rcollection%DquickAccess(4)
  det = rcollection%IquickAccess(1)    
    
  ! The value of \theta in our temporal discretization
  Dtheta = rcollection%DquickAccess(1)
  
  ! The time step size = \nabla t
  Dtstp = rcollection%DquickAccess(2)
  
  ! Real simulation time
  Dt = rcollection%DquickAccess(3)
  
  ! PI number
  dPI = 3.1415926535897932_DP
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rvectorData => RvectorData(1)  
  p_DlocalVector => RvectorData(1)%p_Dentry
  p_DbasTest => RvectorData(1)%p_DbasTest
  
  ! Get the velocity field from the parameters 
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
  ! Calculate the RHS of the velocities
  
  if (det == 0) then
    ! +++++++++++++++++++++++++++++++
    ! NO RE-INITIALIZATION CONSIDERED 
    ! +++++++++++++++++++++++++++++++  
    ! The previous time step known values of \phi
    P_Dphi=> revalVectors%p_RvectorData(3)%p_Ddata
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
      
      dphi  = P_Dphi(icubp,iel,DER_FUNC)
      dPx = P_Dphi(icubp,iel,DER_DERIV2D_X)
      dPy = P_Dphi(icubp,iel,DER_DERIV2D_Y)
          
      ! Outer loop over the DOF's i=1..ndof on our current element,
      ! which corresponds to the (test) basis functions Phi_i:
      do idofe=1,p_rvectorData%ndofTest
      
        ! Fetch the contributions of the (test) basis functions Phi_i
        ! into dbasI
        dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
        dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
        dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
                  
        ! Values of the velocity RHS for the X1 and X2 component
        dval = dbasI*dphi - (1.0_DP-Dtheta)*Dtstp*  &
        (dbasI*dU*dPx + dbasI*dV*dPy) + Dtheta*Dtstp*dphi*    &
        (dU*dbasIx + dV*dbasIy) - Dtheta*(1.0_DP-Dtheta)*Dtstp*Dtstp*  &
        ( dU*dU*dbasIx*dPx + dU*dV*dbasIx*dPy + dV*dU*dbasIy*dPx + &
         dV*dV*dbasIy*dPy)
            
        ! Multiply the values of the basis functions by
        ! the cubature weight and sum up into the local vectors.
        p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
          p_DcubWeight(icubp,iel) * dval
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
  
  else
    ! +++++++++++++++++++++++++++++++
    ! RE-INITIALIZATION IS CONSIDERED 
    ! +++++++++++++++++++++++++++++++
    ! The previous re-initialization loops' known values of \phi
    P_DphiN=> revalVectors%p_RvectorData(3)%p_Ddata
    
    ! The previous time step known values of \phi
    P_Dphi=> revalVectors%p_RvectorData(4)%p_Ddata
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
      
      ! Previous time step 
      dphi  = P_Dphi(icubp,iel,DER_FUNC)
      dPx = P_Dphi(icubp,iel,DER_DERIV2D_X)
      dPy = P_Dphi(icubp,iel,DER_DERIV2D_Y)

      ! Previous nonlinear loop
      dPxN = P_DphiN(icubp,iel,DER_DERIV2D_X)
      dPyN = P_DphiN(icubp,iel,DER_DERIV2D_Y)
      
      dNx = dPxN/(DSQRT(dPxN**2 + dPyN**2))
      dNy = dPyN/(DSQRT(dPxN**2 + dPyN**2))      
          
      ! Outer loop over the DOF's i=1..ndof on our current element,
      ! which corresponds to the (test) basis functions Phi_i:
      do idofe=1,p_rvectorData%ndofTest
      
        ! Fetch the contributions of the (test) basis functions Phi_i
        ! into dbasI
        dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
        dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
        dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
                  
        ! Values of the velocity RHS for the X1 and X2 component
        dval = dbasI*dphi - (1.0_DP-Dtheta)*Dtstp*  &
        (dbasI*dU*dPx + dbasI*dV*dPy) + Dtheta*Dtstp*dphi*    &
        (dU*dbasIx + dV*dbasIy) - Dtheta*(1.0_DP-Dtheta)*Dtstp*Dtstp*  &
        ( dU*dU*dbasIx*dPx + dU*dV*dbasIx*dPy + dV*dU*dbasIy*dPx + &
         dV*dV*dbasIy*dPy)

        dval = dval + mbb*(dNx*dbasIx + dNy*dbasIy)
        
        ! Multiply the values of the basis functions by
        ! the cubature weight and sum up into the local vectors.
        p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
          p_DcubWeight(icubp,iel) * dval
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
  
  end if  
  
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

 
  ! ***************************************************************************

!<subroutine>
  subroutine std_levelset_mat(RmatrixData,rassemblyData,rmatrixAssembly,&
                         npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
  ! This is where the BIG-BANG starts :D
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
  real(DP) :: dbasI, dbasJ, dval, dbasIx, dbasIy, dbasJx, dbasJy
  integer :: iel, icubp, idofe, jdofe
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11 
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2,P_Dphi
  
  ! Real time data
  real(DP) :: Dt, dPI
     
  ! Time-dependent parameters
  real(DP) :: Dtstp, Dtheta
    
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dPx, dPy
    
  ! The re-initialization scale and,
  ! Re-initialization parameter, whether it is active or not
  ! in the current time step
  real(DP) :: mbb, dNx, dNy
  integer :: det
  
  mbb = rcollection%DquickAccess(4)
  det = rcollection%IquickAccess(1)
  
  ! The value of \theta in our temporal discretization
  Dtheta = rcollection%DquickAccess(1)
  
  ! The time step size = \nabla t
  Dtstp = rcollection%DquickAccess(2)
  
  ! Real simulation time
  Dt = rcollection%DquickAccess(3)
  
  ! PI number
  dPI = 3.1415926535897932_DP
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
  
  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
    
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
    
  ! Get the velocity field from the parameters 
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! We do the calculation in a block-by-block manner. All the
  ! relevant blocks will be calculated in the same loop over the
  ! elements.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  ! ++++++++++++++++++++
  ! Calculate blocks A11
  ! ++++++++++++++++++++
  ! Loop over the elements in the current set.
  
  if (det == 0) then
    ! +++++++++++++++++++++++++++++++
    ! NO RE-INITIALIZATION CONSIDERED 
    ! +++++++++++++++++++++++++++++++
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
            
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
        dval = p_DcubWeight(icubp,iel) * (   dbasJ*dbasI  + &
        Dtheta*Dtstp*(dbasJx*dU*dbasI + dbasJy*dV*dbasI)  )    
        p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
                        
        end do ! idofe
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
  
  else

    ! +++++++++++++++++++++++++++++++
    ! RE-INITIALIZATION IS CONSIDERED 
    ! +++++++++++++++++++++++++++++++
    ! Previous re-initialization loop
    P_Dphi=> revalVectors%p_RvectorData(3)%p_Ddata 
    
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
      
      dPx = P_Dphi(icubp,iel,DER_DERIV2D_X)
      dPy = P_Dphi(icubp,iel,DER_DERIV2D_Y)
      
      dNx = dPx/(DSQRT(dPx**2 + dPy**2))
      dNy = dPy/(DSQRT(dPx**2 + dPy**2))
      
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
        dval = p_DcubWeight(icubp,iel) * (   dbasJ*dbasI  + &
        Dtheta*Dtstp*(dbasJx*dU*dbasI + dbasJy*dV*dbasI)  )    
        p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval

        dval = dval + mbb*p_DcubWeight(icubp,iel) * (   dNx**2*dbasJx*dbasIx +&
        dNx*dNy*dbasJx*dbasIy + dNy*dNx*dbasJy*dbasIx + dNy**2*dbasJy*dbasIy   )
        
        p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
                        
        end do ! idofe
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
  end if  
  
  end subroutine

  !****************************************************************************


!<subroutine>
  subroutine std_levelset_rhs(rvectorData,rassemblyData,rvectorAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the RHS vector in a block-by-block procedures.
  ! The rest of the BIG-BANG happens to occure here :D
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
  real(DP) :: dbasI,dbasIx,dbasIy, dval
  integer :: iel, icubp, idofe
  real(DP), dimension(:,:), pointer :: p_DlocalVector
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaVectorData), pointer :: p_rvectorData

  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2,P_Dphi,P_DphiN
   
  ! Time-dependent parameters
  real(DP) :: Dtstp, Dtheta
    
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dPx, dPy, dphi, dPxN, dPyN
    
  ! Real time data
  real(DP) :: Dt, dPI    
    
  ! The re-initialization scale and,
  ! Re-initialization parameter, whether it is active or not
  ! in the current time step
  real(DP) :: mbb, dNx, dNy
  integer :: det
  
  mbb = rcollection%DquickAccess(4)
  det = rcollection%IquickAccess(1)    
    
  ! The value of \theta in our temporal discretization
  Dtheta = rcollection%DquickAccess(1)
  
  ! The time step size = \nabla t
  Dtstp = rcollection%DquickAccess(2)
  
  ! Real simulation time
  Dt = rcollection%DquickAccess(3)
  
  ! PI number
  dPI = 3.1415926535897932_DP
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rvectorData => RvectorData(1)  
  p_DlocalVector => RvectorData(1)%p_Dentry
  p_DbasTest => RvectorData(1)%p_DbasTest
  
  ! Get the velocity field from the parameters 
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
  ! Calculate the RHS of the velocities
  
  if (det == 0) then
    ! +++++++++++++++++++++++++++++++
    ! NO RE-INITIALIZATION CONSIDERED 
    ! +++++++++++++++++++++++++++++++  
    ! The previous time step known values of \phi
    P_Dphi=> revalVectors%p_RvectorData(3)%p_Ddata
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
      
      dphi  = P_Dphi(icubp,iel,DER_FUNC)
      dPx = P_Dphi(icubp,iel,DER_DERIV2D_X)
      dPy = P_Dphi(icubp,iel,DER_DERIV2D_Y)
          
      ! Outer loop over the DOF's i=1..ndof on our current element,
      ! which corresponds to the (test) basis functions Phi_i:
      do idofe=1,p_rvectorData%ndofTest
      
        ! Fetch the contributions of the (test) basis functions Phi_i
        ! into dbasI
        dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                  
        ! Values of the velocity RHS for the X1 and X2 component
        dval = dbasI*dphi - (1.0_DP-Dtheta)*Dtstp*  &
        (dbasI*dU*dPx + dbasI*dV*dPy)
            
        ! Multiply the values of the basis functions by
        ! the cubature weight and sum up into the local vectors.
        p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
          p_DcubWeight(icubp,iel) * dval
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
  
  else
    ! +++++++++++++++++++++++++++++++
    ! RE-INITIALIZATION IS CONSIDERED 
    ! +++++++++++++++++++++++++++++++
    ! The previous re-initialization loops' known values of \phi
    P_DphiN=> revalVectors%p_RvectorData(3)%p_Ddata
    
    ! The previous time step known values of \phi
    P_Dphi=> revalVectors%p_RvectorData(4)%p_Ddata
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

      ! Velocity/derivatives field in this cubature point
      dU = p_Du1(icubp,iel,DER_FUNC)
!      dU = dU * DCOS(dPI * Dt/5.0_DP)
      
      dV = p_Du2(icubp,iel,DER_FUNC)
!      dV = dV * DCOS(dPI * Dt/5.0_DP)
      
      ! Previous time step 
      dphi  = P_Dphi(icubp,iel,DER_FUNC)
      dPx = P_Dphi(icubp,iel,DER_DERIV2D_X)
      dPy = P_Dphi(icubp,iel,DER_DERIV2D_Y)

      ! Previous nonlinear loop
      dPxN = P_DphiN(icubp,iel,DER_DERIV2D_X)
      dPyN = P_DphiN(icubp,iel,DER_DERIV2D_Y)
      
      dNx = dPxN/(DSQRT(dPxN**2 + dPyN**2))
      dNy = dPyN/(DSQRT(dPxN**2 + dPyN**2))      
          
      ! Outer loop over the DOF's i=1..ndof on our current element,
      ! which corresponds to the (test) basis functions Phi_i:
      do idofe=1,p_rvectorData%ndofTest
      
        ! Fetch the contributions of the (test) basis functions Phi_i
        ! into dbasI
        dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
        dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
        dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
                  
        ! Values of the velocity RHS for the X1 and X2 component
        dval = dbasI*dphi - (1.0_DP-Dtheta)*Dtstp*  &
        (dbasI*dU*dPx + dbasI*dV*dPy)

        dval = dval + mbb*(dNx*dbasIx + dNy*dbasIy)
        
        ! Multiply the values of the basis functions by
        ! the cubature weight and sum up into the local vectors.
        p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
          p_DcubWeight(icubp,iel) * dval
        
      end do ! jdofe

      end do ! icubp
    
    end do ! iel
  
  end if  
  
  end subroutine

end module

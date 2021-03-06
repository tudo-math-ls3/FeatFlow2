!##############################################################################
!# ****************************************************************************
!# <name> LS_NS_VVP_Time_MG2D </name>
!# ****************************************************************************
!# <purpose>
!# This module solves the 2D Navier-Stokes (NS) eq. using LSFEM.
!#
!# The second-order elliptic NS equations are reformulated into 
!# first-order equations using the definition of the vorticity:
!#   vorticity:  w = curl(u).
!#   see the following link for the definition of vorticity:
!#   http://www.student.math.uwaterloo.ca/
!#                 ~amat361/Fluid%20Mechanics/topics/vorticity.htm
!# The resulting system is called Velocity-Vorticity-Pressure (VVP).
!# The problem is solved in a coupled manner for the solution of:
!#   1- velocity components   u1, u2 (in 2D)
!#   2- vorticity function    w
!#   3- pressure              p
!# variables.
!#
!# The nonlinear term is first linearized using Newton/Fixed-point method.
!# The LSFEM formulation then applied which yiedls a symmetric-
!# positive definite linear system of equations.
!# This routine uses the Multigrid as linear solver.
!# The discretisation uses the block assembly method to evaluate the
!# mattrix all-in-one.
!# </purpose>
!#
!# Author:        Masoud Nickaeen
!# First Version: July 18, 2012
!# Last Update:   Jan. 27, 2013
!##############################################################################

module LS_NS_VVP_Time_MG2D

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
  
  implicit none

!<types>

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! An object specifying the discretisation of the Lumped Mass matrix
    type(t_blockDiscretisation) :: rdiscMass
    
    ! A Lumped Mass matrix for that specific level.
    type(t_matrixBlock) :: rmatrixMass

    ! A variable describing the discrete boundary conditions.    
    type(t_discreteBC) :: rdiscreteBC

    ! Temporary vector in the size of the RHS/solution vector on that level.
    type(t_vectorBlock) :: rtempVector
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), pointer :: p_rprojection => null()

    ! Cubature information structure which defines the cubature formula.
    type(t_scalarCubatureInfo) :: rcubatureInfo
  
  end type
  
!</types>

contains
  
  !****************************************************************************

!<subroutine>
  subroutine ls_vvp_time_mg2d
  
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
  ! An array of problem levels for the multigrid solver
  type(t_level), dimension(:), pointer :: Rlevels    
        
  ! An object for saving the domain:
  type(t_boundary) :: rboundary
    
  ! A couple of block vectors.
  type(t_vectorBlock) :: rvector,rvector_old,rvector_oldT,rrhs

  ! Temporary scalar vector; used for calculating the nonlinear matrix
  ! on lower levels / projecting the solution from higher to lower levels.
  type(t_vectorScalar), pointer :: p_rtempVectorSc => null()
    
  ! Max. Min. level to solve the problem
  integer :: NLMAX, NLMIN
    
  ! Loop index
  integer :: i

  ! Parameters used in time loop
  integer :: itime,NMaxTime
  real(DP) :: dTimeStep
    
  ! Parameters used in nonliner loop
  integer :: inl,NLN_Max
  real(DP) :: dNLEpsi
  
  ! Convergence parameter, either by error or by NLN_Max
  logical :: converged,diverged

  ! Nonlinear loop control
  logical :: det
    
  ! Collection structure for callback routines
  type(t_collection) :: rcollection
  
  ! All parameter of the LSFEM solver
  type(t_parlist) :: rparams
  
  ! Ok, let's start.
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 0)-Read all the parameters from data file and initialize the collection.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call ls_initParams(rparams,NMaxTime,dTimeStep,NLN_MAX,dNLEpsi,rcollection)
  

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! a)-Read the domain, the mesh and refine it. All about computational GRID.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_grid(Rlevels,rboundary,rparams,rcollection,NLMAX,NLMIN)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! b)- Set up a discretisation and cubature info structure which tells 
  ! the code which finite element to use.
  ! Also, create a 4*4 block matrix structure.
  ! Initialize the structure of the deferred velocities to be use in the
  !  nonlinear matrix assembly routine. This is done for all grid levels.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_structure(Rlevels,rboundary,rparams)

  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! c)- Initialization of Boundary conditions, and the solution/RHS vectors.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize the discrete boundary conditions
  do i = NLMIN, NLMAX
    call ls_BCs_Dirichlet_one(Rlevels(i)%rdiscretisation,rboundary,&
                                                  Rlevels(i)%rDiscreteBC)
  end do
  
  ! Initialize the RHS and fine level solution vectors
  call ls_Init_RhsndSolution(Rlevels,rrhs,&
                         rvector_old, rvector_oldT,rvector,rparams,NLMAX)
  
  ! Initialize the memory required for interlevel projection of the
  !  solution vectors on all levels except the finest one.
  call ls_Init_CoarseSolution(Rlevels,rrhs,p_rtempVectorSc,NLMAX,NLMIN)

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The time loop, the outermost loop, starts here.
  ! The loop performs maximum 'NMaxTime' times.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do itime=1,NMaxTime
     
    call output_separator(OU_SEP_MINUS)
    call output_line ('Time step '//trim(sys_siL(itime,6))// &
                      '     Time '//trim(sys_sdL(itime*dTimeStep,5)))
    call output_lbrk ()
    
    ! Copy the initial/previous solution to the current solution
    call lsysbl_copyVector (rvector_old,rvector_oldT)
    
    det = .false.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! d)- Here the non-linear loop starts.
    ! The loop performs a maximum of NLN_Max iterations. The stopping criterion
    ! is controled with dNLEpsi.
    ! 
    ! Every loop performs the following series of operations:
    !   1- System matrix assembly (requires the evaluation of 
    !                              nonlinear deferred velocities)
    !   1-1 And, jump stabilization set up if required
    !   2- RHS assembly (nonlinear deferred velocities must be released 
    !                    right after this step!!)
    !   3- Boundary conditions implementation, excluding the one time
    !       calculation of descretised Dirichlet boundary conditions which is
    !       done earlier in 'ls_BCs_Dirichlet_One' subroutine
    !   4- Solver setup, solution of the final system, solver release
    !   5- Check for the non-linear loop convergence/divergence
    !   6- Update initial guess 'if (.not. converged) .and. (.not. diverged)'
    !      (all matrix and RHS vectors must be cleaned, 
    !       zero-valued, after this!!)
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    do inl=1,NLN_Max
    
      ! +++++++++++++++++++++++++
      ! 1- System Matrix Assembly
      ! +++++++++++++++++++++++++
      call ls_MatAssembly(Rlevels,rvector_old,p_rtempVectorSc,&
                                  rcollection,rparams,NLMIN,NLMAX,inl)
      
      
      ! +++++++++++++++
      ! 2- RHS Assembly
      ! +++++++++++++++
      call ls_RHS_Assembly(rrhs,rvector_old,rvector_oldT,rcollection,&
                                     Rlevels(NLMAX)%rcubatureInfo,inl) 
          
          
      ! ++++++++++++++++++++++++++++++++
      ! 3- Implement Boundary Conditions
      ! ++++++++++++++++++++++++++++++++     
      ! 3-1 Dirichlet BCs.
      call ls_BCs_Dirichlet(Rlevels,rrhs,rvector,rvector_old,&
                                            rvector_oldT,NLMAX,NLMIN)
      
      ! 3-2 Neuman, Outflow and Etc. BCs.
      call ls_BCs_NOE(Rlevels,rrhs,rboundary,rparams,NLMAX,NLMIN)
      
      ! 3-3 Zero Mean Pressure constraint
      call ls_BCs_ZMP(Rlevels,rrhs,rboundary,inl,rparams,NLMAX,NLMIN)     
       
      ! +++++++++++++++++++
      ! 4- Solve the System
      ! +++++++++++++++++++
      call ls_Solver_linear(Rlevels,rrhs,rvector,rparams)
      
      
      ! ++++++++++++++++++++++++++++++++++++++
      ! 5- Check for Convergence or Divergence
      ! ++++++++++++++++++++++++++++++++++++++
      call ls_con_di_verge(converged,diverged,rvector,&
                                   rvector_old,NLN_Max,inl,dNLEpsi)
      
      
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 6- Update the Initial Solution, clear matrix, vectors
      !  or Terminate the Loop
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
      call ls_update_solution(Rlevels,converged,diverged,&
                       rvector,rvector_old,rrhs,NLMAX,NLMIN,det)
      if (det) exit
          
    end do  ! inl
    
    ! Write the time-dependent data
    call ls_postprocess(rboundary,Rlevels(NLMAX)%rmatrix,&
       rvector,Rlevels(NLMAX)%rtriangulation,&
       Rlevels(NLMAX)%rcubatureInfo,&
       Rlevels(NLMAX)%rdiscretisation,rparams, itime)    
    
  end do  ! itime

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! e)- Post-processing starts here. Writing the solution to file ...
  ! export GMV/VTK files and calculate drag/lift forces ...
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_postprocess(rboundary,Rlevels(NLMAX)%rmatrix,&
       rvector,Rlevels(NLMAX)%rtriangulation,&
       Rlevels(NLMAX)%rcubatureInfo,Rlevels(NLMAX)%rdiscretisation,rparams)
  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! f)- Clean up, free all the memory which is used for variables.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_cleanup(rvector,rvector_old,rvector_oldT,rrhs,Rlevels,&
                             p_rtempVectorSc,rboundary,NLMAX,NLMIN,rparams)

  end subroutine


  !****************************************************************************
  
!<subroutine>
  subroutine ls_initParams(rparams,NMaxTime,dTimeStep,NLN_MAX,dNLEpsi,&
                                                                rcollection)

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

  ! Time loop Max. number of iterations
  integer, intent(out) :: NMaxTime
  
  ! Time step size
  real(DP), intent(out) :: dTimeStep
  
  ! Collection structure for callback routines  
  type(t_collection), intent(out) :: rcollection
 !</output>

 !</subroutine>

  ! Local variables
  ! Kinematic Viscosity
  real(DP) :: dnu
  
  ! Continuity equation scaling factor
  !  default = 1.0_DP  no scaling
  real(DP) :: alpha
  
  ! Other scaling factors
  integer :: scPhysic, scADN  
  
  ! Linearization method parameters
  integer :: LinTyp, FPIter
  
  ! Time-dependent parameters
  real(DP) :: dTheta
  integer :: detTimeScale
  
  
  ! reading data from the *.dat file 
  call parlst_init(rparams)
  call parlst_readfromfile (rparams, "./data/lsvvp.dat") 
     
  ! Viscosity parameter: noo = 1/Re
  call parlst_getvalue_double (rparams, 'GFPROPER', 'dnu', dnu, 1.0_DP)   
  
  ! Nonlinear loop Max. number of iterations
  call parlst_getvalue_int (rparams, 'NLLOOP', 'NLN_MAX', NLN_MAX, 3)
  
  ! Nonlinear loop stopping criterion
  call parlst_getvalue_double (rparams, 'NLLOOP', 'dNLEpsi', dNLEpsi, 1.0_DP)
  
  ! Continuity equation scaling factor
  call parlst_getvalue_double (rparams, 'GFPROPER', 'alpha', alpha, 1.0_DP)
  
  ! Linearization method
  call parlst_getvalue_int (rparams, 'NLLOOP', 'LinTyp', LinTyp, 0)
  call parlst_getvalue_int (rparams, 'NLLOOP', 'FPIter', FPIter, 0)
  
  ! Time-dependent parameters
  call parlst_getvalue_int (rparams, 'TIME', 'detTimeScale', detTimeScale, 0)
  call parlst_getvalue_double (rparams, 'TIME', 'dTheta', dTheta, 0.5_DP)
  call parlst_getvalue_double (rparams, 'TIME', 'dTimeStep', dTimeStep, 0.1_DP)
  call parlst_getvalue_int (rparams, 'TIME', 'NMaxTime', NMaxTime, 10)
  
  ! Other scalings, scale least-squares functionals or not.
  call parlst_getvalue_int (rparams, 'GFPROPER', 'scPhysic', scPhysic, 0)
  call parlst_getvalue_int (rparams, 'GFPROPER', 'scADN', scADN, 0)
   
  
  ! Initializing the collection
  ! Put kinematic viscosity there, to be used in nonlinear assembly
  rcollection%DquickAccess(1) = dnu
 
  ! Physical scaling
  if (scPhysic .eq. 1) then
      rcollection%DquickAccess(2) = 1.0_DP/(dnu*dnu)
  else
      rcollection%DquickAccess(2) = 1.0_DP
  end if
  
  ! ADN theory scaling
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
  
  ! Assigne the theta scheme parameter
  rcollection%DquickAccess(5) = dTheta
  
  ! Also, pass the time step to the collection
  rcollection%DquickAccess(6) = dTimeStep
  
  ! Check whether we have to scale based on the time scales
  if (detTimeScale .eq. 1) then
    rcollection%DquickAccess(7) = dTimeStep/dnu
  else
    rcollection%DquickAccess(7) = 1.0_DP
  end if  
  
  end subroutine
  

  !****************************************************************************
  
!<subroutine>
  subroutine ls_grid(Rlevels,rboundary,rparams,rcollection,NLMAX,NLMIN)

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

    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
    
    ! Max. level to solve the problem
    integer, intent(out) :: NLMAX,NLMIN
 !</output>

 !<inputoutput>
    ! Collection structure for callback routines    
    type(t_collection), intent(inout) :: rcollection
 !</inputoutput>


 !</subroutine>

    ! Local variables
    ! Path to the mesh    
    character(LEN=SYS_STRLEN) :: sfile
    
    ! Loop index
    integer :: i
    
    ! Mesh statistics determin
    integer :: detPMesh, detTri
    
    ! Mesh perturbance percentage
    real(DP) :: dPert
        
    ! We want to solve our NS problem on level...
    call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
    call parlst_getvalue_int (rparams, 'MESH', 'NLMIN', NLMIN, 3)
    call parlst_getvalue_int (rparams, 'MESH', 'detTri', detTri, 0)
        
    ! Allocate memory for all levels
    allocate(Rlevels(NLMIN:NLMAX))   
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call parlst_getvalue_string (rparams, 'MESH', &
                      'sFilenamePathMesh',sfile, """""", bdequote=.true.)      
    call boundary_read_prm(rboundary, trim(sfile)//'.prm')

    ! Now read in the basic triangulation, put it into NLMIN triangulation
    call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation,&
                                trim(sfile)//'.tri', rboundary)
    
    ! Refine the mesh up to the coarse grid level
    call tria_quickRefine2LevelOrdering (NLMIN-1,&
                       Rlevels(NLMIN)%rtriangulation,rboundary)

    ! Changing Quads to Tris if requested
    ! Each quadrilateral is divided into two triangles
    if (detTri == 1)  call tria_rawGridToTri (Rlevels(NLMIN)%rtriangulation)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
                                                             rboundary)

    ! Perturb the mesh on the coarse level so that all other meshes
    !  in the hierarchy have the same structure
    call parlst_getvalue_double (rparams, 'MESH', 'dPert', dPert, 0.0_DP)
    if (dPert .gt. 0.0_DP) then
        call meshmod_disturbMesh(Rlevels(NLMIN)%rtriangulation,dPert)
    end if
    
    ! Now refine the grid for the fine levels.
    do i = NLMIN+1, NLMAX

      ! Refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
                                  Rlevels(i)%rtriangulation,rboundary)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
                                                            rboundary)

    end do
   
    ! Compress the level hierarchy.
    ! Share the vertex coordinates of all levels, so the coarse grid coordinates
    ! are 'contained' in the fine grid coordinates. The effect is:
    ! 1.) Save some memory
    ! 2.) Every change in the fine grid coordinates also affects the coarse
    !     grid coordinates and vice versa.
    do i=NLMAX-1,NLMIN,-1
      call tria_compress2LevelOrdHierarchy (Rlevels(i+1)%rtriangulation,&
          Rlevels(i)%rtriangulation)
    end do   
    
    ! Print mesh statistics or not
    call parlst_getvalue_int (rparams, 'MESH', 'detPMesh', detPMesh, 0)
    
    ! Print mesh information
    if (detPMesh .eq. 1) then
        call output_lbrk ()
        call output_line ('Mesh statistics:')
        call output_lbrk ()
        do i=NLMIN,NLMAX
        call tria_infoStatistics (Rlevels(i)%rtriangulation,&
            i .eq. NLMIN,i)
        end do
   end if
   
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_structure(Rlevels,rboundary,rparams)

 !<description>  
    ! Set up a discretisation and cubature info structure. This tells
    ! the code which finite element to use.
    ! Also, create a block matrix structure.
 !</description>                                

 !<input>
    ! An object for saving the domain:
    type(t_boundary), intent(in) :: rboundary

    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams 
 !</input>

 !<inputoutput>
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels  
    
 !</inputoutput>

 
!</subroutine>
 
    ! Local variables
    ! String variable
    character(len=SYS_STRLEN) :: sstring
        
    ! Type of finite elements
    integer(I32) :: Velm, Pelm, Welm    

    ! Type of cubature rule for numerical integration
    integer(I32) :: ccubType
    
    ! Jump stabilization parameters
    integer :: detVJump, detWJump, detPJump
    
    ! Loop index
    integer :: i
    
    ! Min, Max levels
    integer :: NLMIN, NLMAX
    
    ! Print discretization statistics or not
    integer :: detPDisc
    
    ! We want to solve our NS problem on level...
    call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
    call parlst_getvalue_int (rparams, 'MESH', 'NLMIN', NLMIN, 3)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 4 blocks in the
    ! solution vector.
    ! Do this for all levels
    do i = NLMIN, NLMAX
      call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation,4,&
                                 Rlevels(i)%rtriangulation, rboundary)
    end do


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


    ! Here we set up one discretisation structure for the 
    ! 1st component of vector variable
    do i = NLMIN, NLMAX
      call spdiscr_initDiscr_simple (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                Velm, Rlevels(i)%rtriangulation, rboundary)
    end do
                    
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd component (Y-velocity). This needs no additional memory,
    ! as both structures will share the same dynamic information afterwards.
    do i = NLMIN, NLMAX
        call spdiscr_duplicateDiscrSc(&
               Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                     Rlevels(i)%rdiscretisation%RspatialDiscr(2))
    end do
    
    ! For the pressure (3rd component), we set up a separate discretisation
    ! structure, as this 'MAY' use different finite elements for trial and test
    ! functions.
    ! Read the finite element for Pressure
    call parlst_getvalue_string (rparams, 'MESH', 'Pelm', sstring)
    Pelm = elem_igetID(sstring)

    do i = NLMIN, NLMAX
        call spdiscr_deriveSimpleDiscrSc (&
                Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                    Pelm, Rlevels(i)%rdiscretisation%RspatialDiscr(3))
    end do
        
    ! And for Vorticity, a separate doscretisation structure as well.
    ! Read the finite element for Vorticity
    call parlst_getvalue_string (rparams, 'MESH', 'Welm', sstring)
    Welm = elem_igetID(sstring)  

    do i = NLMIN, NLMAX    
        call spdiscr_deriveSimpleDiscrSc (&
                Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                    Welm, Rlevels(i)%rdiscretisation%RspatialDiscr(4))
    end do
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use on all levels
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use.
    ! Get the integration rule from the data file
    call parlst_getvalue_string (rparams,'MESH','ccubType', sstring)
    ccubType = cub_igetID(sstring)        
    do i = NLMIN, NLMAX    
       call spdiscr_createDefCubStructure(&  
       Rlevels(i)%rdiscretisation%RspatialDiscr(1),Rlevels(i)%rcubatureInfo,&
                                                                   ccubType)                                                                   
    end do    
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a block matrix structure for the system matrix on all levels
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    do i = NLMIN, NLMAX
        ! Initialise the block matrix with default values based on
        ! the discretisation.
        call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscretisation,&
                                                      Rlevels(i)%rmatrix)
        
        ! Now as the discretisation is set up, we can start to generate
        ! the structure of the system matrix which is to solve.
        !
        ! The global system looks like this, a full 4*4 block matrix
        ! which is dymmetric and positive definite.
        !
        !    ( A11 A12 A13 A14 )
        !    ( A21 A22 A23 A24 )
        !    ( A31 A32 A33 A34 )
        !    ( A41 A42 A43 A44 )
        !    
        ! Create the matrix structure of the X-velocity. Block A11,
        !
        ! Let's check if we have to set up jump stabilization
        ! If so, we need to define the matrix structure accordingly
        call parlst_getvalue_int (rparams, 'JUMP', 'detVJump', detVJump, 0)    

        ! Velocity jump stabilization
        if (detVJump .eq. 1) then  
            call bilf_createMatrixStructure (&
            Rlevels(i)%rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
            Rlevels(i)%rmatrix%RmatrixBlock(1,1),cconstrType=BILF_MATC_EDGEBASED)
        else 
            call bilf_createMatrixStructure (&
            Rlevels(i)%rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
            Rlevels(i)%rmatrix%RmatrixBlock(1,1))    
        end if
                                       
        ! Block A12,
        call bilf_createMatrixStructure (&
              Rlevels(i)%rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
                 Rlevels(i)%rmatrix%RmatrixBlock(1,2))    
        ! Block A13
        call bilf_createMatrixStructure (&
              Rlevels(i)%rdiscretisation%RspatialDiscr(3), LSYSSC_MATRIX9, &
              Rlevels(i)%rmatrix%RmatrixBlock(1,3),&
              Rlevels(i)%rdiscretisation%RspatialDiscr(1))       
        ! Block A14
        call bilf_createMatrixStructure (&
              Rlevels(i)%rdiscretisation%RspatialDiscr(4), LSYSSC_MATRIX9, &
              Rlevels(i)%rmatrix%RmatrixBlock(1,4),&
              Rlevels(i)%rdiscretisation%RspatialDiscr(1))  
                                                                        
                                                                        
        ! Use X-velocity structure for the Y-velocity. Block A22,
        call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
            Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        ! Block A23,
        call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,3),&
            Rlevels(i)%rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        ! Block A24,
        call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,4),&
            Rlevels(i)%rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)


        ! Create the matrix structure of the Pressure. Block A33,
        ! Let's check if we have to set up jump stabilization
        ! If so, we need to define the matrix structure accordingly
        call parlst_getvalue_int (rparams, 'JUMP', 'detPJump', detPJump, 0)    

        ! Velocity jump stabilization
        if (detPJump .eq. 1) then
          call bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(3), LSYSSC_MATRIX9, &
          Rlevels(i)%rmatrix%RmatrixBlock(3,3),cconstrType=BILF_MATC_EDGEBASED)        
        else
          call bilf_createMatrixStructure (&
                 Rlevels(i)%rdiscretisation%RspatialDiscr(3), LSYSSC_MATRIX9, &
                    Rlevels(i)%rmatrix%RmatrixBlock(3,3))
        end if
        
        ! Block A34,
        call bilf_createMatrixStructure (&
          Rlevels(i)%rdiscretisation%RspatialDiscr(4), LSYSSC_MATRIX9, &
             Rlevels(i)%rmatrix%RmatrixBlock(3,4), &
                Rlevels(i)%rdiscretisation%RspatialDiscr(3))


        ! Create the matrix structure of the Vorticity. Block A44,
        ! 
        ! Let's check if we have to set up jump stabilization
        ! If so, we need to define the matrix structure accordingly
        call parlst_getvalue_int (rparams, 'JUMP', 'detWJump', detWJump, 0)    

        ! Velocity jump stabilization
        if (detWJump .eq. 1) then     
            call bilf_createMatrixStructure (&
            Rlevels(i)%rdiscretisation%RspatialDiscr(4), LSYSSC_MATRIX9, &
            Rlevels(i)%rmatrix%RmatrixBlock(4,4),cconstrType=BILF_MATC_EDGEBASED)
        else 
            call bilf_createMatrixStructure (&
            Rlevels(i)%rdiscretisation%RspatialDiscr(4), LSYSSC_MATRIX9, &
            Rlevels(i)%rmatrix%RmatrixBlock(4,4))    
        end if
            
            
        ! Create the structure for the transpoed matrices by transposing the structre
        ! of the relevant matrices. Block A21,
        call lsyssc_transposeMatrix (&
            Rlevels(i)%rmatrix%RmatrixBlock(1,2), &
            Rlevels(i)%rmatrix%RmatrixBlock(2,1),LSYSSC_TR_STRUCTURE)
        
        
        ! Block A31,
        call lsyssc_transposeMatrix (&
            Rlevels(i)%rmatrix%RmatrixBlock(1,3), &
            Rlevels(i)%rmatrix%RmatrixBlock(3,1),LSYSSC_TR_STRUCTURE)
        ! Block A32,
        call lsyssc_duplicateMatrix (&
        Rlevels(i)%rmatrix%RmatrixBlock(3,1),&
        Rlevels(i)%rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
            
            
        ! Block A41,
        call lsyssc_transposeMatrix (&
            Rlevels(i)%rmatrix%RmatrixBlock(1,4), &
            Rlevels(i)%rmatrix%RmatrixBlock(4,1),LSYSSC_TR_STRUCTURE)
        ! Block A42,
        call lsyssc_duplicateMatrix (&
        Rlevels(i)%rmatrix%RmatrixBlock(4,1),&
        Rlevels(i)%rmatrix%RmatrixBlock(4,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        ! Block A43    
        call lsyssc_transposeMatrix (&
            Rlevels(i)%rmatrix%RmatrixBlock(3,4), &
            Rlevels(i)%rmatrix%RmatrixBlock(4,3),LSYSSC_TR_STRUCTURE)
        
        ! Now re-assign the block discretisation structure to all matrices
        call lsysbl_assignDiscrDirectMat (Rlevels(i)%rmatrix, &
                                          Rlevels(i)%rdiscretisation)
        
        ! Allocate memory for the matrix
        call lsysbl_allocEmptyMatrix (Rlevels(i)%rmatrix,LSYSSC_SETM_ZERO)
        
    end do

    ! Print discretization statistics
    call parlst_getvalue_int (rparams, 'MESH', 'detPDisc', detPDisc, 0)
    if (detPDisc .eq. 1) then
      call output_lbrk ()
      call output_line ('Discretisation statistics:')
      do i=NLMIN,NLMAX
        call output_lbrk ()
        call output_line ('Level '//sys_siL(i,5))
        call dof_infoDiscrBlock (Rlevels(i)%rdiscretisation,.false.)
        call output_lbrk ()
      end do
    end if


    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create the deferred velocity vector structure on all levels except
    !  the finest level.
    ! Initialize the projection structure on all levels except
    !  the coarset one.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Initialise the deferred velocities structures.
    do i=NLMIN,NLMAX
      
      if (i .lt. NLMAX) then
        call lsysbl_createVecBlockByDiscr (Rlevels(i)%rdiscretisation,&
                                          Rlevels(i)%rtempVector,.true.)
      end if
      
    end do
    
    ! Initialise the standard interlevel projection structures.
    do i = NLMIN+1, NLMAX

        ! Allocate the projection structure
        allocate(Rlevels(i)%p_rprojection)

        ! Initialise the projection based on the discretisation.
        call mlprj_initProjectionDiscr (Rlevels(i)%p_rprojection,&
                                        Rlevels(i)%rdiscretisation)

        ! Initialise the projection structure with data from the INI/DAT
        ! files. This allows to configure prolongation/restriction.
        call ls_getProlRest (Rlevels(i)%p_rprojection,&
                                          rparams ,'PROLREST')

    end do
    
  end subroutine


! ***************************************************************************

!<subroutine>

  subroutine ls_getProlRest (rprojection, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>  
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(in) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  type(t_interlevelProjectionBlock), intent(inout) :: rprojection
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: i1
    real(DP) :: d1
    
    call parlst_querysection(rparamList, sname, p_rsection)
        
    ! Prolongation/restriction order for velocity components
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:2)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:2)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:2)%iinterpolationOrder = i1
    end if
    
    ! Prolongation/restriction variant for velocity components
    ! in case of Q1~ discretisation
    call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    if (i1 .ne. -1) then
      rprojection%RscalarProjection(:,1:2)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:2)%irestVariant  = i1
    end if
    
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction
    call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    if (i1 .ne. 1) then
      rprojection%RscalarProjection(:,1:2)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:2)%irestARIndicatorEX3Y  = i1
    end if

    ! Aspect-ratio bound for switching to constant prolongation/restriction
    ! in case of Q1~ discretisation with extended prolongation/restriction
    call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    if (d1 .ne. 20.0_DP) then
      rprojection%RscalarProjection(:,1:2)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:2)%drestARboundEX3Y  = d1
    end if

  end subroutine
    
    
  !****************************************************************************    

!<subroutine>
  subroutine ls_BCs_Dirichlet_One(rdiscretisation,rboundary,rdiscreteBC)
                                
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
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! edge 2 of boundary component 1.
    call boundary_createregion(rboundary,1,2,rboundaryregion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newdirichletbconrealbd (rdiscretisation,1,&
                                       rboundaryregion,rdiscretebc,&
                                       getboundaryvalues_2d)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

!    ! Edge 5 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,5,rboundaryRegion)
!    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND    
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!
!    ! Edge 6 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,6,rboundaryRegion)
!    rboundaryRegion%iproperties = 2**1-2**1
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 of boundary component 1.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)


!    ! Edge 5 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,5,rboundaryRegion)
!    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND    
!    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!
!    ! Edge 6 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,6,rboundaryRegion)
!    rboundaryRegion%iproperties = 2**1-2**1    
!    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)

!    ! Pressure BCs.
!    ! Edge 2 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!    rboundaryRegion%iproperties = 2**1-2**1
!    call bcasm_newDirichletBConRealBD (rdiscretisation,3,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)


!    ! Flow around cylinder
!    ! X-velocity
!    ! Edge 1 of boundary component 2. That is it.
!    call boundary_createRegion(rboundary,2,1,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!
!    ! Y-velocity
!    ! Edge 1 of boundary component 2. That is it.
!    call boundary_createRegion(rboundary,2,1,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)


    ! Proot 2006
    ! X-velocity
    ! Edge 1 of boundary component 2. That is it.
    call boundary_createRegion(rboundary,2,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    call boundary_createRegion(rboundary,2,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    call boundary_createRegion(rboundary,2,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    call boundary_createRegion(rboundary,2,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! Y-velocity
    ! Edge 1 of boundary component 2. That is it.
    call boundary_createRegion(rboundary,2,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    call boundary_createRegion(rboundary,2,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    call boundary_createRegion(rboundary,2,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    call boundary_createRegion(rboundary,2,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)                                       
                   
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Init_RhsndSolution(Rlevels,rrhs,rvector_old,rvector_oldT,&
                                               rvector,rparams,NLMAX)
                
 !<description>  
  ! Initializing the RHS and the solution vector
 !</description>                

 !<output>
  ! Block vectors
  type(t_vectorBlock), intent(out) :: rvector_old,rvector_oldT,rvector,rrhs
 !</output>
 
 !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams 

  ! Max. level to solve the problem
  integer, intent(in) :: NLMAX
 !</input>

 !<inputoutput>
  ! An array of problem levels for the multigrid solver
  type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>

!</subroutine>

  ! local variables
  integer :: nlinit, ilev, NEQ

  ! Initial value for the 1st step of nonliner loop
  real(DP) :: dinit_vect(4)

  ! Path to the data file which has the initial solution
  character(LEN=SYS_STRLEN) :: sfile, sstring, sarray
  
  ! Projection structure
  type(t_interlevelProjectionBlock) :: rprojection
  
  ! Temporary vectors
  type(t_vectorBlock), target :: rvector1,rvector2
  type(t_vectorScalar) :: rvectorTemp  
  
  ! Let's start
  ! Create a RHS and a solution vectors based on the discretisation.
  ! Fill with zero.
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                                                      rrhs,.true.)
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                                                   rvector,.true.)
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                                              rvector_oldT,.true.)
  
  ! Determine how to setup initial nonlinear solution
  call parlst_getvalue_int (rparams, 'ISOLUTION', 'nlinit', nlinit, 0)  
  
  if (nlinit .eq. 0) then
    
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                                                   rvector_old,.true.)    
    
    ! Initialize the solution vector(s) with constant values
    call parlst_getvalue_string (rparams, 'ISOLUTION', 'initValues',&
                               sstring, '0.0_DP 0.0_DP 0.0_DP 0.0_DP')
    read (sstring,*) dinit_vect(1), dinit_vect(2), dinit_vect(3), & 
                                                         dinit_vect(4)
    
    ! Scale the sub-vectors to initialize the nonlineaer iteration loop
    call lsyssc_clearVector (rvector_old%RvectorBlock(1),dinit_vect(1))
    call lsyssc_clearVector (rvector_old%RvectorBlock(2),dinit_vect(2))
    call lsyssc_clearVector (rvector_old%RvectorBlock(3),dinit_vect(3))
    call lsyssc_clearVector (rvector_old%RvectorBlock(4),dinit_vect(4))    

  else      
  
    ! Ignor the initial values, read from file
    !
    ! Read from data file the initial solution level
    call parlst_getvalue_int (rparams,'ISOLUTION',&
                                       'iinitialSolutionLevel',ilev,0)

    ! First creat a block vector structure and zero-valued it
    ! it must have the size of the initial solution
    call lsysbl_createVectorBlock (Rlevels(ilev)%rdiscretisation,&
                                                      rvector1,.true.) 
    call parlst_getvalue_string (rparams, 'ISOLUTION', &
                 'sFilenamePathResult',sfile, """""", bdequote=.true.)
    call vecio_readBlockVectorHR (rvector1, sarray, .true.,&
                                                     0, sfile, .true.)

    ! If the vector is on level < NLMAX,
    !  we have to bring it to level NLMAX
    do while (ilev .lt. NLMAX)

      ! Initialise a vector for the higher level and a 
      !  prolongation structure.
      call lsysbl_createVectorBlock (&
                  Rlevels(ilev+1)%rdiscretisation,rvector2,.false.)

      call mlprj_initProjectionVec (rprojection,rvector2)

      ! Prolongate to the next higher level.
      NEQ = mlprj_getTempMemoryVec (rprojection,rvector1,rvector2)
      if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
      call mlprj_performProlongation (rprojection,rvector1, &
                          rvector2,rvectorTemp)
      if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)

      ! Swap rvector1 and rvector2. Release the coarse grid vector.
      call lsysbl_swapVectors (rvector1,rvector2)
      call lsysbl_releaseVector (rvector2)

      call mlprj_doneProjection (rprojection)

      ! rvector1 is now on level ilev+1
      ilev = ilev+1

    end do

    ! Copy the resulting vector rvector1 to the output vector.
    call lsysbl_copyVector (rvector1,rvector_old)

    ! Release the temp vector
    call lsysbl_releaseVector (rvector1)
                      
  end if
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Init_CoarseSolution(Rlevels,rrhs,p_rtempVectorSc,NLMAX,NLMIN)
                                
 !<description>  
    ! Calculating the memory required for the interlevel projections.
 !</description>                                
 
 !<input>
    ! Level info.
    integer, intent(in) :: NLMAX,NLMIN 

    ! RHS block vector
    type(t_vectorBlock), intent(in) :: rrhs
 !</input>

 !<output>
    ! The temporary vector to be used in interlevel projections 
    type(t_vectorScalar), pointer :: p_rtempVectorSc
 !</output>
    
 !<inputoutput>
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>

!</subroutine>

    ! local variables
    integer :: imaxmem, i

    ! How much memory is necessary for performing the level change?
    ! We ourself must build nonlinear matrices on multiple levels and have
    ! to interpolate the solution vector from finer level to coarser ones.
    ! We need temporary memory for this purpose...

    imaxmem = 0
    do i=NLMIN+1,NLMAX
        ! Pass the system metrices on the coarse/fine grid to
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryDirect (&
            Rlevels(i)%p_rprojection,&
            Rlevels(i-1)%rdiscretisation%RspatialDiscr(1:rrhs%nblocks),&
            Rlevels(i)%rdiscretisation%RspatialDiscr(1:rrhs%nblocks)))
    end do

    ! Set up a scalar temporary vector that we need for building up nonlinear
    ! matrices.
    allocate(p_rtempVectorSc)
    call lsyssc_createVector (p_rtempVectorSc,imaxmem,.false.)
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_MatAssembly(Rlevels,rvector_old,p_rtempVectorSc,&
                                 rcollection,rparams,NLMIN,NLMAX,inl)
                                
 !<description>  
    ! Initializing the solution vectors on all levels by calculating the 
    ! memory required to the interlevel projections.
 !</description>                                
 
 !<input>
    ! Level info.
    integer, intent(in) :: NLMAX,NLMIN
    
    ! Current nonlinear loop iteration
    integer, intent(in) :: inl
    
    ! RHS block vector
    type(t_vectorBlock), intent(in), target :: rvector_old
    
    ! Collection structure for callback routines    
    type(t_collection) :: rcollection
    
    ! The temporary vector to be used in interlevel projections 
    type(t_vectorScalar), pointer :: p_rtempVectorSc
    
    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams    
 !</input>
    
 !<inputoutput>
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>

!</subroutine>

    ! Local variables
    ! Level indicator
    integer :: ilev
    
    ! Collection of vectors to evaluate in nonlinear terms
    type(t_fev2Vectors) :: revalVectors
    
    ! Temporary block matrix and vectors
    type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
    type(t_vectorScalar), pointer :: p_rvectorTemp
    type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse
    
    type(t_matrixBlock) :: rmatrix
    
    ! A pointer to the projection structure in each level
    type(t_interlevelProjectionBlock), pointer :: p_rprojection
    
    ! A filter chain for the linear solver
    type(t_filterChain), dimension(1), target :: RfilterChain


    ! Linearization Scheme
    ! If it is a mixed type, check when to shift to
    ! Newton's
    if (rcollection%IquickAccess(2) .eq. 3) then
      if (inl .gt. rcollection%IquickAccess(3)) then
        ! It's time to shift to Newton's method 
        rcollection%DquickAccess(4) = 1.0_DP
      end if
    end if
   
    ! Get the temporary vector from the collection.
    ! Our 'parent' prepared there how to interpolate the solution on the
    ! fine grid to coarser grids.
    p_rvectorTemp => p_rtempVectorSc

    ! On all levels, we have to set up the nonlinear system matrix,
    ! so that the linear solver can be applied to it.
    
    nullify(p_rmatrix)

    do ilev=NLMAX,NLMIN,-1
    
      ! Get the matrix on the current level.
      ! Shift the previous matrix to the pointer of the fine grid matrix.
      p_rmatrixFine => p_rmatrix
      p_rmatrix => Rlevels(ilev)%rmatrix
    
      ! On the highest level, we use rvector_old as solution to build the
      ! nonlinear matrix. On lower levels, we have to create a solution
      ! on that level from a fine-grid solution before we can use
      ! it to build the matrix!
      if (ilev .eq. NLMAX) then
      
        p_rvectorCoarse => rvector_old
        
      else
        ! We have to discretise a level hierarchy and are on a level < NLMAX.
        
        ! Get the projection structure for this level.
        p_rprojection => Rlevels(ilev+1)%p_rprojection

        ! Get the temporary vector on level i. Will receive the solution
        ! vector on that level.
        p_rvectorCoarse => Rlevels(ilev)%rtempVector
        
        ! Get the solution vector on level i+1. This is either the temporary
        ! vector on that level, or the solution vector on the maximum level.
        if (ilev .lt. NLMAX-1) then
          p_rvectorFine => Rlevels(ilev+1)%rtempVector
        else
          p_rvectorFine => rvector_old
        end if

        ! Interpolate the solution from the finer grid to the coarser grid.
        ! The interpolation is configured in the interlevel projection
        ! structure we got from the collection.
        call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                         p_rvectorFine,p_rvectorTemp)

      end if
    
      ! Preparing a collection of vectors to evaluate nonlinear terms
      call ls_vec_collection(revalVectors,p_rvectorCoarse)

      ! Assemble the whole system matrix on each level    
      call bma_buildMatrix (p_rmatrix,BMA_CALC_STANDARD,ls_ns2DT_Matrix,&
             rcubatureInfo=Rlevels(ilev)%rcubatureInfo,rcollection=rcollection, &
             revalVectors=revalVectors)

      ! Set up jump stabilization if there is
      call ls_jump(p_rmatrix,rparams,rcollection)
       
      ! Release the vector structure used in linearization
      call fev2_releaseVectorList(revalVectors)
        
    end do
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_RHS_Assembly(rrhs,rvector_old,rvector_oldT,rcollection,&
                                                          rcubatureInfo,inl)
                                
 !<description>  
    ! Initializing the solution vectors on all levels by calculating the 
    ! memory required to the interlevel projections.
 !</description>                                
 
 !<input>
    ! Cubature information structure which defines the cubature formula.
    type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo 
    
    ! Current nonlinear loop iteration
    integer, intent(in) :: inl     
 !</input> 
 
 !<inputoutput>
    ! RHS/solution block vectors
    type(t_vectorBlock), intent(inout) :: rrhs, rvector_old, rvector_oldT
    
    ! Collection structure for callback routines    
    type(t_collection) :: rcollection  
 !</inputoutput>
    
!</subroutine>


    ! Local variables    
    ! Collection of vectors to evaluate in nonlinear terms
    type(t_fev2Vectors) :: revalVectors
    
    ! Preparing a collection of vectors to evaluate nonlinear terms
    call ls_vec_collection(revalVectors,rvector_old,rvector_oldT)
      
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,ls_ns2DT_rhs,&
         rcubatureInfo=rcubatureInfo,rcollection=rcollection, &
         revalVectors=revalVectors)    
    
    ! Release the vector structure used in linearization
    call fev2_releaseVectorList(revalVectors)
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_Dirichlet(Rlevels,rrhs,rvector,rvector_old,rvector_oldT,&
                                                                  NLMAX,NLMIN)
                                
 !<description>  
    ! Implementing BCs to the matrix and solution/RHS vectors.
 !</description>                                

 !<input>
    integer, intent(in) :: NLMAX,NLMIN 
 !</input>
 
 !<inputoutput>
    ! Block vectors 
    type(t_vectorBlock), intent(inout) :: rvector,rvector_old,rvector_oldT,rrhs
     
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>

 
!</subroutine>
    
    ! Local variables
    ! Loop index
    integer :: i

    do i = NLMIN, NLMAX
    
        ! Assign the boundary conditions to the matrix on all levels.
        call lsysbl_assignDiscreteBC(Rlevels(i)%rmatrix,Rlevels(i)%rdiscreteBC)
        ! Implement the filter
        call matfil_discreteBC (Rlevels(i)%rmatrix)
        
    end do
    
    ! Assign the boundary conditions to the vectors ONLY on the finest level.
    call lsysbl_assignDiscreteBC(rrhs,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector_old,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector_oldT,Rlevels(NLMAX)%rdiscreteBC)
    
    ! Implement the filter  
    call vecfil_discreteBCrhs (rrhs)
    call vecfil_discreteBCsol (rvector)
    call vecfil_discreteBCsol (rvector_old)
    call vecfil_discreteBCsol (rvector_oldT)
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_ZMP(Rlevels,rrhs,rboundary,inl,rparams,NLMAX,NLMIN)
                                
 !<description>  
    ! Implementing the Zero Mean Pressure Constraint.
 !</description>                                

 !<input>
    ! An object for saving the domain:
    type(t_boundary), intent(in) :: rboundary
        
    ! Current nonlinear loop iteration
    integer, intent(in) :: inl
    
    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams

    ! Min, Max levels
    integer, intent(in) :: NLMIN, NLMAX
 !</input>
 
 !<inputoutput>
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
     
    ! RHS block vector
    type(t_vectorBlock), intent(inout) :: rrhs
 !</inputoutput>

 
!</subroutine>
    
    ! Local variables
    ! Loop index and Zero Mean pressure parameters
    integer :: i, detZMV, irow
    integer, dimension(1) :: Irows
    
    ! Pressure element type
    integer(I32) :: Pelm
    
    ! String variable
    character(len=SYS_STRLEN) :: sstring

    ! Do we have the Zero Mean pressure constraint
    call parlst_getvalue_int (rparams, 'ZMV', 'detZMV', detZMV, 0)
    
    select case (detZMV)
   
    case (1)
      ! Lumped mass matrix technique 
      ! Check if this is the first nonlinear iteration,
      !  then make lumped mass matrix on all levels
      if (inl .eq. 1) then
      
        ! Read the finite element for the Pressure
        call parlst_getvalue_string (rparams, 'MESH', 'Pelm', sstring)
        Pelm = elem_igetID(sstring)
      
       do i=NLMIN,NLMAX
         ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Setup the Mass matrix on all levels
         ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Set up a block discretisation structure that has 1 block.
         call spdiscr_initBlockDiscr (Rlevels(i)%rdiscMass,1,&
                                       Rlevels(i)%rtriangulation, rboundary)

         ! For the pressure mass matrix, we set up a discretisation
         call spdiscr_initDiscr_simple (&
                  Rlevels(i)%rdiscMass%RspatialDiscr(1),Pelm, &
                                       Rlevels(i)%rtriangulation, rboundary)
         
         call lsysbl_createMatBlockByDiscr (Rlevels(i)%rdiscMass,&
                                                     Rlevels(i)%rmatrixMass)    
         call bilf_createMatrixStructure (&
                Rlevels(i)%rdiscMass%RspatialDiscr(1),LSYSSC_MATRIX9, &
                                   Rlevels(i)%rmatrixMass%RmatrixBlock(1,1))   
         
         ! Allocate memory for the matrix
         call lsysbl_allocEmptyMatrix (Rlevels(i)%rmatrixMass,&
                                                           LSYSSC_SETM_ZERO)

         ! Build up the mass matrix
         call bma_buildMatrix (Rlevels(i)%rmatrixMass,BMA_CALC_STANDARD,&
                      ls_Mass,rcubatureInfo=Rlevels(i)%rcubatureInfo)
         
         ! Making the mass matrix lumped with,
         !  old structure    .false.
         !  diagonal matrix  .true.
         call lsyssc_lumpMatrixScalar (&
          Rlevels(i)%rmatrixMass%RmatrixBlock(1,1), LSYSSC_LUMP_DIAG,.true.)
      
       end do
      
      end if
    
      ! read the Zero Mean pressure row
      call parlst_getvalue_int (rparams, 'ZMV', 'irow', irow, 1)
    
      ! Modify the RHS here
      ! Setting a zero on the row number 'irow' of the pressure RHS vector    
      call vecfil_OneEntryZero(rrhs,3,irow)    
    
      ! Modify the pressure matrix here
      do i=NLMIN,NLMAX   
        ! Set the values of the row number 'irow' of the system matrix
        !  to the diagonal values of the lumped mass matrix
        call mmod_replaceLineByLumpedMass (Rlevels(i)%rmatrix%RmatrixBlock(3,3),&
                              irow,Rlevels(i)%rmatrixMass%RmatrixBlock(1,1))
      end do
    
    case (2)
      ! L^2_0 shifting technique
      ! Do nothing!!
   
    case (3)
      ! One pressure DOF = 0
      ! read the Zero Mean pressure row
      call parlst_getvalue_int (rparams, 'ZMV', 'irow', irow, 1)
      Irows = (/irow/)
      
      ! Modify the RHS here
      ! Setting a zero on the row number 'irow' of the pressure RHS vector    
      call vecfil_OneEntryZero(rrhs,3,irow)
     
      ! Modify the pressure matrix here
      do i=NLMIN,NLMAX   
       ! Set a '1' on the main diagonal of the row number 'irow' 
       !  of pressure block of the rmatrix on all levels
       !  and zero elsewhere on that row.
       call mmod_replaceLinesByUnitBlk (Rlevels(i)%rmatrix,3,Irows)
      end do 

    case default
      ! No zero mean pressure constraint required
        
    end select
       
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_NOE(Rlevels,rrhs,rboundary,rparams,NLMAX,NLMIN)
                                
 !<description>
    ! Set up Neumann, Outflow and Etc. BCs.
 !</description>

 !<inputoutput>
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>
 
  !<input>
    ! Level information
    integer, intent(in) :: NLMAX, NLMIN
  
    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams
    
    ! An object for saving the domain:
    type(t_boundary), intent(in) :: rboundary
    
    ! RHS vector
    ! May be required later, for instance to impose 
    ! non-homogeneous Neumann BC, which requires the RHS 
    ! modifications
    type(t_vectorBlock), intent(inout) :: rrhs
  !</input>
 
!</subroutine>

    ! Local variables
    ! Do we have outflow BC
    integer :: detOFlow
    
    ! Which edge do we apply the outflow BC
    integer :: edge

    ! An object to save the desired edge
    type(t_boundaryRegion) :: rboundaryRegion    
    
    ! Line integration rule
    integer(I32) :: ccubType
    
    ! String variable
    character(len=SYS_STRLEN) :: sstring

    ! Kinematic Viscosity
    real(DP) :: dnu
    
    ! Loop index
    integer :: i
    
    ! Bilinear form
    type(t_bilinearForm) :: rform
    
    ! Determine whether we have outflow BC
    call parlst_getvalue_int (rparams,'BCS','detOFlow',detOFlow,0)

    if (detOFlow .eq. 1) then

        ! Get the line integration rule
        call parlst_getvalue_string (rparams,'BCS','ccubType', sstring)
        ccubType = cub_igetID(sstring)
        
        ! Which edge of the triangulation do we apply the outflow BC
        call parlst_getvalue_int (rparams,'BCS','edge',edge,2)         
        ! Create a boundary region
        ! Warning: We assume that this edge is on the first boundary component!
        call boundary_createRegion(rboundary,1,edge,rboundaryRegion) 
        ! Exclude the start and end points from the outflow
        rboundaryRegion%iproperties = 2**1-2**1
        
        ! Get the kinematic viscosity: noo = 1/Re
        call parlst_getvalue_double (rparams, 'GFPROPER', 'dnu', dnu, 1.0_DP)        
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !  For the moment we ONLY support natural boundary condition on
        !  the vertical outflow boundaries of a 2D problem. So,
        !        1     du1              1     du2
        !  -P + --- * ---- = 0   (1),  --- * ---- = 0   (2)
        !        Re    dx               Re    dx
        !  is implemented.  This condition is implemented in a
        !  weak manner, using the L^2 norm of the equations over the outflow
        !  boundary. Therefore, we have the following extra functionals 
        !  added to our LSFEM:
        ! 
        !   ||-q + 1/Re dv1/dx ||^2_{\Gamma_o}  (1)
        !   || 1/Re dv2/dx ||^2_{\Gamma_o}      (2)
        !
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! First, functional (1) has 4 components in a variational 
        ! setting which comes out of the Least-squares formulation:
        !
        ! (p,q) - 1/Re (p, dv1/dx) - 1/Re (du1/dx, q) + (1/Re)^2 (du1/dx, dv1/dx)
        !  *1*           *2*                 *3*                  *4*
        !
        ! Creating entries of the bilinear form for term *1*
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_FUNC
        
        ! Here we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = 1.0_DP
        
        do i=NLMIN, NLMAX
            ! Assemble the matrix, the new matrix entries are added 
            ! to the existing entries  .false.  
            call bilf_buildMatrixScalarBdr2D (rform, ccubType, .false., &
                                  Rlevels(i)%rmatrix%RmatrixBlock(3,3), &
                                          rboundaryRegion=rboundaryRegion)
        end do

        ! Creating entries of the bilinear form for term *2*
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_FUNC
        rform%Idescriptors(2,1) = DER_DERIV_X
        
        ! Here we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = -dnu
        
        do i=NLMIN, NLMAX
            ! Assemble the matrix, the new matrix entries are added 
            ! to the existing entries  .false.  
            call bilf_buildMatrixScalarBdr2D (rform, ccubType, .false., &
                                  Rlevels(i)%rmatrix%RmatrixBlock(1,3), &
                                          rboundaryRegion=rboundaryRegion)
        end do
        
        ! Creating entries of the bilinear form for term *3*
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_FUNC
        
        ! Here we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = -dnu
        
        do i=NLMIN, NLMAX
            ! Assemble the matrix, the new matrix entries are added 
            ! to the existing entries  .false.  
            call bilf_buildMatrixScalarBdr2D (rform, ccubType, .false., &
                                  Rlevels(i)%rmatrix%RmatrixBlock(3,1), &
                                          rboundaryRegion=rboundaryRegion)
        end do


        ! Creating entries of the bilinear form for term *4*
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        
        ! Here we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = dnu**2
        
        do i=NLMIN, NLMAX
            ! Assemble the matrix, the new matrix entries are added 
            ! to the existing entries  .false.  
            call bilf_buildMatrixScalarBdr2D (rform, ccubType, .false., &
                                  Rlevels(i)%rmatrix%RmatrixBlock(1,1), &
                                          rboundaryRegion=rboundaryRegion)
        end do

        ! Now the 2nd functional (2), this has only 1 component in a variational 
        ! setting which comes out of the Least-squares formulation:
        !
        !  (1/Re)^2 (du2/dx, dv2/dx)
        !           *1*
        !
        ! Creating entries of the bilinear form for term *1*
        rform%itermCount = 1
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        
        ! Here we have constant coefficients:
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = dnu**2
        
        do i=NLMIN, NLMAX
            ! Assemble the matrix, the new matrix entries are added 
            ! to the existing entries  .false.  
            call bilf_buildMatrixScalarBdr2D (rform, ccubType, .false., &
                                  Rlevels(i)%rmatrix%RmatrixBlock(2,2), &
                                          rboundaryRegion=rboundaryRegion)
        end do

    end if
    
  end subroutine


  !****************************************************************************
 
!<subroutine>
  subroutine ls_Solver_linear(Rlevels,rrhs,rvector,rparams)
                                
 !<description>
    ! Set up a linear solver, solve the problem, release the solver.
 !</description>

 !<inputoutput>
    ! Solution Vector    
    type(t_vectorBlock), intent(inout) :: rvector
    
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>
 
  !<input>
    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams
    
     ! RHS vector
    type(t_vectorBlock), intent(inout) :: rrhs  
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

    ! Solver node and other smoothers and oreconditioners
    type(t_linsolNode), pointer :: p_rsolverNodeM
    type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditionerC
    type(t_linsolNode), pointer :: p_rsmoother, p_rpreconditioner
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo    

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(2), target :: RfilterChain
    
    ! Level info.
    integer :: NLMAX,NLMIN

    ! Loop index
    integer :: i

    ! Linear Solver parameters
    integer :: LINsolverType, LINsolver, nmaxItsLS, ioutputLevelLS
    real(DP) :: depsRelLS
    
    ! Zero Mean pressure parameters
    integer :: detZMV, irow
    
    ! Multigrid integers
    integer :: ctMG, nsmMG, nmaxItsMG, nminItsMG, ioutputLevelMG
    integer :: istoppingCriterionMG, smoothMG, PrecSmoothMG
    integer :: CGsolverMG, NItCGsolverMG, ioutputLevelCG
    
    ! Multigrid reals
    real(DP) :: depsRelMG, depsAbsMG, DampPrecSmoothMG, depsRelCGsolverMG
    
    ! Read in some parameters
    ! Level info.
    call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
    call parlst_getvalue_int (rparams, 'MESH', 'NLMIN', NLMIN, 3)
    
    ! Multigrid, cycle type
    call parlst_getvalue_int (rparams, 'MULTI', 'ctMG', ctMG, 0)
    ! Multigrid, number of pre/post smoother steps
    call parlst_getvalue_int (rparams, 'MULTI', 'nsmMG', nsmMG, 4)
    ! Multigrid, max. number of iterations
    call parlst_getvalue_int (rparams, 'MULTI', 'nmaxItsMG', nmaxItsMG, 10)
    ! Multigrid, min. number of iterations
    call parlst_getvalue_int (rparams, 'MULTI', 'nminItsMG', nminItsMG, 2)
    ! Multigrid, output level
    call parlst_getvalue_int (rparams, 'MULTI', 'ioutputLevelMG',ioutputLevelMG, 0)
    ! Multigrid, stopping criterion
    call parlst_getvalue_int (rparams, 'MULTI', 'istoppingCriterionMG', &
                                                  istoppingCriterionMG, 0)
    ! Multigrid, type of stopping criterion
    call parlst_getvalue_int (rparams, 'MULTI', 'ctMG', ctMG, 0)
    ! Multigrid, relative error
    call parlst_getvalue_double (rparams, 'MULTI', 'depsRelMG', depsRelMG, 0.001_DP)         
    ! Multigrid, absolute error
    call parlst_getvalue_double (rparams, 'MULTI', 'depsAbsMG', depsAbsMG, 0.00001_DP)


    ! Now we have to build up the level information for multigrid.
    ! Create a Multigrid-solver. Attach the filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    ! Dirichlet boundary condition filter
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
    
    ! Do we have the Zero Mean pressure constraint
    call parlst_getvalue_int (rparams, 'ZMV', 'detZMV', detZMV, 0)
    select case (detZMV)
    case (1,3)
        ! Lumped mass matrix technique .OR.
        ! One pressure DOF = 0 technique
        call parlst_getvalue_int (rparams, 'ZMV', 'irow', irow, 1)
        RfilterChain(2)%ifilterType = FILTER_ONEENTRY0
        RfilterChain(2)%iblock = 3  ! pressure block
        RfilterChain(2)%irow = irow
    case (2)
        ! L^2_0 shifting technique
        RfilterChain(2)%ifilterType = FILTER_TOL20
        RfilterChain(2)%itoL20component = 3  ! pressure block
    end select
    
    call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)

    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
    !     Set up a coarse grid solver.
    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+    
    ! The coarse grid in multigrid is always grid 1
    call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)

    ! Choose Coarse Grid Solver
    call parlst_getvalue_int (rparams, 'MULTI', 'CGsolverMG', CGsolverMG, 2)
    
    if (CGsolverMG .eq. 1) then
      ! We set up UMFPACK as coarse grid solver
      call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
      
    else
    
      ! We set up an iterative solver (The same as smoother)
      !   as coarse grid solver.
      nullify(p_rpreconditionerC)
      
      !!! Preconditioner
      call parlst_getvalue_int (rparams, 'MULTI', 'PrecSmoothMG', PrecSmoothMG, 1)
      call parlst_getvalue_double (rparams, 'MULTI', 'DampPrecSmoothMG', DampPrecSmoothMG, 1.0_DP)
      
      select case (PrecSmoothMG)
      case (1)
        ! Vanka preconditioner.
        call linsol_initVANKA (p_rpreconditionerC,DampPrecSmoothMG)      
      case (2)
        ! SSOR preconditioner.
        call linsol_initSSOR (p_rpreconditionerC, DampPrecSmoothMG, .true.)
      case (3)
        ! Jacobi preconditioner  
        call linsol_initJacobi (p_rpreconditionerC,DampPrecSmoothMG)
      case default
        ! No preconditioner is required.
      end select
      
      !!! Coarse grid solver
      call parlst_getvalue_int (rparams, 'MULTI', 'smoothMG', smoothMG, 1)
      select case (smoothMG)
      case (1)
        call linsol_initCG (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditionerC,RfilterChain)
      case (2)
        call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rpreconditionerC,RfilterChain)
      end select
      
      ! Some other coarse grid properties
      call parlst_getvalue_int (rparams, 'MULTI', 'NItCGsolverMG', NItCGsolverMG, 5)
      call parlst_getvalue_double (rparams, 'MULTI', 'depsRelCGsolverMG', depsRelCGsolverMG, 0.00001_DP)    
      call parlst_getvalue_int (rparams, 'MULTI', 'ioutputLevelCG', ioutputLevelCG, -1)
      
      ! Number of iteration of the coarse grid solver
      p_rlevelInfo%p_rcoarseGridSolver%nmaxIterations = NItCGsolverMG
      ! Relative error of the coarse grid solver
      p_rlevelInfo%p_rcoarseGridSolver%depsRel = depsRelCGsolverMG
      ! Output level of the coarse grid solver
      p_rlevelInfo%p_rcoarseGridSolver%ioutputLevel = ioutputLevelCG
               
    end if
 
    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
    ! Now set up the other levels...    
    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
    call parlst_getvalue_int (rparams, 'MULTI', 'PrecSmoothMG', PrecSmoothMG, 1)
    call parlst_getvalue_double (rparams, 'MULTI', 'DampPrecSmoothMG', DampPrecSmoothMG, 1.0_DP)    
    call parlst_getvalue_int (rparams, 'MULTI', 'smoothMG', smoothMG, 1)
    
    do i = NLMIN+1, NLMAX
    
      !!! We set up preconditioner of the smoother
      nullify(p_rpreconditioner)
   
      select case (PrecSmoothMG)
      case (1)
        ! Vanka preconditioner.
        call linsol_initVANKA (p_rpreconditioner,DampPrecSmoothMG)      
      case (2)
        ! SSOR preconditioner.
        call linsol_initSSOR (p_rpreconditioner, DampPrecSmoothMG, .true.)
      case default
        ! No preconditioner is required.
      end select
      
      !!! Smoother
      select case (smoothMG)
      case (1)
        call linsol_initCG (p_rsmoother,p_rpreconditioner,RfilterChain)
      case (2)
        call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,RfilterChain)
      end select    

      call linsol_convertToSmoother(p_rsmoother, nsmMG, 1.0_DP)

      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-NLMIN+1,p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
    end do  


    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
    ! Linear Solver Set up    
    !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
    call parlst_getvalue_int (rparams, 'LINSOL', 'LINsolverType', LINsolverType, 1)
    if (LINsolverType .eq. 1) then
    
        ! Stand-alone Multigrid solver
        ! Attach the system matrices to the solver.
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

        ! Some more solver parameters
        p_rsolverNode%nmaxIterations = nmaxItsMG
        p_rsolverNode%nminIterations = nminItsMG
        p_rsolverNode%depsRel = depsRelMG
        p_rsolverNode%depsAbs = depsAbsMG
        p_rsolverNode%istoppingCriterion = istoppingCriterionMG
        p_rsolverNode%P_RSUBNODEMULTIGRID2%ICYCLE = ctMG
        p_rsolverNode%ioutputLevel = ioutputLevelMG
        
        ! Finally solve the system
        call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)
        
        ! Check for the linear solver divergence 
        if (p_rsolverNode%DFinalDefect .gt. 1E8) then
            call output_lbrk()
            call output_line ('Linear Loop is diverged :(')
            call output_line (&
             'Try to modify the initial solution OR linear solver properties.')
            call output_lbrk()
            call sys_halt()
         end if
        
        ! Release solver data and structure
        call linsol_doneData (p_rsolverNode)
        call linsol_doneStructure (p_rsolverNode)
        
        ! Release the solver node and all subnodes attached to it (if at all):
        call linsol_releaseSolver (p_rsolverNode)
      
        ! Release temporary vector
        call lsysbl_releaseVector (rtempBlock)
    
    else
      
        ! Multigrid-preconditioned (CG\BiCGStab) solver
        ! Use our multigrid solver node, 'p_RsolverNode', as a preconditioner
        call parlst_getvalue_int (rparams, 'LINSOL', 'LINsolver', LINsolver, 1)
        
        if (LINsolver .eq. 1) then
          call linsol_initCG(p_RsolverNodeM,p_RsolverNode,RfilterChain)
        else
          call linsol_initBiCGStab(p_RsolverNodeM,p_RsolverNode,RfilterChain)  
        end if
         
        ! Attach the system matrices to the solver.
        allocate(Rmatrices(NLMIN:NLMAX))
        do i = NLMIN, NLMAX
          call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
              Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end do
        
        call linsol_setMatrices(p_RsolverNodeM,Rmatrices(NLMIN:NLMAX))

        ! We can release Rmatrices immediately -- as long as we do not
        ! release Rlevels(i)%rmatrix!
        do i=NLMIN,NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)
        
        ! Initialise structure/data of the solver. This allows the
        ! solver to allocate memory / perform some precalculation
        ! to the problem.
        call linsol_initStructure (p_rsolverNodeM, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) then
          call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
          call sys_halt()
        end if
        
        call linsol_initData (p_rsolverNodeM, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) then
          call output_line("Matrix singular!",OU_CLASS_ERROR)
          call sys_halt()
        end if

        ! Some more solver parameters
        p_rsolverNode%nmaxIterations = nmaxItsMG
        p_rsolverNode%nminIterations = nminItsMG
        p_rsolverNode%depsRel = depsRelMG
        p_rsolverNode%depsAbs = depsAbsMG
        p_rsolverNode%istoppingCriterion = istoppingCriterionMG
        p_rsolverNode%P_RSUBNODEMULTIGRID2%ICYCLE = ctMG
        p_rsolverNode%ioutputLevel = ioutputLevelMG
        
        ! Main linear solver properties
        ! Max. number of iterations
        call parlst_getvalue_int (rparams, 'LINSOL', 'nmaxItsLS', nmaxItsLS, 10)
        ! Output level
        call parlst_getvalue_int (rparams, 'LINSOL', 'ioutputLevelLS',ioutputLevelLS, 0)
        ! Relative error
        call parlst_getvalue_double (rparams, 'LINSOL', 'depsRelLS', depsRelLS, 0.001_DP)         
        
        p_rsolverNodeM%nmaxIterations = nmaxItsLS
        p_rsolverNodeM%depsRel = depsRelLS
        p_rsolverNodeM%ioutputLevel = ioutputLevelLS

        ! Finally solve the system
        call linsol_solveAdaptively (p_rsolverNodeM,rvector,rrhs,rtempBlock)
        
        ! Check for the linear solver divergence 
        if (p_rsolverNodeM%DFinalDefect .gt. 1E8) then
            call output_lbrk()
            call output_line ('Linear Loop is diverged :(')
            call output_line (&
             'Try to modify the initial solution OR linear solver properties.')
            call output_lbrk()
            call sys_halt()
         end if
        
        ! Release solver data and structure
        call linsol_doneData (p_rsolverNodeM)
        call linsol_doneStructure (p_rsolverNodeM)
        
        ! Release the solver node and all subnodes attached to it (if at all):
        call linsol_releaseSolver (p_rsolverNodeM)
      
        ! Release temporary vector
        call lsysbl_releaseVector (rtempBlock)
       
      
    end if
       
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_vec_collection(revalVectors,rvector_old,rvector_oldT)
                                
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
   ! Solution vector in the current nonliner/time iteration  
    type(t_vectorBlock), intent(inout) :: rvector_old
    
    ! Solution vector due to the time discretization
    ! This is optional because the subroutine may also be called during 
    !  the matrix assembly.
    type(t_vectorBlock), intent(inout), optional :: rvector_oldT
  !</input>
 
!</subroutine>
 
    ! The routine <verb>fev2_addVectorToEvalList</verb> allows to define
    ! the evaluation of functions and 1st/2nd derivatives.
    
    ! First we add the nonlinear vectors
    call fev2_addVectorToEvalList(revalVectors,&
         rvector_old%RvectorBlock(1),1)   ! u1       (1)
    call fev2_addVectorToEvalList(revalVectors,&
         rvector_old%RvectorBlock(2),1)   ! u2       (2)
    
    ! Then, we add the time discretization vector if it is passed
    if (present(rvector_oldT)) then
        call fev2_addVectorToEvalList(revalVectors,&
         rvector_oldT%RvectorBlock(1),1)   ! u1_T    (3)
        call fev2_addVectorToEvalList(revalVectors,&
         rvector_oldT%RvectorBlock(2),1)   ! u2_T    (4)
        call fev2_addVectorToEvalList(revalVectors,&
         rvector_oldT%RvectorBlock(3),1)   ! p_T     (5)
        call fev2_addVectorToEvalList(revalVectors,&
         rvector_oldT%RvectorBlock(4),1)   ! w_T     (6)
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
    integer :: detVJump
    integer :: detWJump
    integer :: detPJump
    real(DP) :: dJumpV, dJumpStarV, deojEdgeExpV
    real(DP) :: dJumpW, dJumpStarW, deojEdgeExpW
    real(DP) :: dJumpP
    
    ! Let's check if we realy have to set up jump stabilization
    call parlst_getvalue_int (rparams, 'JUMP', 'detVJump', detVJump, 0)
    call parlst_getvalue_int (rparams, 'JUMP', 'detWJump', detWJump, 0)    
    call parlst_getvalue_int (rparams, 'JUMP', 'detPJump', detPJump, 0) 
    
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
        rjumpStabil%ccubType = CUB_G3_1D

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
        rjumpStabil%ccubType = CUB_G3_1D

        ! Call the jump stabilisation technique for the vorticity.
        call conv_jumpStabilisation2d (rjumpStabil, CONV_MODMATRIX, &
                                            rmatrix%RmatrixBlock(4,4)) 
        
    end if    
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Vorticity jump stabilization
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    if (detPJump .eq. 1) then
        call parlst_getvalue_double (rparams, 'JUMP', 'dJumpP', &
                                                    dJumpP, 0.01_DP)
                                                    
        ! Call the jump stabilisation technique for the pressure.
        call jstab_calcReacJumpStabilisation (rmatrix%RmatrixBlock(3,3),&
                  dJumpP,1.0_DP,CUB_G3_1D,1.0_DP)        
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
    real(DP) :: Dres(4),Dresv(4),Dres_rel(4)
    integer, dimension(4) :: Cnorms
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
    Dres_rel(3) = Dres(3)/Dresv(3)
    Dres_rel(4) = Dres(4)/Dresv(4)
    
    ! Convergence check
    converged = .false.
    isum = 0   
    ! Iteration number control
    if (inl .eq. NLN_Max) then
        converged = .true.
    else
        ! Norm control
        do i=1,4
            if (Dres_rel(i) .lt. dNLEpsi) then
               isum = isum + 1
            end if
        end do
        if (isum .eq. 4) then
            converged = .true.
        end if
    end if  

    ! Divergence check
    diverged = .false.
    diverged = .not. (Dres_rel(1) .lt. 1E8 .and. Dres_rel(2) .lt. 1E8 &
                   .and. Dres_rel(3) .lt. 1E8 .and. Dres_rel(4) .lt. 1E8)

  
    ! Release the block vector
    call lsysbl_releaseVector (rdiff)


    ! Some output data
    if (inl .eq. 1) then
        call output_line ('Iter. ' //' U1 Rel. Err. ' //' U2 Rel. Err. ' &
        //' P  Rel. Err. ' //' W  Rel. Err. ')
        call output_line ('----------------------------------------------'//&
            '----------------')
        call output_line (sys_siL(inl, 5) //'  '&
        //trim(sys_sdEL(Dres_rel(1),6))//'  '&
        //trim(sys_sdEL(Dres_rel(2),6))//'  '&
        //trim(sys_sdEL(Dres_rel(3),6))//'  '&
        //trim(sys_sdEL(Dres_rel(4),6)))   
    else
        call output_line (sys_siL(inl, 5) //'  '&
        //trim(sys_sdEL(Dres_rel(1),6))//'  '&
        //trim(sys_sdEL(Dres_rel(2),6))//'  '&
        //trim(sys_sdEL(Dres_rel(3),6))//'  '&
        //trim(sys_sdEL(Dres_rel(4),6)))   
        if ( (mod(inl,10) .eq. 0) .and. (inl .ne. NLN_Max) &
            .and. (.not. converged) .and. (.not. diverged)) then
            call output_lbrk()
            call output_line ('Iter. ' &
            //' U1 Rel. Err. ' //' U2 Rel. Err. ' //' P  Rel. Err. ' &
            //' W  Rel. Err. ')
            call output_line ('--------------------------------------'//&
                '------------------------')
        end if
    end if

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_update_solution(Rlevels,converged,diverged,&
                                rvector,rvector_old,rrhs,NLMAX,NLMIN,det)
                                
 !<description>  
    ! 
 !</description>                                

 !<inputoutput>
    logical, intent(inout) :: det
    
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
      
    ! Convergence divergence parameters
    logical, intent(inout) :: converged,diverged
    
   ! Solution vectors in the current/previous nonliner iterations
   ! and the RHS vector
    type(t_vectorBlock), intent(inout) :: rvector,rvector_old,rrhs    
    
 !<\inputoutput>
 
  !<input>
    ! Nonlinear loop's maximum/current number of iterations
    integer, intent(in) :: NLMAX,NLMIN
  !</input>
 
!</subroutine>

    
    ! Local variables
    ! Loop index
    integer :: i
    
    ! Let's check the convergence or divergence status
    if ((.not. converged) .and. (.not. diverged)) then
    
      ! Copy the current solution to the old solution
      call lsysbl_copyVector (rvector,rvector_old)
      !*** Clear all data in matrix and RHS ***!
      call lsysbl_clearVector (rrhs)
      do i = NLMIN, NLMAX
        call lsysbl_clearMatrix (Rlevels(i)%rmatrix)
      end do
    else
    
      if (diverged) then
        call output_lbrk()
        call output_line ('Nonlinear Loop is diverged :(')
        call output_lbrk()
        call sys_halt()
      end if
      
      ! This will stop the nonlinear loop
      !  but since the loop is converged, continues
      !  with postprocessing. 
      det = .true.
      
    end if
 
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
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! A block vector which contains projected data
    type(t_vectorBlock) :: rprjVector

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
    integer :: detGMC, nlevels
    real(DP) :: Dfluxi, Dfluxo, Dgmc
    real(DP) :: Dfluxo5, Dfluxo1,Dfluxo15, Dfluxo2,Dfluxo21, Dfluxo22
    real(DP), dimension(2,2) :: Dcoords
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
      call parlst_getvalue_int (rparams, 'MESH', 'Vtild', Vtild, 0)
      call parlst_getvalue_int (rparams, 'MESH', 'Ptild', Ptild, 0)
      call parlst_getvalue_int (rparams, 'MESH', 'Wtild', Wtild, 0)

      if ( (Vtild .eq. 1) .or. (Ptild .eq. 1) .or. (Wtild .eq. 1)) then
              
          ! make a new discretization structure for the projected data
          ! and modify its sub-structures if required
          call spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
          
          if (Vtild .eq. 1) then
              call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
              EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(1))

              call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
              EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(2))
          endif

          if (Ptild .eq. 1) then
              call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(3), &
              EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(3))                 
          endif        

          if (Wtild .eq. 1) then
              call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(4), &
              EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(4))                 
          endif 
       
          ! Now set up a new solution vector based on this discretisation,
          ! allocate memory.
          call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)

          ! Then take our original solution vector and convert it according to the
          ! new discretisation:
          call spdp_projectSolution (rvector,rprjVector)

          ! Discretise the boundary conditions according to the new discretisation
          ! Create a t_discreteBC structure where we store all discretised boundary
          ! conditions.
          call bcasm_initDiscreteBC(rprjDiscreteBC)

          ! Prepare the discrete BCs. data
          call ls_BCs_Dirichlet_One(rprjDiscretisation,rboundary,rprjDiscreteBC)

          ! Hang the pointer into the vector.
          rprjVector%p_rdiscreteBC => rprjDiscreteBC

          ! Send the vector to the boundary-condition implementation filter.
          ! This modifies the vector according to the discrete boundary
          ! conditions.
          call vecfil_discreteBCsol (rprjVector) 
          
          select case (ExporType)
            case (0)
          
              ! Start UCD export to VTK file:
              call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.vtk')

              ! Write Pressure
              call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

              ! Write velocity field
              call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
              call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
              call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

              ! Write Vorticity
              call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)
          
            case (1)
              
              ! Start UCD export to GMV file:
              call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

              ! Write Pressure
              call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
              
              ! Write velocity field       
              call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
              call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
              call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
              
              ! Write Vorticity
              call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)         

            case (2)
              
              ! Start UCD export to GMV file, Binary GMV
              call ucd_startBGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

              ! Write Pressure
              call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
              
              ! Write velocity field       
              call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
              call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
              call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
              
              ! Write Vorticity
              call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)

            case default
              call output_line ('Invalid visualisation output type.', &
                                OU_CLASS_ERROR,OU_MODE_STD,'ls_postprocess')
              call sys_halt()
          end select
                    
          ! Release the temporary projected vector
          call lsysbl_releaseVector (rprjVector)

          ! Release our discrete version of the projected boundary conditions
          call bcasm_releaseDiscreteBC (rprjDiscreteBC)
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
              call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

              ! Write velocity field
              call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
              call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
              call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

              ! Write Vorticity
              call lsyssc_getbase_double (rvector%RvectorBlock(4),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'w',UCD_VAR_STANDARD,p_Ddata)
          
          case (1)
              ! Start UCD export to GMV file:
              call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

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
          
          case (2)
              ! Start UCD export to GMV file, Binary GMV
              call ucd_startBGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                  trim(sucddir)//'/nslsfem.'//trim(sys_si0L(itime,4))//'.gmv')  

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
            
            case default
              call output_line ('Invalid visualisation output type.', &
                                OU_CLASS_ERROR,OU_MODE_STD,'ls_postprocess')
              call sys_halt()
          
          end select  
      
      end if ! end of real or projected data condition
        
      ! Write the file to disc, that is it.
      call ucd_write (rexport)
      call ucd_release (rexport)
!     
!      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!      ! Calculate drag-/lift coefficients on the 2nd boundary component.
!      ! This is for the benchmark problem: flow around cylinder!
!      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
!
!      ! Determine whether to calculate flow around cylinder parameters or not
!      call parlst_getvalue_int (rparams, 'POST', 'LiftDragASO', LiftDragASO, 0)
!      
!      if (LiftDragASO .eq. 1) then
!        
!          call boundary_createRegion (rboundary,2,0, rboundaryRegion)
!          rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
!          
!          ! Kinematic viscosity, noo = 1/Re
!          call parlst_getvalue_double (rparams, 'GFPROPER', 'dnu', dnu, 1.0_DP)  
!                   
!          call ppns2D_bdforces_uniform (rvector,rboundaryRegion,Dforces,CUB_G3_1D,&
!                  dnu,2.0_DP)
!          call output_lbrk()
!          call output_line ('Body forces (Line Integration)')
!          call output_line ('-----------------------------')
!          call output_line ('Drag/Lift')
!          call output_line (trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
!                            //trim(sys_sdEP(Dforces(2),15,6)))
!
!          call ppns2D_bdforces_uniform (rvector,rboundaryRegion,Dforces,CUB_G3X3,&
!                  dnu,2.0_DP)
!          call output_lbrk()
!          call output_line ('Body forces (Volume Integration)')
!          call output_line ('--------------------------------')
!          call output_line ('Drag/Lift')
!          call output_line (trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
!                            //trim(sys_sdEP(Dforces(2),15,6)))
!
!      end if      
!      
!      
      
    else
      ! ##########################################
      ! This is a full post-processing subroutine 
      ! ##########################################   
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write the final result in a data file. This can be later read as an
      !  initial solution for the non-linear loop.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! Determine whether to write the final solution in a data file
      call parlst_getvalue_int (rparams, 'POST', 'detWriteResult', detWriteResult, 0)   
      if (detWriteResult .eq. 1) then
      
          ! Write the final solution on level NLMAX to a data file
          call parlst_getvalue_string (rparams, 'POST', &
                 'sFilenamePathResult',sfile, "", bdequote=.true.)
            
          call vecio_writeBlockVectorHR (rvector, 'SOLUTION', .true.,&
                                             0, sfile, '(E22.15)')
      end if

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
      !				     \int_{\Gamma_i}{n.v}
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! Check whether to calculate the GMC or not.
      call parlst_getvalue_int (rparams, 'POST', 'detGMC',detGMC,0)

      if (detGMC .eq. 1) then
      
        ! Obtain the number of elements
        call parlst_getvalue_int (rparams, 'POST', 'nlevels',nlevels,0)
        
        if (nlevels == 0) then
          ! Calculate the number of elements automatically
          
          ! Input line coordinates
          call parlst_getvalue_string (rparams, 'POST', 'inputGMC',sstring)
          read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
                             Dcoords(1,2), Dcoords(2,2)
          ! Input line flux
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                                     Dcoords(1:2,2),Dfluxi)
          call output_lbrk()
          call output_line ('flux input(%)-0.082')
          call output_line ('-------------------')
          call output_line (trim(sys_sdEP(Dfluxi,17,10)))
          
          ! Better to use the exact value of the inflow fluxes, rather than
          !  calculating it numerically
          ! Flow Around Cylinder
           Dfluxi = -0.082_DP
          ! Poiseuelle Flow
          ! Dfluxi = -1.0_DP/6.0_DP
          
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
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))       
          

          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 2.15_DP
          Dcoords(1,2) = Dcoords(1,2) + 2.15_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                                     Dcoords(1:2,2),Dfluxo5)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))      
        
        else
          ! Calculate the number of elements from nlevels
          
          ! Input line coordinates
          call parlst_getvalue_string (rparams, 'POST', 'inputGMC',sstring)
          read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
                             Dcoords(1,2), Dcoords(2,2)
          ! Input line flux
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxi,nlevels=nlevels)
          call output_lbrk()
          call output_line ('flux input(%)-0.082')
          call output_line ('-------------------')
          call output_line (trim(sys_sdEP(Dfluxi,17,10)))
          
          ! Better to use the exact value of the inflow fluxes, rather than
          !  calculating it numerically
          ! Flow Around Cylinder
          Dfluxi = -0.082_DP
          ! Poiseuelle Flow
          ! Dfluxi = -1.0_DP/6.0_DP
          
          Dfluxi = abs(Dfluxi)
          
          ! Output line coordinates
          call parlst_getvalue_string (rparams, 'POST', 'outputGMC',sstring)
          read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
                             Dcoords(1,2), Dcoords(2,2)
          ! Output line flux
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                      Dcoords(1:2,2),Dfluxo,nlevels=nlevels)
          
          Dfluxo = abs(Dfluxo)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))       
          

          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
          Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))       
          

          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
          Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))           
          

          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
          Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))           
          

          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
          Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))
          

          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 0.1_DP
          Dcoords(1,2) = Dcoords(1,2) + 0.1_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))           
          
          
          ! Output line flux
          Dcoords(1,1) = Dcoords(1,1) + 0.1_DP
          Dcoords(1,2) = Dcoords(1,2) + 0.1_DP
          call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
          
          Dfluxo5 = abs(Dfluxo5)
          
          ! The GMC is then calculated as
          Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                                                     
          
          ! Print the GMC value
          call output_lbrk()
          call output_line ('Global Mass Conservation(%)')
          call output_line (&
          '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
          call output_line (trim(sys_sdEP(Dgmc,16,7)))           
          
        end if
          
      end if
      
    end if
    
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_cleanup(rvector,rvector_old,rvector_oldT,rrhs,Rlevels,&
                      p_rtempVectorSc,rboundary,NLMAX,NLMIN,rparams)
 
 !<description>  
    ! Release all the memory used in our calculations.
 !</description> 

  !<input>
    integer, intent(in) :: NLMAX, NLMIN

    ! All parameters in LSFEM solver
    type(t_parlist), intent(in) :: rparams
  !</input>  

  !<inputoutput>
    ! An array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels
      
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! A couple of block vectors.
    type(t_vectorBlock) :: rvector,rvector_old,rvector_oldT,rrhs
    
    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar), pointer :: p_rtempVectorSc
  !</inputoutput>

!</subroutine>

    ! Local variables
    ! Loop index
    integer :: i
    
    ! Determine if zero mean pressure constraint was active 
    integer :: detZMV 
    
    
    ! Zero mean pressure active or not
    call parlst_getvalue_int (rparams, 'ZMV', 'detZMV', detZMV, 0)
    
    ! Now, clean up so that all the memory is available again.    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rvector_old)
    call lsysbl_releaseVector (rvector_oldT)
    call lsysbl_releaseVector (rrhs)
    
    do i = NLMAX, NLMIN, -1
        call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
        if (detZMV .eq. 1) call lsysbl_releaseMatrix (Rlevels(i)%rmatrixMass)
        if (i .ne. NLMAX) then
          call lsysbl_releaseVector (Rlevels(i)%rtempVector)
        end if
    end do

    ! Release temporary scalar vector
    if (p_rtempVectorSc%NEQ .ne. 0) call lsyssc_releaseVector(p_rtempVectorSc)

    ! Release the cubature formula
    do i = NLMAX, NLMIN, -1
        call spdiscr_releaseCubStructure (Rlevels(i)%rcubatureInfo)
    end do
    
    ! Release our discrete version of the boundary conditions
    do i = NLMAX, NLMIN, -1    
        call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
    end do

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    do i = NLMAX, NLMIN, -1
        call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
        if (detZMV .eq. 1) call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscMass)
    end do
    
    ! Release the triangulation.
    do i = NLMAX, NLMIN, -1
        call tria_done (Rlevels(i)%rtriangulation)
    end do
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

  end subroutine
 
 
  ! ***************************************************************************

!<subroutine>
  subroutine ls_ns2DT_Matrix(RmatrixData,rassemblyData,rmatrixAssembly,&
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
      
    ! Time-dependent parameters
    real(DP) :: Dtstp, Dtheta, Dtsc, Dtt   
    
    ! A handle to the element volumes
    integer :: ihandle, i, ielreal
  
    ! Velocity values/derivatives in cubature points 
    real(DP) :: dU, dV, dUx, dUy, dVx, dVy
      
    ! Viscosity
    dnu = rcollection%DquickAccess(1)
    
    ! Continuity equation scaling factor
    !  default = 1.0_DP    no scaling
    alpha = rcollection%DquickAccess(3)

    ! Linearization Scheme
    Dfpn = rcollection%DquickAccess(4)
  
    ! The value of \theta in our temporal discretization
    Dtheta = rcollection%DquickAccess(5)
    
    ! The time step size = \nabla t
    Dtstp = rcollection%DquickAccess(6)
    
    ! The time-dependent scaling factor = T/nu
    Dtsc = rcollection%DquickAccess(7)
    
    ! A very common expression is \theta * \nabla t
    Dtt = Dtheta * Dtstp
  
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
            dval = p_DcubWeight(icubp,iel) * (  Dtsc*Dtt*Dtt*Dphy*(dU*dU * dbasJx*dbasIx + &
                   dU*dV * dbasJx*dbasIy + dV*dU * dbasJy*dbasIx + &
                   dV*dV * dbasJy*dbasIy + Dfpn*(dU*dUx * dbasJx*dbasI + &
                   dV*dUx * dbasJy*dbasI + dUx*dU * dbasJ*dbasIx + &
                   dUx*dV * dbasJ*dbasIy + dUx*dUx * dbasJ*dbasI + &
                   dVx*dVx * dbasJ*dbasI)) + Dadn(iel)*(dbasJy*dbasIy + alpha*dbasJx*dbasIx)  )
            dval = dval + p_DcubWeight(icubp,iel) * Dtsc*Dphy*(   dbasJ*dbasI + Dtt*(  &
                    dbasJ*dU*dbasIx + dbasJ*dV*dbasIy + dbasJx*dU*dbasI + dbasJy*dV*dbasI +&
                    Dfpn*(dbasJ*dUx*dbasI + dbasJ*dUx*dbasI)  )   )
                   
            p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
            
            ! A22
            dval = p_DcubWeight(icubp,iel) * (  Dtsc*Dtt*Dtt*Dphy*(dU*dU * dbasJx*dbasIx + &
                   dU*dV * dbasJx*dbasIy + dV*dU * dbasJy*dbasIx + &
                   dV*dV * dbasJy*dbasIy + Dfpn*(dU*dVy * dbasJx*dbasI + &
                   dV*dVy * dbasJy*dbasI + dVy*dU * dbasJ*dbasIx + &
                   dVy*dV * dbasJ*dbasIy + dUy*dUy * dbasJ*dbasI + &
                   dVy*dVy * dbasJ*dbasI)) + Dadn(iel)*(alpha*dbasJy*dbasIy + dbasJx*dbasIx)  )
            dval = dval + p_DcubWeight(icubp,iel) * Dtsc*Dphy*(   dbasJ*dbasI + Dtt*(  &
                    dbasJ*dU*dbasIx + dbasJ*dV*dbasIy + dbasJx*dU*dbasI + dbasJy*dV*dbasI +&
                    Dfpn*(dbasJ*dVy*dbasI + dbasJ*dVy*dbasI)  )   )
            
            p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
            
            ! A12
            dval = p_DcubWeight(icubp,iel) * (  Dtsc*Dtt*Dtt*Dfpn*Dphy*(dU*dVx * dbasJx*dbasI + &
                   dV*dVx * dbasJy*dbasI + dUy*dU * dbasJ*dbasIx + &
                   dUy*dV * dbasJ*dbasIy + dUy*dUx * dbasJ*dbasI + &
                   dVy*dVx * dbasJ*dbasI) + Dadn(iel)*(-dbasJx*dbasIy + alpha*dbasJy*dbasIx)  )
            dval = dval + p_DcubWeight(icubp,iel) * Dtsc*Dphy*Dtt*Dfpn*(  dbasJ*dVx*dbasI + dbasJ*dUy*dbasI  )                   
                             
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
            dval = p_DcubWeight(icubp,iel) * Dtsc*Dtt*Dtt*Dphy*(  dU * dbasJx*dbasIx + &
                   dV * dbasJx*dbasIy + Dfpn*(dUx * dbasJx*dbasI + dVx * dbasJy*dbasI)  )
            dval = dval + p_DcubWeight(icubp,iel) * Dtsc*Dtt*Dphy*(dbasJx*dbasI)
            p_DlocalMatrixA13(jdofe,idofe,iel) = p_DlocalMatrixA13(jdofe,idofe,iel) + dval
            
            ! A31
            p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
                     
            ! A23
            dval = p_DcubWeight(icubp,iel) * Dtsc*Dtt*Dtt*Dphy*(  dU * dbasJy*dbasIx + &
                   dV * dbasJy*dbasIy + Dfpn*(dUy * dbasJx*dbasI + dVy * dbasJy*dbasI)  )
            dval = dval + p_DcubWeight(icubp,iel) * Dtsc*Dtt*Dphy*(dbasJy*dbasI)
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
            dval = p_DcubWeight(icubp,iel) * (    Dtsc*Dtt*Dtt*dnu*Dphy*(  dU * dbasJy*dbasIx + &
                   dV * dbasJy*dbasIy + Dfpn*(dUx * dbasJy*dbasI - dVx * dbasJx*dbasI)  ) + &
                   Dadn(iel)*(dbasJ*dbasIy)    )
            dval = dval + p_DcubWeight(icubp,iel) * Dtsc*Dtt*dnu*Dphy*(dbasJy*dbasI)
            p_DlocalMatrixA14(jdofe,idofe,iel) = p_DlocalMatrixA14(jdofe,idofe,iel) + dval

            ! A41
            p_DlocalMatrixA41(idofe,jdofe,iel) = p_DlocalMatrixA41(idofe,jdofe,iel) + dval

            ! A24
            dval = p_DcubWeight(icubp,iel) * (    Dtsc*Dtt*Dtt*dnu*Dphy*(  -dU * dbasJx*dbasIx - &
                   dV * dbasJx*dbasIy + Dfpn*(dUy * dbasJy*dbasI - dVy * dbasJx*dbasI)  ) - &
                   Dadn(iel)*(dbasJ*dbasIx)    )
            dval = dval - p_DcubWeight(icubp,iel) * Dtsc*Dtt*dnu*Dphy*(dbasJx*dbasI)
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
            dval = p_DcubWeight(icubp,iel) * Dtsc*Dtt*Dtt*Dphy*(dbasJx*dbasIx + dbasJy*dbasIy)
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
            dval = Dtsc*Dtt*Dtt*dnu * Dphy* p_DcubWeight(icubp,iel) * ( dbasJy*dbasIx - dbasJx*dbasIy )
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
            dval = p_DcubWeight(icubp,iel) * ( Dtsc*Dtt*Dtt*dnu*dnu*Dphy*(dbasJy*dbasIy + dbasJx*dbasIx) +&
                   Dadn(iel)*(dbasJ*dbasI) )
            p_DlocalMatrixA44(jdofe,idofe,iel) = p_DlocalMatrixA44(jdofe,idofe,iel) + dval
            
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel    
  
  end subroutine

  !****************************************************************************


!<subroutine>
  subroutine ls_ns2DT_rhs(rvectorData,rassemblyData,rvectorAssembly,&
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
    real(DP) :: dbasI,dbasIx,dbasIy, dval1, dval2, dval3, dval4, dnu
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
    real(DP), dimension(:,:), pointer :: p_DlocalVector3, p_DlocalVector4
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest1,p_DbasTest3,p_DbasTest4
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData3
    type(t_bmaVectorData), pointer :: p_rvectorData4


    ! Known velocity data, nonlinear deferred velocities
    real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

    ! Time discretization deferred vectors
    real(DP), dimension(:,:,:), pointer :: p_Du3,p_Du4,p_Du5,p_Du6

    ! Velocity values/derivatives in cubature points,
    !  nonlinear deferred velocities
    real(DP) :: dU, dV, dUx, dUy, dVx, dVy

    ! Time discretization deferred vectors
    real(DP) :: dUt, dVt, dUxt, dUyt, dVxt, dVyt
    real(DP) :: dPxt, dPyt, dWxt, dWyt

    real(DP) :: Dphy, Dfpn
    
    ! Time-dependent parameters
    real(DP) :: Dtheta, Dtsc, Dtt, Dtt1, Dttt, Dtstp
    
    ! Viscosity
    dnu = rcollection%DquickAccess(1)

    ! Assigne the value of the Physical scaling
    ! This value is set in subroutine 'ls_initParams' 
    ! If no scaling is required, this value is set to 1.0_DP
    Dphy = rcollection%DquickAccess(2)
  
    ! Linearization Scheme
    Dfpn = rcollection%DquickAccess(4)
  
    ! The value of \theta in our temporal discretization
    Dtheta = rcollection%DquickAccess(5)
    
    ! The time step size = \nabla t
    Dtstp = rcollection%DquickAccess(6)
    
    ! The time-dependent scaling factor = T/nu
    Dtsc = rcollection%DquickAccess(7)
    
    ! A very common expression is
    !  (\theta * \nabla t) * T/nu * 1/(nu^2)
    Dtt = Dtheta * Dtstp * Dtsc * Dphy
    
    ! A very common expression is
    !  ((1 - \theta) * \nabla t) * T/nu * 1/(nu^2)
    Dtt1 = (1.0_DP-Dtheta) * Dtstp * Dtsc * Dphy
    
    ! A very common expression is
    !  (\theta * \nabla t) * ((1 - \theta) * \nabla t) * T/nu * 1/(nu^2)
    Dttt = Dtheta * Dtstp * (1.0_DP-Dtheta) * Dtstp * Dtsc * Dphy
  
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
    
    
    ! Get the diferred velocity field from the parameters
    p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
    p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
    
    ! Also, get the time discretization diferred vectors
    p_Du3 => revalVectors%p_RvectorData(3)%p_Ddata
    p_Du4 => revalVectors%p_RvectorData(4)%p_Ddata
    p_Du5 => revalVectors%p_RvectorData(5)%p_Ddata
    p_Du6 => revalVectors%p_RvectorData(6)%p_Ddata
    
    
    
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
        
        ! Now the time discretization deferred quantities
        dUt = p_Du3(icubp,iel,DER_FUNC)
        dVt = p_Du4(icubp,iel,DER_FUNC)
        dUxt = p_Du3(icubp,iel,DER_DERIV2D_X)
        dVxt = p_Du4(icubp,iel,DER_DERIV2D_X)
        dUyt = p_Du3(icubp,iel,DER_DERIV2D_Y)
        dVyt = p_Du4(icubp,iel,DER_DERIV2D_Y)
        
        dPxt = p_Du5(icubp,iel,DER_DERIV2D_X)
        dPyt = p_Du5(icubp,iel,DER_DERIV2D_Y)
        
        dWxt = p_Du6(icubp,iel,DER_DERIV2D_X)
        dWyt = p_Du6(icubp,iel,DER_DERIV2D_Y)
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData1%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest1(idofe,DER_FUNC,icubp,iel)
          dbasIx = p_DbasTest1(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTest1(idofe,DER_DERIV2D_Y,icubp,iel)
                              
          ! Values of the velocity RHS for the X1 and X2 component
          dval1 = 0.0_DP + Dtheta*Dtstp*Dtt*Dfpn*(   (dU*dUx + dV*dUy) * dU * dbasIx + &
                  (dU*dUx + dV*dUy) * dV* dbasIy + &
                  (dU*dUx + dV*dUy) * dUx * dbasI + &
                  (dU*dVx + dV*dVy) * dVx * dbasI   )
          dval1 = dval1 + Dfpn*Dtt*((dU*dUx + dV*dUy) * dbasI) - &
                  Dtt1*dnu*dWyt*dbasI - Dttt*dnu*dWyt*(dU*dbasIx+dV*dbasIy)-&
                  Dttt*dnu*Dfpn*(dWyt*dbasI*dUx - dWxt*dbasI*dVx) - &
                  Dtt1*(dUt*dUxt+dVt*dUyt)*dbasI - &
                  Dttt*(dUt*dUxt+dVt*dUyt)*(dU*dbasIx+dV*dbasIy) - &
                  Dttt*Dfpn*((dUt*dUxt+dVt*dUyt)*dbasI*dUx + (dUt*dVxt+dVt*dVyt)*dbasI*dVx) -&
                  Dtt1*dPxt*dbasI - Dttt*dPxt*(dU*dbasIx+dV*dbasIy) - &
                  Dttt*Dfpn*(dPxt*dbasI*dUx+dPyt*dbasI*dVx) + Dtsc*dUt*dbasI + &
                  Dtt*dUt*(dU*dbasIx+dV*dbasIy) + Dtt*Dfpn*(dUt*dbasI*dUx+dVt*dbasI*dVx)
                  
          dval2 = 0.0_DP + Dtheta*Dtstp*Dtt*Dfpn*(   (dU*dVx + dV*dVy) * dU * dbasIx + &
                  (dU*dVx + dV*dVy) * dV * dbasIy + &
                  (dU*dUx + dV*dUy) * dUy * dbasI + &
                  (dU*dVx + dV*dVy) * dVy * dbasI   )
          dval2 = dval2 + Dfpn*Dtt*((dU*dVx + dV*dVy) * dbasI) + &
                  Dtt1*dnu*dWxt*dbasI + Dttt*dnu*dWxt*(dU*dbasIx+dV*dbasIy)-&
                  Dttt*dnu*Dfpn*(dWyt*dbasI*dUy - dWxt*dbasI*dVy) - &
                  Dtt1*(dUt*dVxt+dVt*dVyt)*dbasI - &
                  Dttt*(dUt*dVxt+dVt*dVyt)*(dU*dbasIx+dV*dbasIy) - &
                  Dttt*Dfpn*((dUt*dUxt+dVt*dUyt)*dbasI*dUy + (dUt*dVxt+dVt*dVyt)*dbasI*dVy) -&
                  Dtt1*dPyt*dbasI - Dttt*dPyt*(dU*dbasIx+dV*dbasIy) - &
                  Dttt*Dfpn*(dPxt*dbasI*dUy+dPyt*dbasI*dVy) + Dtsc*dVt*dbasI + &
                  Dtt*dVt*(dU*dbasIx+dV*dbasIy) + Dtt*Dfpn*(dUt*dbasI*dUy+dVt*dbasI*dVy)
                  
          ! Multiply the values of the basis functions by
          ! the cubature weight and sum up into the local vectors.
          p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval1
          p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval2
          
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
        
        ! Now the time discretization deferred quantities
        dUt = p_Du3(icubp,iel,DER_FUNC)
        dVt = p_Du4(icubp,iel,DER_FUNC)
        dUxt = p_Du3(icubp,iel,DER_DERIV2D_X)
        dVxt = p_Du4(icubp,iel,DER_DERIV2D_X)
        dUyt = p_Du3(icubp,iel,DER_DERIV2D_Y)
        dVyt = p_Du4(icubp,iel,DER_DERIV2D_Y)
        
        dPxt = p_Du5(icubp,iel,DER_DERIV2D_X)
        dPyt = p_Du5(icubp,iel,DER_DERIV2D_Y)
        
        dWxt = p_Du6(icubp,iel,DER_DERIV2D_X)
        dWyt = p_Du6(icubp,iel,DER_DERIV2D_Y)
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData3%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasIx = p_DbasTest3(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTest3(idofe,DER_DERIV2D_Y,icubp,iel)
                              
          ! Values of the velocity RHS for pressure
          dval3 = 0.0_DP + Dtheta*Dtstp*Dtt*Dfpn*(  (dU*dUx + dV*dUy) * dbasIx + &
                  (dU*dVx + dV*dVy) * dbasIy  )
          dval3 = dval3 - Dttt*dnu*(dWyt*dbasIx-dWxt*dbasIy) - &
                  Dttt*((dUt*dUxt+dVt*dUyt)*dbasIx+(dUt*dVxt+dVt*dVyt)*dbasIy)-&
                  Dttt*(dPxt*dbasIx+dPyt*dbasIy) + Dtt*(dUt*dbasIx+dVt*dbasIy)               
                  
          ! Multiply the values of the basis functions by
          ! the cubature weight and sum up into the local vectors.
          p_DlocalVector3(idofe,iel) = p_DlocalVector3(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval3
          
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
        
        ! Now the time discretization deferred quantities
        dUt = p_Du3(icubp,iel,DER_FUNC)
        dVt = p_Du4(icubp,iel,DER_FUNC)
        dUxt = p_Du3(icubp,iel,DER_DERIV2D_X)
        dVxt = p_Du4(icubp,iel,DER_DERIV2D_X)
        dUyt = p_Du3(icubp,iel,DER_DERIV2D_Y)
        dVyt = p_Du4(icubp,iel,DER_DERIV2D_Y)
        
        dPxt = p_Du5(icubp,iel,DER_DERIV2D_X)
        dPyt = p_Du5(icubp,iel,DER_DERIV2D_Y)
        
        dWxt = p_Du6(icubp,iel,DER_DERIV2D_X)
        dWyt = p_Du6(icubp,iel,DER_DERIV2D_Y)        
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData4%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasIx = p_DbasTest4(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTest4(idofe,DER_DERIV2D_Y,icubp,iel)
                              
          ! Values of the velocity RHS for vorticity
          dval4 = 0.0_DP + Dtheta*Dtstp*Dtt*Dfpn*dnu * (  (dU*dUx + dV*dUy) * dbasIy - &
                  (dU*dVx + dV*dVy) * dbasIx  )
          dval4 = dval4 - Dttt*dnu*dnu*(dWyt*dbasIy+dWxt*dbasIx) - &
                  Dttt*dnu*((dUt*dUxt+dVt*dUyt)*dbasIy-(dUt*dVxt+dVt*dVyt)*dbasIx)-&
                  Dttt*dnu*(dPxt*dbasIy-dPyt*dbasIx) + Dtt*dnu*(dUt*dbasIy-dVt*dbasIx)
                  
          ! Multiply the values of the basis functions by
          ! the cubature weight and sum up into the local vectors.
          p_DlocalVector4(idofe,iel) = p_DlocalVector4(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval4
          
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

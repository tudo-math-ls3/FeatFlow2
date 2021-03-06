!##############################################################################
!# ****************************************************************************
!# <name> LS_SVP_NNew_MG2D </name>
!# ****************************************************************************
!# <purpose>   
!# This module solves the 2D non-Newtonian fluid flows using LSFEM.
!#
!# The second-order elliptic non-Newtonian fluid flow equations
!# are reformulated into first-order system of equations using the
!# definition of the stresses:
!#   stress: -pI + \nu(Dii(u)) \nabla(\bu)
!# The resulting system in this case is Stress-Velocity-Pressure (SVP).
!# The problem is solved in a coupled manner for the solution of:
!#   - velocity component  u1, u2 (in 2D)
!#   - pressure            p
!#   -stress component     s1, s2, s3 (in 2D)
!# variables.
!#
!# The viscosity is not constant and depends on the second invariant, Dii(u),
!# of the deformation rate tensor, D(u), for instance through power law model.
!# Here we have:
!#    D(u) = 1/2 ( grad(u) + grad(u)^T )
!#    Dii(u) = 1/2 ( 2D(u):2D(u) ) = 1/2 tr(Dii(u)^2)
!#  
!#  Viscosity defines with the POWER LAW:
!#                               r/2-1
!#    { \nu(Dii(u)) = \nu_0 Dii(u)^     
!#    { where  ( r>1 , \nu_0>0 )
!# 
!#  and with the CARREAU LAW:
!#                                                              r/2-1
!#     { \nu(Dii(u)) = \nu_inf + (\nu_0-\nu_inf) (1+\landa Dii(u))^     
!#     { where  ( \landa>0, r>1 , \nu_0>\nu_inf>=0 )
!#
!# The nonlinear terms are first linearized using Newton/Fixed-point method.
!# The LSFEM formulation then applied which yiedls a symmetric-
!# positive definite linear system of equations. This routine uses the
!# multigrid-preconditioned CG as linear solver.
!# The discretisation uses the block assembly method to evaluate
!# the system matrices all-in-one.
!#
!# </purpose>
!#
!# Author:    Masoud Nickaeen
!# First Version: Apr  22, 2013
!# Last Update:   Jun  25, 2013
!# 
!##############################################################################

module LS_SVP_NNew_MG2D

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
  use statistics
  use boundaryintegral
  use basicgeometry
  use perfconfig
  use elementpreprocessing
  use transformation
  
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
  subroutine ls_svp_nn_mg2d
  
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
  
  ! A couple of block vectors. These will be filled
  ! with data for the linear solver.
  type(t_vectorBlock) :: rvector,rvector_old,rrhs

  ! Defect calculation
  type(t_vectorBlock) :: rdefect 

  ! Temporary scalar vector; used for calculating the nonlinear matrix
  ! on lower levels / projecting the solution from higher to lower levels.
  type(t_vectorScalar), pointer :: p_rtempVectorSc => null()
  
  ! Max. Min. level to solve the problem
  integer :: NLMAX, NLMIN
  
  ! Loop index
  integer :: i
  
  ! Parameters used in nonliner loop
  integer :: inl,NLN_Max
  real(DP) :: dNLEpsi

  ! Damping parameter in nonlinear deferred correction loop
  real(DP) :: domega
  
  ! Convergence parameter, either by error or by NLN_Max
  logical :: converged,diverged

  ! Nonlinear loop control
  logical :: det = .false.
  
  ! Collection structure for callback routines
  type(t_collection) :: rcollection
  
  ! All parameter of the LSFEM solver
  type(t_parlist) :: rparams

  ! Timer objects for stopping time
  type(t_timer) :: rtimerTotal
  type(t_timer) :: rtimerSolver

  ! Ok, let's start.
  ! Initialise the timers by zero:
  call stat_clearTimer(rtimerTotal)
  call stat_clearTimer(rtimerSolver)

  ! Start the total-time Timer
  call stat_startTimer(rtimerTotal)
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 0)-Read all the parameters from data file and initialize the collection.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call ls_initParams(rparams,NLN_MAX,dNLEpsi,domega,rcollection)
  

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! a)-Read the domain, the mesh and refine it. All about computational GRID.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_grid(Rlevels,rboundary,rparams,rcollection,NLMAX,NLMIN)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! b)- Set up a discretisation and cubature info structure which tells 
  ! the code which finite element to use.
  ! Also, create a 6*6 block matrix structure.
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
      Rlevels(i)%rDiscreteBC,rcollection)
  end do
  
  ! Initialize the RHS and fine level solution vector
  call ls_Init_RhsndSolution(Rlevels,rrhs,&
          rvector_old, rvector,rparams,NLMAX)
  
  ! Initialize the memory required for interlevel projection of the
  !  solution vectors on all levels except the finest one.
  call ls_Init_CoarseSolution(Rlevels,rrhs,p_rtempVectorSc,NLMAX,NLMIN)
  
        
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! d)- Here the non-linear loop starts.
  ! The loop performs a maximum of NLN_Max iterations. The stopping criterion
  ! is controled with dNLEpsi.
  ! 
  ! Every loop performs the following series of operations:
  !   1- System matrix assembly (requires the evaluation of 
  !        nonlinear deferred velocities)
  !   1-1 And, jump stabilization set up if required
  !   2- RHS assembly (nonlinear deferred velocities must be released 
  !      right after this step!!)
  !   3- Boundary conditions implementation, excluding the one time
  !   calculation of descretized Dirichlet boundary conditions which is
  !   done earlier in 'ls_BCs_Dirichlet_One' subroutine
  !   4- Solver setup, solution of the final system, solver release
  !   5- Check for the non-linear loop convergence/divergence
  !   6- Update initial guess 'if (.not. converged) .and. (.not. diverged)'
  !   (all matrix and RHS vectors must be cleaned, 
  !   zero-valued, after this!!)
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Start Nonlinear Solver Timer
  call stat_startTimer(rtimerSolver)
  
  do inl=1,NLN_Max
  
    ! +++++++++++++++++++++++++
    ! 1- System Matrix Assembly
    ! +++++++++++++++++++++++++
    call ls_MatAssembly(Rlevels,rvector_old,p_rtempVectorSc,&
            rcollection,rparams,NLMIN,NLMAX,inl,0)

            
    ! +++++++++++++++
    ! 2- RHS Assembly
    ! +++++++++++++++
    call ls_RHS_Assembly(rrhs,rvector_old,rcollection,&
             Rlevels(NLMAX)%rcubatureInfo,inl) 
     

    ! ++++++++++++++++++++++
    ! 3- Boundary Conditions
    ! ++++++++++++++++++++++   
    ! 3-1 Dirichlet BCs.
    call ls_BCs_Dirichlet(Rlevels,NLMAX,NLMIN,rrhs,rvector,rvector_old)

    ! 3-2 Zero Mean Pressure constraint
    call ls_BCs_ZMP(Rlevels,rrhs,rboundary,inl,rparams,NLMAX,NLMIN,0)   
    

    ! +++++++++++++++++++
    ! 4- Solve the System
    ! +++++++++++++++++++
    ! 4-1 Calculate the linearized system defect
    call ls_Defect(Rlevels(NLMAX)%rmatrix, rvector_old,rrhs,rdefect)
    
    ! 4-2 Take care of the Newton scheme additional terms, if any! 
    call ls_MatAssembly(Rlevels,rvector_old,p_rtempVectorSc,&
          rcollection,rparams,NLMIN,NLMAX,inl,1)

    ! 4-3 Dirichlet BCs. preconditioner matrix ONLY
    call ls_BCs_Dirichlet(Rlevels,NLMAX,NLMIN)

    ! 4-4 Zero Mean Pressure constraint
    call ls_BCs_ZMP(Rlevels,rdefect,rboundary,inl,rparams,NLMAX,NLMIN,1)
 
    ! 4-5 Solve the system
    call ls_Solver_linear(Rlevels,rvector,rdefect,rparams)
    
    
    ! ++++++++++++++++++++++++++++++++++++++
    ! 5- Check for Convergence or Divergence
    ! ++++++++++++++++++++++++++++++++++++++
    call ls_con_di_verge(converged,diverged,rvector,&
          rvector_old,NLN_Max,inl,dNLEpsi,rdefect,domega)

    
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 6- Update the Initial Solution, clear matrix, vectors
    !  or Terminate the Loop
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++
    call ls_update_solution(Rlevels,converged,diverged,&
         rvector,rvector_old,rrhs,NLMAX,NLMIN,det)
    if (det) exit
    
  end do  ! inl
  
  ! Stop Nonlinear Solver Timer
  call stat_stopTimer(rtimerSolver)

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! e)- Post-processing starts here. Writing the solution to file ...
  ! export GMV/VTK files and calculate drag/lift forces ...
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_postprocess(rboundary,Rlevels(NLMAX)%rmatrix,&
     rvector,Rlevels(NLMAX)%rtriangulation,Rlevels(NLMAX)%rcubatureInfo,&
     Rlevels(NLMAX)%rdiscretisation,rparams,rcollection)
  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! f)- Clean up, free all the memory which is used for variables.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_cleanup(rvector,rvector_old,rrhs,Rlevels,p_rtempVectorSc,&
             rboundary,NLMAX,NLMIN,rparams)
  ! Stop the total-time Timer
  call stat_stopTimer(rtimerTotal)
  ! Print some statistics
  call output_lbrk ()
  call output_lbrk ()
  call output_line ("Time for Total Solution: "//&
  trim(sys_sdL(rtimerTotal%delapsedReal,10)))
  
  call output_line ("Time for Nonlinear Solver: "//&
  trim(sys_sdL(rtimerSolver%delapsedReal,10)))
  
  end subroutine


  !****************************************************************************
  
!<subroutine>
  subroutine ls_initParams(rparams,NLN_MAX,dNLEpsi,domega,rcollection)

 !<description>  
  ! In this subroutine, the collection structrure is initialized.
 !</description>        

 !<output>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(out) :: rparams 

  ! Nonlinear loop stopping criteria and damping parameter
  real(DP), intent(out) :: dNLEpsi, domega

  ! Nonlinear loop Max. number of iterations
  integer, intent(out) :: NLN_MAX
  
  ! Collection structure for callback routines  
  type(t_collection), intent(out) :: rcollection
 !</output>

 !</subroutine>

  ! Local variables  
  ! Continuity equation scaling factor
  !  default = 1.0_DP  no scaling
  real(DP) :: alpha
   
  ! Linearization method parameters
  integer :: LinTyp, FPIter
  
  ! reading data from the *.dat file 
  call parlst_init(rparams)
  call parlst_readfromfile (rparams, "./data/lssvp_NonNew.dat") 

  ! The problem to solve, by default cavity
  call parlst_getvalue_int (rparams, 'GFPROPER', 'Problem', &
  rcollection%IquickAccess(9), 0) 
   
  ! Nonlinear loop Max. number of iterations
  call parlst_getvalue_int (rparams, 'NLLOOP', 'NLN_MAX', NLN_MAX, 3)
  
  ! Nonlinear loop stopping criteria
  ! Serves as a measure of relative residual.
  call parlst_getvalue_double (rparams, 'NLLOOP', 'dNLEpsi', dNLEpsi, &
        0.001_DP)  
  
  ! Initializing the collections
  ! Scale least-squares functionals or not.
  call parlst_getvalue_double (rparams, 'GFPROPER', 'alpha', alpha, 1.0_DP)
  ! Continuity equation scaling
  !  \alpha * || \nabla \cdot \bu ||_0
  rcollection%DquickAccess(3) = alpha


  ! Physical scaling
  call parlst_getvalue_int (rparams, 'GFPROPER', 'scPhysic', &
            rcollection%IquickAccess(4), 1) 
  
  ! Linearization method
  call parlst_getvalue_int (rparams, 'NLLOOP', 'LinTyp', LinTyp, 0)
  call parlst_getvalue_int (rparams, 'NLLOOP', 'FPIter', FPIter, 0)

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

  
  ! Viscosity model, Powe Law by default
  call parlst_getvalue_int (rparams, 'GFPROPER', 'viscosity_model', &
  rcollection%IquickAccess(1), 1) 
  
  ! Viscosity model constants
  if (rcollection%IquickAccess(1) .eq. 1) then
  
  ! Power Law 
    call parlst_getvalue_double (rparams, 'GFPROPER', 'PL_nu0', &
                    rcollection%DquickAccess(17),1.0_DP)

    call parlst_getvalue_double (rparams, 'GFPROPER', 'PL_r', &
                    rcollection%DquickAccess(18),1.0_DP)
  else
  
  ! Carreau Law   
    call parlst_getvalue_double (rparams, 'GFPROPER', 'CL_nu0', &
                    rcollection%DquickAccess(17),1.0_DP)

    call parlst_getvalue_double (rparams, 'GFPROPER', 'CL_r', &
                    rcollection%DquickAccess(18),1.0_DP)  

    call parlst_getvalue_double (rparams, 'GFPROPER', 'CL_nuinf', &
                    rcollection%DquickAccess(19),1.0_DP)

    call parlst_getvalue_double (rparams, 'GFPROPER', 'CL_landa', &
                    rcollection%DquickAccess(20),1.0_DP) 
  
  end if
     
  ! Determine the RHS function in NS
  ! By default, no body force
  call parlst_getvalue_int (rparams, 'RHS', 'detRHS', &
    rcollection%IquickAccess(8), 0) 
  
  if (rcollection%IquickAccess(8) == 1) then 
    ! The constant parameter in the pressure analytic polynomial
    !     OR
    ! the coefficient of surface tension, sigma  
    call parlst_getvalue_double (rparams, 'RHS', 'dC', &
      rcollection%DquickAccess(9), 1.0_DP)
  end if

  ! The damping parameter in nonlinear deferred correction loop
  call parlst_getvalue_double (rparams, 'NLLOOP', 'omega', &
                               domega, 1.0_DP)
  
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
  !   grid coordinates and vice versa.
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
  integer(I32) :: Velm, Pelm, Selm  

  ! Type of cubature rule for numerical integration
  integer(I32) :: ccubType
  
  ! Jump stabilization parameters
  integer :: detVJump, detSJump, detPJump
  
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
  ! a block discretisation structure that specifies 6 blocks in the
  ! solution vector. Do this for all levels
  do i = NLMIN, NLMAX
    call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation,6,&
                 Rlevels(i)%rtriangulation, rboundary)
  end do

  ! rdiscretisation%RspatialDiscr is a list of scalar
  ! discretisation structures for every component of the solution vector.
  ! We have a solution vector with 6 components:
  !  Component 1 = X-velocity
  !  Component 2 = Y-velocity
  !  Component 3 = Pressure
  !  Component 4 = Stress1 --> \sigma_{xx}
  !  Component 5 = Stress2 --> \sigma_{xy}
  !  Component 6 = Stress3 --> \sigma_{yy}
  
  ! We set up one discretisation structure for the velocity...
  ! Read the finite element for velocities
  call parlst_getvalue_string (rparams, 'MESH', 'Velm', sstring)
  Velm = elem_igetID(sstring)

  ! Here we set up one discretisation structure for the 
  ! 1st component of the velocity vector
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
  ! structure. Read the finite element for Pressure
  call parlst_getvalue_string (rparams, 'MESH', 'Pelm', sstring)
  Pelm = elem_igetID(sstring)

  do i = NLMIN, NLMAX
    call spdiscr_deriveSimpleDiscrSc (&
        Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          Pelm, Rlevels(i)%rdiscretisation%RspatialDiscr(3))
  end do
    
  ! And for the Stresse, a separate discretisation structure as well.
  ! Read the finite element for the Stresses
  call parlst_getvalue_string (rparams, 'MESH', 'Selm', sstring)
  Selm = elem_igetID(sstring)  

  do i = NLMIN, NLMAX  
    call spdiscr_deriveSimpleDiscrSc (&
        Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
          Selm, Rlevels(i)%rdiscretisation%RspatialDiscr(4))
  end do
  
  ! ...and copy this structure also to the discretisation structure
  ! of the 5th and 6th components.
  do i = NLMIN, NLMAX
    call spdiscr_duplicateDiscrSc(&
         Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
           Rlevels(i)%rdiscretisation%RspatialDiscr(5))

    call spdiscr_duplicateDiscrSc(&
         Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
           Rlevels(i)%rdiscretisation%RspatialDiscr(6))
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
    ! The global system looks like this, a full 6*6 block matrix
    ! which is symmetric and positive definite.
    !
    !  ( A11 A12 A13 A14 A15 A16 )
    !  ( A21 A22 A23 A24 A25 A26 )
    !  ( A31 A32 A33 A34 A35 A36 )
    !  ( A41 A42 A43 A44 A45 A46 )
    !  ( A51 A52 A53 A54 A55 A56 )
    !  ( A61 A62 A63 A64 A65 A66 )
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
    ! Block A15
    call bilf_createMatrixStructure (&
        Rlevels(i)%rdiscretisation%RspatialDiscr(5), LSYSSC_MATRIX9, &
        Rlevels(i)%rmatrix%RmatrixBlock(1,5),&
        Rlevels(i)%rdiscretisation%RspatialDiscr(1)) 
    ! Block A16
    call bilf_createMatrixStructure (&
        Rlevels(i)%rdiscretisation%RspatialDiscr(6), LSYSSC_MATRIX9, &
        Rlevels(i)%rmatrix%RmatrixBlock(1,6),&
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
    ! Block A25,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,5),&
      Rlevels(i)%rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    ! Block A26,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,6),&
      Rlevels(i)%rmatrix%RmatrixBlock(2,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)

    ! Create the matrix structure of the Pressure. Block A33,
    ! Let's check if we have to set up jump stabilization
    ! If so, we need to define the matrix structure accordingly
    call parlst_getvalue_int (rparams, 'JUMP', 'detPJump', detPJump, 0)  

    ! Pressure jump stabilization
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
    ! Block A35,
    call bilf_createMatrixStructure (&
      Rlevels(i)%rdiscretisation%RspatialDiscr(5), LSYSSC_MATRIX9, &
       Rlevels(i)%rmatrix%RmatrixBlock(3,5), &
        Rlevels(i)%rdiscretisation%RspatialDiscr(3))
    ! Block A36,
    call bilf_createMatrixStructure (&
      Rlevels(i)%rdiscretisation%RspatialDiscr(6), LSYSSC_MATRIX9, &
       Rlevels(i)%rmatrix%RmatrixBlock(3,6), &
        Rlevels(i)%rdiscretisation%RspatialDiscr(3))

    ! Create the matrix structure of the Stress1.
    ! Block A44,
    ! Let's check if we have to set up jump stabilization
    ! If so, we need to define the matrix structure accordingly
    call parlst_getvalue_int (rparams, 'JUMP', 'detSJump', detSJump, 0)  

    ! Velocity jump stabilization
    if (detSJump .eq. 1) then   
      call bilf_createMatrixStructure (&
      Rlevels(i)%rdiscretisation%RspatialDiscr(4), LSYSSC_MATRIX9, &
      Rlevels(i)%rmatrix%RmatrixBlock(4,4),cconstrType=BILF_MATC_EDGEBASED)
    else 
      call bilf_createMatrixStructure (&
      Rlevels(i)%rdiscretisation%RspatialDiscr(4), LSYSSC_MATRIX9, &
      Rlevels(i)%rmatrix%RmatrixBlock(4,4))  
    end if

    ! Block A55,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(4,4),&
      Rlevels(i)%rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    ! Block A66,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(4,4),&
      Rlevels(i)%rmatrix%RmatrixBlock(6,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE) 
      
    ! Block A45,
    call bilf_createMatrixStructure (Rlevels(i)%rdiscretisation%RspatialDiscr(4),&
      LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(4,5))   
    ! Block A46,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(4,5),&
      Rlevels(i)%rmatrix%RmatrixBlock(4,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)       
                
    ! Block A56,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(4,5),&
      Rlevels(i)%rmatrix%RmatrixBlock(5,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)         
         
    
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
    

    ! Block A51,
    call lsyssc_transposeMatrix (&
      Rlevels(i)%rmatrix%RmatrixBlock(1,5), &
      Rlevels(i)%rmatrix%RmatrixBlock(5,1),LSYSSC_TR_STRUCTURE)
    ! Block A52,
    call lsyssc_duplicateMatrix (&
    Rlevels(i)%rmatrix%RmatrixBlock(5,1),&
    Rlevels(i)%rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    ! Block A53  
    call lsyssc_transposeMatrix (&
      Rlevels(i)%rmatrix%RmatrixBlock(3,5), &
      Rlevels(i)%rmatrix%RmatrixBlock(5,3),LSYSSC_TR_STRUCTURE)
    ! Block A54  
    call lsyssc_transposeMatrix (&
      Rlevels(i)%rmatrix%RmatrixBlock(4,5), &
      Rlevels(i)%rmatrix%RmatrixBlock(5,4),LSYSSC_TR_STRUCTURE)


    ! Block A61,
    call lsyssc_transposeMatrix (&
      Rlevels(i)%rmatrix%RmatrixBlock(1,6), &
      Rlevels(i)%rmatrix%RmatrixBlock(6,1),LSYSSC_TR_STRUCTURE)
    ! Block A62,
    call lsyssc_duplicateMatrix (&
    Rlevels(i)%rmatrix%RmatrixBlock(6,1),&
    Rlevels(i)%rmatrix%RmatrixBlock(6,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    ! Block A63  
    call lsyssc_transposeMatrix (&
      Rlevels(i)%rmatrix%RmatrixBlock(3,6), &
      Rlevels(i)%rmatrix%RmatrixBlock(6,3),LSYSSC_TR_STRUCTURE)
    ! Block A64  
    call lsyssc_transposeMatrix (&
      Rlevels(i)%rmatrix%RmatrixBlock(4,6), &
      Rlevels(i)%rmatrix%RmatrixBlock(6,4),LSYSSC_TR_STRUCTURE)
    ! Block A65,
    call lsyssc_duplicateMatrix (&
    Rlevels(i)%rmatrix%RmatrixBlock(6,4),&
    Rlevels(i)%rmatrix%RmatrixBlock(6,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    
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
  subroutine ls_BCs_Dirichlet_One(rdiscretisation,rboundary,rdiscreteBC,&
                              rcollection)
                
 !<description>  
  ! In this subroutine we discretise the boundary conditions and
  ! prepare them to be applied to the matrix/RHS/sol in the 
  ! nonlinear iteration loop.
 !</description>                

 !<output>
  ! A set of variables describing the analytic and discrete boundary
  ! conditions.
  type(t_discreteBC), intent(out), target :: rdiscreteBC
 !</output>

 !<input>
  ! An object for saving the domain:
  type(t_boundary), intent(in) :: rboundary

  ! An object specifying the discretisation.
  ! This contains also information about trial/test functions,...
  type(t_blockDiscretisation), intent(in) :: rdiscretisation 
  
  ! Collection structure for callback routines  
  type(t_collection), intent(inout) :: rcollection
 !</input>

!</subroutine>

  ! Local variables
  type(t_boundaryRegion) :: rboundaryRegion
  
  ! Fictitious BC
  integer, dimension(1) :: Iequations  

  ! Create a t_discreteBC structure where we store all discretised boundary
  ! conditions.
  call bcasm_initDiscreteBC(rdiscreteBC)

  select case (rcollection%IquickAccess(9))
  
  case (0)
    ! Reg. Cavity Flow
    ! edge 1 of boundary component 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     
    ! edge 2 of boundary component 1.
    call boundary_createregion(rboundary,1,2,rboundaryregion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newdirichletbconrealbd (rdiscretisation,1,&
                     rboundaryregion,rdiscretebc,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     

    ! Edge 1 of boundary component 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 2 of boundary component 1.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection) 
                           
  case (1)
    ! FAC, zero stress outflow 
    ! edge 1 of boundary component 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     

    ! Edge 2 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection) 
                 
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     

    ! Edge 1 of boundary component 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 2 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection) 
                               
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection) 

    ! cylinder
    ! X-velocity
    ! Edge 1 of boundary component 2. That is it.
    call boundary_createRegion(rboundary,2,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Y-velocity
    ! Edge 1 of boundary component 2. That is it.
    call boundary_createRegion(rboundary,2,1,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
  
  case (9)
    ! Fully developed power law flow
    ! edge 1 of boundary component 1.
    ! Symmetry line, the shear stress is zero here
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     
    ! edge 2 of boundary component 1.
    ! outflow, zero normal stress
    call boundary_createregion(rboundary,1,2,rboundaryregion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART
    call bcasm_newdirichletbconrealbd (rdiscretisation,4,&
                     rboundaryregion,rdiscretebc,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 3 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 4 of boundary component 1.
    ! inflow profile
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     

    ! Edge 1 of boundary component 1.
    ! Syymetric line, vertical velocity is zero
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 2 of boundary component 1.
    ! This is not required, just for better convergence!
    ! outflow
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 3 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 4 of boundary component 1.
    ! inflow profile
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
  
  case (10)
    ! Sudden expansion 
    ! edge 1 of boundary component 1.
    ! Symmetry line, the shear stress is zero here
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                     

    ! edge 2 of boundary component 1.
    ! outflow, zero normal stress
    call boundary_createregion(rboundary,1,2,rboundaryregion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART
    call bcasm_newdirichletbconrealbd (rdiscretisation,4,&
                     rboundaryregion,rdiscretebc,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 3 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 4 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 5 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,5,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
    
    ! Edge 6 of boundary component 1.
    ! inflow profile
    call boundary_createRegion(rboundary,1,6,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 1 of boundary component 1.
    ! Syymetric line, vertical velocity is zero
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    rboundaryRegion%iproperties = 2**1-2**1
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 2 of boundary component 1.
    ! This is not required, just for better convergence!
    ! outflow
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
                 
    ! Edge 3 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 4 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 5 of boundary component 1.
    ! no-slip on upper wall
    call boundary_createRegion(rboundary,1,5,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)

    ! Edge 6 of boundary component 1.
    ! inflow profile
    call boundary_createRegion(rboundary,1,6,rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHEND
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                     rboundaryRegion,rdiscreteBC,&
              getBoundaryValues_2D,rcollection=rcollection)
              
  case default
    ! Un-known problem
    call output_line ("Unknown problem.", OU_CLASS_WARNING, OU_MODE_STD, &
            "ls_BCs_Dirichlet_One")
    call sys_halt()   
  end select

!  ! Proot 2006
!  ! X-velocity
!  ! Edge 1 of boundary component 2. That is it.
!  call boundary_createRegion(rboundary,2,1,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!  call boundary_createRegion(rboundary,2,2,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!  call boundary_createRegion(rboundary,2,3,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!  call boundary_createRegion(rboundary,2,4,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!
!  ! Y-velocity
!  ! Edge 1 of boundary component 2. That is it.
!  call boundary_createRegion(rboundary,2,1,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!  call boundary_createRegion(rboundary,2,2,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!  call boundary_createRegion(rboundary,2,3,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)
!  call boundary_createRegion(rboundary,2,4,rboundaryRegion)
!  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
!                     rboundaryRegion,rdiscreteBC,&
!                     getBoundaryValues_2D)                     
                     

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Init_RhsndSolution(Rlevels,rrhs,rvector_old,&
                           rvector,rparams,NLMAX)
                
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
  real(DP) :: dinit_vect(6)

  ! Path to the data file which has the initial solution
  character(LEN=SYS_STRLEN) :: sfile, sstring, sarray
  
  ! Projection structure
  type(t_interlevelProjectionBlock) :: rprojection
  
  ! Temporary vectors
  type(t_vectorBlock), target :: rvector1,rvector2
  type(t_vectorScalar) :: rvectorTemp  
  
  ! Let's start
  ! Create a RHS and a solution vector based on the discretisation.
  ! Fill with zero.
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                                rrhs,.true.)
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                               rvector,.true.)   
  
  ! Determine how to setup initial nonlinear solution
  call parlst_getvalue_int (rparams, 'ISOLUTION', 'nlinit', nlinit, 0)  
  
  if (nlinit .eq. 0) then
    
    call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                             rvector_old,.true.)  
    
    ! Initialize the solution vector(s) with constant values
    call parlst_getvalue_string (rparams, 'ISOLUTION', 'initValues',&
                   sstring, '0.0_DP 0.0_DP 0.0_DP 0.0_DP 0.0_DP 0.0_DP')
    read (sstring,*) dinit_vect(1), dinit_vect(2), dinit_vect(3), & 
                               dinit_vect(4),dinit_vect(5),dinit_vect(6)
    
    ! Scale the sub-vectors to initialize the nonlineaer iteration loop
    call lsyssc_clearVector (rvector_old%RvectorBlock(1),dinit_vect(1))
    call lsyssc_clearVector (rvector_old%RvectorBlock(2),dinit_vect(2))
    call lsyssc_clearVector (rvector_old%RvectorBlock(3),dinit_vect(3))
    call lsyssc_clearVector (rvector_old%RvectorBlock(4),dinit_vect(4))    
    call lsyssc_clearVector (rvector_old%RvectorBlock(5),dinit_vect(5))
    call lsyssc_clearVector (rvector_old%RvectorBlock(6),dinit_vect(6))     
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
                 rcollection,rparams,NLMIN,NLMAX,inl,imode)
                
 !<description>  
  ! Initializing the solution vectors on all levels by calculating the 
  ! memory required to the interlevel projections.
 !</description>                
 
 !<input>
  ! Level info.
  integer, intent(in) :: NLMAX,NLMIN
  
  ! Current nonlinear loop iteration
  integer, intent(in) :: inl

  ! The calculation mode
  !  = 0 The original matrix
  !  = 1 The modified (Newton scheme) matrix ONLY, i.e. only the 
  !      new bilinear forms of the Newton method, if any, will
  !      be added to the matrix!
  integer, intent(in) :: imode
    
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
    ! structure which is setup earlier.
    call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                     p_rvectorFine,p_rvectorTemp)

    end if
  
    ! Preparing a collection of vectors to evaluate nonlinear terms
    call ls_vec_collection(revalVectors,p_rvectorCoarse)
    
    select case (imode)
    
    case (0)
      ! Calculate the original matrix to be used for the defect calculation
      ! Check if we need to assemble the physically weighted
      ! matrix
      if (rcollection%IquickAccess(4) .eq. 1) then
      
        ! The weighted case
        ! Assemble the whole system matrix on each level  
        call bma_buildMatrix (p_rmatrix,BMA_CALC_STANDARD,ls_svp2D_Matrix,&
           rcubatureInfo=Rlevels(ilev)%rcubatureInfo,rcollection=rcollection, &
           revalVectors=revalVectors)
      
      else
      
        ! The un-weighted case
        ! Assemble the whole system matrix on each level  
        call bma_buildMatrix (p_rmatrix,BMA_CALC_STANDARD,ls_svp2D_Matrix_un,&
           rcubatureInfo=Rlevels(ilev)%rcubatureInfo,rcollection=rcollection, &
           revalVectors=revalVectors)
      
      end if  ! Weighted\un-weighted

      ! Set up jump stabilization if there is
      call ls_jump(p_rmatrix,rparams,rcollection)
      
    case (1)
      ! Add the extra Newton scheme terms, if any!
      
      ! Linearization Scheme
      ! If it is a mixed type, check when to shift to
      ! Newton's
      if (rcollection%IquickAccess(2) .eq. 3) then
        if (inl .gt. rcollection%IquickAccess(3)) then
        ! It's time to shift to Newton's method 
        rcollection%DquickAccess(4) = 1.0_DP
        end if
      end if
      
      ! Add the extra terms?
      if (rcollection%DquickAccess(4) .eq. 1.0_DP) then
        ! Assemble the extra Newton terms on each level  
        ! Check if we need to assemble the physically weighted
        ! matrix
        if (rcollection%IquickAccess(4) .eq. 1) then
        
          ! The weighted case
          ! Assemble the whole system matrix on each level  
          call bma_buildMatrix (p_rmatrix,BMA_CALC_STANDARD,&
             ls_svp2D_NewtonMatrix,rcubatureInfo=Rlevels(ilev)%rcubatureInfo,&
             rcollection=rcollection,revalVectors=revalVectors)
        
        else
        
          ! The un-weighted case
          ! Assemble the whole system matrix on each level  
          call bma_buildMatrix (p_rmatrix,BMA_CALC_STANDARD,&
             ls_svp2D_NewtonMatrix_un,&
             rcubatureInfo=Rlevels(ilev)%rcubatureInfo,&
             rcollection=rcollection,revalVectors=revalVectors)
        
        end if  ! Weighted\un-weighted
      end if  ! Newton\fixed-point method 
    end select
    
    ! Release the vector structure used in linearization
    call fev2_releaseVectorList(revalVectors)
    
  end do
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_RHS_Assembly(rrhs,rvector_old,rcollection,rcubatureInfo,inl)
                
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
  ! RHS block vector
  type(t_vectorBlock), intent(inout) :: rrhs, rvector_old
  
  ! Collection structure for callback routines  
  type(t_collection) :: rcollection  
 !</inputoutput>
  
!</subroutine>


  ! Local variables  
  ! Collection of vectors to evaluate in nonlinear terms
  type(t_fev2Vectors) :: revalVectors
  
  ! Decide on the RHS force vectors in NS equation
  if (rcollection%IquickAccess(8) == 0) then
    ! No body force, do nothing 
    return
    
  else
    ! Body forces must be taken into account  
    ! Preparing a collection of vectors to evaluate nonlinear terms
    call ls_vec_collection(revalVectors,rvector_old)  
    
    !...
    
    ! Release the vector structure used in linearization
    call fev2_releaseVectorList(revalVectors)    
    
  end if
 
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_Dirichlet(Rlevels,NLMAX,NLMIN,rrhs,rvector,rvector_old)
                
 !<description>  
  ! Implementing BCs to the matrix and solution/RHS vectors.
 !</description>                

 !<input>
  integer, intent(in) :: NLMAX,NLMIN
 !</input>
 
 !<inputoutput>
  ! Block solution vectors 
  type(t_vectorBlock), intent(inout), optional :: rvector,rvector_old

  ! Block RHS vector
  type(t_vectorBlock), intent(inout), optional :: rrhs
   
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
  
  ! The RHS vector if passed
  if (present(rrhs)) then
    ! Assign the boundary conditions to the vectors ONLY on the finest level.
    ! The RHS vector
    call lsysbl_assignDiscreteBC(rrhs,Rlevels(NLMAX)%rdiscreteBC)
    ! Implement the filter  
    call vecfil_discreteBCrhs (rrhs)
  end if  
  
  ! The solution vectors if passed
  if (present(rvector)) then
    call lsysbl_assignDiscreteBC(rvector,Rlevels(NLMAX)%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rvector_old,Rlevels(NLMAX)%rdiscreteBC)
    call vecfil_discreteBCsol (rvector)
    call vecfil_discreteBCsol (rvector_old)
  end if

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_ZMP(Rlevels,rrhs,rboundary,inl,rparams,NLMAX,NLMIN,irhsmode)
                
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

  ! A parameter that separates the two modes of using this routine
  !  .eq. 0  the rhs vector is passed
  !  .eq. 1  the defect vector is passed
  integer, intent(in) :: irhsmode
 !</input>
 
 !<inputoutput>
  ! An array of problem levels for the multigrid solver
  type(t_level), dimension(:), pointer :: Rlevels
   
  ! RHS or defect block vector
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
    if ( (irhsmode .eq. 0) .and. (inl .eq. 1) ) then
    
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
     !  keep the old structure  .false.
     !  change the structure to a diagonal matrix  .true.
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

  case (4)
    ! Combined technique 
    ! Check if this is the first nonlinear iteration,
    !  then make lumped mass matrix ONLY on coarse grid NLMIN
    if (inl .eq. 1) then
    
     ! Read the finite element for the Pressure
     call parlst_getvalue_string (rparams, 'MESH', 'Pelm', sstring)
     Pelm = elem_igetID(sstring)
    
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Setup the Mass matrix on NLMIN
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Set up a block discretisation structure that has 1 block.
     call spdiscr_initBlockDiscr (Rlevels(NLMIN)%rdiscMass,1,&
                     Rlevels(NLMIN)%rtriangulation, rboundary)

     ! For the pressure mass matrix, we set up a discretisation
     call spdiscr_initDiscr_simple (&
          Rlevels(NLMIN)%rdiscMass%RspatialDiscr(1),Pelm, &
                     Rlevels(NLMIN)%rtriangulation, rboundary)
     
     call lsysbl_createMatBlockByDiscr (Rlevels(NLMIN)%rdiscMass,&
                           Rlevels(NLMIN)%rmatrixMass)  
     call bilf_createMatrixStructure (&
        Rlevels(NLMIN)%rdiscMass%RspatialDiscr(1),LSYSSC_MATRIX9, &
                   Rlevels(NLMIN)%rmatrixMass%RmatrixBlock(1,1))   
     
     ! Allocate memory for the matrix
     call lsysbl_allocEmptyMatrix (Rlevels(NLMIN)%rmatrixMass,&
                               LSYSSC_SETM_ZERO)

     ! Build up the mass matrix
     call bma_buildMatrix (Rlevels(NLMIN)%rmatrixMass,BMA_CALC_STANDARD,&
            ls_Mass,rcubatureInfo=Rlevels(NLMIN)%rcubatureInfo)
     
     ! Making the mass matrix lumped with,
     !  keep the old structure  .false.
     !  change the structure to a diagonal matrix  .true.
     call lsyssc_lumpMatrixScalar (&
      Rlevels(NLMIN)%rmatrixMass%RmatrixBlock(1,1), LSYSSC_LUMP_DIAG,.true.)
    
    end if
  
    ! read the Zero Mean pressure row
    call parlst_getvalue_int (rparams, 'ZMV', 'irow', irow, 1)
  
    ! Dooshvari darim Dooshvari :(
    ! Modify the pressure matrix here
    ! Set the values of the row number 'irow' of the system matrix
    !  to the diagonal values of the lumped mass matrix
    call mmod_replaceLineByLumpedMass (Rlevels(NLMIN)%rmatrix%RmatrixBlock(3,3),&
                irow,Rlevels(NLMIN)%rmatrixMass%RmatrixBlock(1,1))

  case default
    ! No zero mean pressure constraint required
    
  end select
     
  end subroutine


  !****************************************************************************
 
!<subroutine>
  subroutine ls_Solver_linear(Rlevels,rvector,rdefect,rparams)
                
 !<description>
  ! Set up a linear solver, solve the problem, release the solver.
 !</description>

 !<inputoutput>
  ! Solution Vector  
  type(t_vectorBlock), intent(inout) :: rvector

  ! Defect Vector
  type(t_vectorBlock), intent(inout) :: rdefect
  
  ! An array of problem levels for the multigrid solver
  type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>
 
  !<input>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams
  !</input>
 
!</subroutine>

  ! Local variables
  ! An array for the system matrix(matrices) during the initialisation of
  ! the linear solver.
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  
  ! A temporary vector for saving the defect
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
  type(t_filterChain), dimension(3), target :: RfilterChain, RfilterChainC
  
  ! Level info.
  integer :: NLMAX,NLMIN

  ! Loop index, filter counter
  integer :: i, nfilter

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
  real(DP) :: depsRelMG, depsAbsMG, DampPrecSmoothMG, depsRelCGsolverMG, DampSmoothMG
  
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
  nfilter = 1
  RfilterChain(nfilter)%ifilterType = FILTER_DISCBCDEFREAL
  
  ! Do we have the Zero Mean pressure constraint
  call parlst_getvalue_int (rparams, 'ZMV', 'detZMV', detZMV, 0)
  select case (detZMV)
  case (1,3)
    ! Lumped mass matrix technique .OR.
    ! One pressure DOF = 0 technique
    nfilter = nfilter + 1
    call parlst_getvalue_int (rparams, 'ZMV', 'irow', irow, 1)
    RfilterChain(nfilter)%ifilterType = FILTER_ONEENTRY0
    RfilterChain(nfilter)%iblock = 3  ! pressure block
    RfilterChain(nfilter)%irow = irow

  case (2)
    ! L^2_0 shifting technique
    nfilter = nfilter + 1
    RfilterChain(nfilter)%ifilterType = FILTER_TOL20
    RfilterChain(nfilter)%itoL20component = 3  ! pressure block

  case (4)
    ! L^2_0 shifting technique
    ! For every level except the coarse grid
    nfilter = nfilter + 1
    RfilterChain(nfilter)%ifilterType = FILTER_TOL20
    RfilterChain(nfilter)%itoL20component = 3  ! pressure block
    
    ! For coarse grid
    ! Lumped mass matrix technique
    nfilter = 1
    call parlst_getvalue_int (rparams, 'ZMV', 'irow', irow, 1)
    RfilterChainC(nfilter)%ifilterType = FILTER_ONEENTRY0
    RfilterChainC(nfilter)%iblock = 3  ! pressure block
    RfilterChainC(nfilter)%irow = irow
    
    ! Dirichlet BCs
    nfilter = nfilter + 1
    RfilterChainC(nfilter)%ifilterType = FILTER_DISCBCDEFREAL
    
    ! L^2_0 shifting technique
    nfilter = nfilter + 1
    RfilterChainC(nfilter)%ifilterType = FILTER_TOL20
    RfilterChainC(nfilter)%itoL20component = 3  ! pressure block
     
  end select
  
  call linsol_initMultigrid2 (p_rsolverNode,NLMAX-NLMIN+1,RfilterChain)

  !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
  !   Set up a coarse grid solver.
  !+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+  
  ! The coarse grid in multigrid is always grid 1
  call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)

  ! Choose Coarse Grid Solver
  call parlst_getvalue_int (rparams, 'MULTI', 'CGsolverMG', CGsolverMG, 2)
  
  if (CGsolverMG .eq. 1) then
    ! We set up UMFPACK as coarse grid solver
    if (detZMV == 4) then
      call linsol_initUMFPACK4 (p_rpreconditionerC)
      call linsol_initDefCorr (p_rlevelInfo%p_rcoarseGridSolver,&
                p_rpreconditionerC,RfilterChainC)
      p_rlevelInfo%p_rcoarseGridSolver%nmaxIterations = 1
    else
      call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)
    end if
  else
  
    ! We set up an iterative solver (The same as smoother)
    !   as coarse grid solver.
    nullify(p_rpreconditionerC)
    
    !!! Preconditioner
    call parlst_getvalue_int (rparams, 'MULTI', 'PrecSmoothMG', PrecSmoothMG, 1)
    call parlst_getvalue_double (rparams, 'MULTI', 'DampPrecSmoothMG', &
         DampPrecSmoothMG, 1.0_DP)
    
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
      call linsol_initCG (p_rlevelInfo%p_rcoarseGridSolver,&
           p_rpreconditionerC,RfilterChain)
    case (2)
      call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,&
           p_rpreconditionerC,RfilterChain)
    case (3,4)
    ! The Jacobi-type coarse grid solvers!!
    end select
    
    ! Some other coarse grid properties
    call parlst_getvalue_int (rparams, 'MULTI', 'NItCGsolverMG', &
            NItCGsolverMG, 5)
    call parlst_getvalue_double (rparams, 'MULTI', 'depsRelCGsolverMG', &
            depsRelCGsolverMG, 0.00001_DP)
    call parlst_getvalue_int (rparams, 'MULTI', 'ioutputLevelCG', &
            ioutputLevelCG, -1)
    
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
  call parlst_getvalue_double (rparams, 'MULTI', 'DampPrecSmoothMG', &
  DampPrecSmoothMG, 1.0_DP)  
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
    case (3)
      ! Jacobi preconditioner.
      call linsol_initJacobi (p_rpreconditioner, DampPrecSmoothMG)
    case default
      ! No preconditioner is required.
    end select
    
    !!! Smoother
    select case (smoothMG)
    case (1)
      call linsol_initCG (p_rsmoother,p_rpreconditioner,RfilterChain)
    case (2)
      call linsol_initBiCGStab (p_rsmoother,p_rpreconditioner,RfilterChain)
    case (3)
      call parlst_getvalue_double (rparams, 'MULTI', 'DampSmoothMG', &
        DampSmoothMG, 1.0_DP)  
      call linsol_initSSOR (p_rsmoother, DampSmoothMG, .true.)
    case (4)
      call parlst_getvalue_double (rparams, 'MULTI', 'DampSmoothMG', &
        DampSmoothMG, 1.0_DP)  
      call linsol_initJacobi (p_rsmoother, DampSmoothMG)
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
    ! This overwrites rdefect with the correction vector.
    ! Therefore, copy the rdefect for later use, cosmetic stuff :)
    call lsysbl_copyVector (rdefect,rtempBlock)
    
    ! Solve the preconditioned system for the correction vector
    call linsol_precondDefect (p_rsolverNode,rdefect)

    ! Copy the correction vector to the new solution vector
    ! The solution vector will be finally updated with the addition of
    ! the previous solution to it.
     call lsysbl_copyVector (rdefect,rvector)
    
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
    
    ! Release temporary vector, BUT
    ! first retrieve the rdefect vector
    call lsysbl_copyVector (rtempBlock,rdefect)
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
    call parlst_getvalue_int (rparams, 'LINSOL', 'ioutputLevelLS',&
        ioutputLevelLS, 0)
    ! Relative error
    call parlst_getvalue_double (rparams, 'LINSOL', 'depsRelLS', &
        depsRelLS, 0.001_DP)     
    
    p_rsolverNodeM%nmaxIterations = nmaxItsLS
    p_rsolverNodeM%depsRel = depsRelLS
    p_rsolverNodeM%ioutputLevel = ioutputLevelLS

    ! Finally solve the system
    ! This overwrites rdefect with the correction vector.
    ! Therefore, copy the rdefect for later use, cosmetic stuff :)
    call lsysbl_copyVector (rdefect,rtempBlock)
    
    ! Solve the preconditioned system for the correction vector
    call linsol_precondDefect (p_rsolverNodeM,rdefect)

    ! Copy the correction vector to the new solution vector
    ! The solution vector will be finally updated with the addition of
    ! the previous solution to it.
     call lsysbl_copyVector (rdefect,rvector)
    
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
    call linsol_doneData (p_rsolverNodeM)
    call linsol_doneStructure (p_rsolverNodeM)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNodeM)
    
    ! Release temporary vector, BUT
    ! first retrieve the rdefect vector
    call lsysbl_copyVector (rtempBlock,rdefect)
    call lsysbl_releaseVector (rtempBlock)
           
  end if
     
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
  integer :: detSJump
  integer :: detPJump
  real(DP) :: dJumpV, dJumpStarV, deojEdgeExpV
  real(DP) :: dJumpS
  real(DP) :: dJumpP, dJumpStarP, deojEdgeExpP
  
  ! Let's check if we realy have to set up jump stabilization
  call parlst_getvalue_int (rparams, 'JUMP', 'detVJump', detVJump, 0)
  call parlst_getvalue_int (rparams, 'JUMP', 'detSJump', detSJump, 0)  
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
  ! Stress jump stabilization
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (detSJump .eq. 1) then
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpS', &
                          dJumpS, 0.01_DP)                      
    call jstab_reacJumpStabil2d_mod (&
        rmatrix%RmatrixBlock(4,4),dJumpS,1.0_DP,CUB_G3_1D,1.0_DP)                      

    call jstab_reacJumpStabil2d_mod (&
        rmatrix%RmatrixBlock(5,5),dJumpS,1.0_DP,CUB_G3_1D,1.0_DP) 
        
    call jstab_reacJumpStabil2d_mod (&
        rmatrix%RmatrixBlock(6,6),dJumpS,1.0_DP,CUB_G3_1D,1.0_DP)         
                                             
  end if  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pressure jump stabilization
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (detPJump .eq. 1) then
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpP', &
                          dJumpP, 0.01_DP)  
    call parlst_getvalue_double (rparams, 'JUMP', 'dJumpStarP',&
                        dJumpStarP, 0.0_DP)  
    call parlst_getvalue_double (rparams, 'JUMP', 'deojEdgeExpP',&
                        deojEdgeExpP, 2.0_DP)                            
  end if  
  
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_con_di_verge(converged,diverged,rvector,rvector_old,&
               NLN_Max,inl,dNLEpsi,rdefect,domega)
                
 !<description>  
  ! Check the nonlinear loop convergence/divergence.
  ! Print the residuals.
 !</description>                

 !<output>
  ! Convergence parameter, either by error or by NLN_Max
  logical, intent(out) :: converged,diverged
 !</output>


  !<inputoutput>
  ! Solution vectors in the current nonliner iteration and
  ! the defect vector
  type(t_vectorBlock), intent(inout) :: rvector,rdefect
  !</inputoutput>
  
  !<input>
  ! Solution vectors in the previous nonliner iteration
  type(t_vectorBlock), intent(in) :: rvector_old
  
  ! Nonlinear loop's maximum/current number of iterations
  integer, intent(in) :: NLN_Max,inl
  
  ! Nonlinear loop's stopping criteria for the norms
  real(DP), intent(in) :: dNLEpsi
  
  ! The damping parameter in nonlinear deferred correction loop
  real(DP), intent(in) :: domega
  !</input>
 
!</subroutine>
 
  ! Local variables
  real(DP) :: Dres(6),Dresv(6),Dres_rel(6), defnorm(1)
  integer, dimension(6) :: Cnorms, Cdefnorm(1)
  integer :: i,isum
  
  ! Scaling factors
  real(DP) :: cx,cy
  
  ! Difference between the current and the previous
  ! vectors in nonlinear loop
  type(t_vectorBlock) :: rdiff

  ! We update the solution vector in a 
  ! deferred correction loop
  !  x^n+1 = x^n - w A^-1 (Ax^n - b)
  cx = -domega
  call lsysbl_vectorLinearComb (rvector_old,rvector,1.0_DP,cx)
  
  ! Euclidian vector norm: (vector,vector) 0
  ! $l_2$-norm: 1/sqrt(NEQ) * (vector,vector) 2
  ! max-norm 3
  Cnorms(:) = 2
  
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
  Dres_rel(5) = Dres(5)/Dresv(5)
  Dres_rel(6) = Dres(6)/Dresv(6)
  
  ! Convergence check
  converged = .false.
  isum = 0   
  ! Iteration number control
  if (inl .eq. NLN_Max) then
    converged = .true.
  else
    ! Norm control
    do i=1,6
      if (Dres_rel(i) .lt. dNLEpsi) then
         isum = isum + 1
      end if
    end do
    if (isum .eq. 6) then
      converged = .true.
    end if
  end if  

  ! Divergence check
  diverged = .false.
  diverged = .not.( Dres_rel(1) .lt. 1E8 .and. Dres_rel(2) .lt. 1E8 &
             .and. Dres_rel(3) .lt. 1E8 .and. Dres_rel(4) .lt. 1E8  &
             .and. Dres_rel(5) .lt. 1E8 .and. Dres_rel(6) .lt. 1E8  )
 
  
  ! Calculate the L^2 norm of the defect
  Cdefnorm(1) = 2
  call lsysbl_vectorNormBlock (rdefect,Cdefnorm,defnorm)


  ! Some output data
  if (inl .eq. 1) then
      call output_line ('Iter. ' &
      //' U1 Rel. Err. ' //' U2 Rel. Err. ' //' P  Rel. Err. ' &
      //' S1 Rel. Err. ' //' S2 Rel. Err. ' //' S3 Rel. Err. ' &
      //' Defect  Err. ')
      call output_line ('--------------------------------------'//&
          '---------------------------------------------------------'//&
          '---------')
    call output_line (sys_siL(inl, 5) //'  '&
    //trim(sys_sdEL(Dres_rel(1),6))//'  '&
    //trim(sys_sdEL(Dres_rel(2),6))//'  '&
    //trim(sys_sdEL(Dres_rel(3),6))//'  '&
    //trim(sys_sdEL(Dres_rel(4),6))//'  '&
    //trim(sys_sdEL(Dres_rel(5),6))//'  '&
    //trim(sys_sdEL(Dres_rel(6),6))//'  '&
    //trim(sys_sdEL(defnorm(1),6)))   
  else
    call output_line (sys_siL(inl, 5) //'  '&
    //trim(sys_sdEL(Dres_rel(1),6))//'  '&
    //trim(sys_sdEL(Dres_rel(2),6))//'  '&
    //trim(sys_sdEL(Dres_rel(3),6))//'  '&
    //trim(sys_sdEL(Dres_rel(4),6))//'  '&
    //trim(sys_sdEL(Dres_rel(5),6))//'  '&
    //trim(sys_sdEL(Dres_rel(6),6))//'  '&
    //trim(sys_sdEL(defnorm(1),6))) 
    if ( (mod(inl,10) .eq. 0) .and. (inl .ne. NLN_Max) &
      .and. (.not. converged) .and. (.not. diverged)) then
      call output_lbrk()
      call output_line ('Iter. ' &
      //' U1 Rel. Err. ' //' U2 Rel. Err. ' //' P  Rel. Err. ' &
      //' S1 Rel. Err. ' //' S2 Rel. Err. ' //' S3 Rel. Err. ' &
      //' Defect  Err. ')
      call output_line ('--------------------------------------'//&
          '---------------------------------------------------------'//&
          '---------')
    end if
  end if

  ! Release the difference vector
  call lsysbl_releaseVector (rdiff)     
  ! Release the defect vector
  call lsysbl_releaseVector (rdefect)
  
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
  
  ! Solution vectors in the previous nonliner iteration
  ! and the RHS vector
  type(t_vectorBlock), intent(inout) :: rvector_old,rrhs
 !<\inputoutput>
 
  !<input>
  ! Nonlinear loop's maximum/current number of iterations
  integer, intent(in) :: NLMAX,NLMIN
  
  ! Solution vectors in the current nonliner iteration
  type(t_vectorBlock), intent(in) :: rvector
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
    !  but since the loop is converged, continue
    !  with postprocessing. 
    det = .true.
    
  end if
 
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_postprocess(rboundary,rmatrix,rvector,rtriangulation,&
              rcubatureInfo,rdiscretisation,rparams,rcollection)
                
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
  
  ! Collection structure for callback routines  
  type(t_collection), intent(inout) :: rcollection  
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
  integer :: KEnergy, detKEnergy, Div
  integer :: Vtild, Ptild, Stild
  
  ! Kinematic viscosity noo = 1/Re  
  real(DP) :: dnu
     
   ! Path to the data file which has the initial solution
  character(LEN=SYS_STRLEN) :: sfile
 
  ! Output block for UCD output to GMV/VTK file
  type(t_ucdExport) :: rexport
  character(len=SYS_STRLEN) :: sucddir
  real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
  integer :: NLMAX, NLMIN
  
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


  ! Norm Calculations
  integer :: L2U,L2P,L2S,H1U,H1P,StbblError
  real(DP) :: dintvalue, dintvalue1
  type(t_fev2Vectors) :: revalVectors
  ! Cubature information structure for static bubble
  type(t_scalarCubatureInfo) :: rcubatureInfo2

  ! We have solved our NS problem on level...
  call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 5)
  call parlst_getvalue_int (rparams, 'MESH', 'NLMIN', NLMIN, 3)
   
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
  ! Calculate drag-/lift coefficients on the 2nd boundary component.
  ! This is for the benchmark problem: flow around cylinder!
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

  ! Determine whether to calculate flow around cylinder parameters or not
  call parlst_getvalue_int (rparams, 'POST', 'LiftDragASO', LiftDragASO, 0)
  
  if (LiftDragASO .eq. 1) then
    
    call boundary_createRegion (rboundary,2,0, rboundaryRegion)
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
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

    
    ! Calculate the Lift and Drag coefficients using the
    !  stress tensor calculated.
    ! Print out the negatives of the values computed in the
    ! following routines as the routines work with the normal
    ! vector pointing into the opposite direction as needed.
    call bdint_normalFlux2D (rvector%RvectorBlock(4),&
        rvector%RvectorBlock(5),CUB_G3_1D,Dforces(1), rboundaryRegion)
    call bdint_normalFlux2D (rvector%RvectorBlock(5),&
        rvector%RvectorBlock(6),CUB_G3_1D,Dforces(2), rboundaryRegion)
        
    call output_lbrk()
    call output_line ('Coefficients (Direct calculation)')
    call output_line ('--------------------------------')
    call output_line ('Drag/Lift')
    call output_line (trim(sys_sdEP(-Dforces(1)*500.0_DP,15,6)) // ' / '&
      //trim(sys_sdEP(-Dforces(2)*500.0_DP,15,6)))        
    call output_lbrk()
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
  call parlst_getvalue_int (rparams, 'MESH', 'Stild', Stild, 0)

  if ( (Vtild .eq. 1) .or. (Ptild .eq. 1) .or. (Stild .eq. 1)) then
      
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

    if (Stild .eq. 1) then
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(4), &
      EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(4))         
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(5), &
      EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(5))
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(6), &
      EL_Q1, CUB_G3_2D, rprjDiscretisation%RspatialDiscr(6))
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
    call ls_BCs_Dirichlet_One(rprjDiscretisation,rboundary,rprjDiscreteBC,&
                              rcollection)

    ! Hang the pointer into the vector.
    rprjVector%p_rdiscreteBC => rprjDiscreteBC

    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector) 
    
    if (ExporType .eq. 0) then
    
      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nslsfem-'//trim(sys_si0L(NLMAX,1))//&
        '-'//trim(sys_si0L(NLMIN,1))//'.vtk')

      ! Write Pressure
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

      ! Write velocity field
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

      ! Write Stresses
      call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sx',UCD_VAR_STANDARD,p_Ddata)

      call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sxy',UCD_VAR_STANDARD,p_Ddata)

      call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sy',UCD_VAR_STANDARD,p_Ddata)      
    
    else
      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nslsfem-'//trim(sys_si0L(NLMAX,1))//&
        '-'//trim(sys_si0L(NLMIN,1))//'.gmv')

      ! Write Pressure
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
      
      ! Write velocity field     
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
      
      ! Write Stresses
      call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sx',UCD_VAR_STANDARD,p_Ddata)     

      call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sxy',UCD_VAR_STANDARD,p_Ddata) 
      
      call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sy',UCD_VAR_STANDARD,p_Ddata)       
    
    end if
    
    ! Release the temporary projected vector
    call lsysbl_releaseVector (rprjVector)

    ! Release our discrete version of the projected boundary conditions
    call bcasm_releaseDiscreteBC (rprjDiscreteBC)
    ! Release the projected discretisation structure and 
    ! all spatial discretisation structures in it.
    call spdiscr_releaseBlockDiscr(rprjDiscretisation)
  
  else ! real data will be used
  
    if (ExporType .eq. 0) then
    
      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nslsfem-'//trim(sys_si0L(NLMAX,1))//&
        '-'//trim(sys_si0L(NLMIN,1))//'.vtk')

      ! Write Pressure
      call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

      ! Write velocity field
      call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

      ! Write Stresses
      call lsyssc_getbase_double (rvector%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sx',UCD_VAR_STANDARD,p_Ddata)
      
      call lsyssc_getbase_double (rvector%RvectorBlock(5),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sxy',UCD_VAR_STANDARD,p_Ddata)
      
      call lsyssc_getbase_double (rvector%RvectorBlock(6),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sy',UCD_VAR_STANDARD,p_Ddata)            
    
    else
      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/nslsfem-'//trim(sys_si0L(NLMAX,1))//&
        '-'//trim(sys_si0L(NLMIN,1))//'.gmv')

      ! Write Pressure
      call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
      
      ! Write velocity field     
      call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)
      
      ! Write Stresses
      call lsyssc_getbase_double (rvector%RvectorBlock(4),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sx',UCD_VAR_STANDARD,p_Ddata)     

      call lsyssc_getbase_double (rvector%RvectorBlock(5),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sxy',UCD_VAR_STANDARD,p_Ddata)
    
      call lsyssc_getbase_double (rvector%RvectorBlock(6),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'Sy',UCD_VAR_STANDARD,p_Ddata)    
    
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
  ! Calculating the velocity divergence L^2-error.
  !   Div = \int{ (u_x + v_y)^2 dX}
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  ! Determine whether to calculate the Divergence L^2-Error
  call parlst_getvalue_int (rparams, 'POST', 'Div', Div, 0)  
    
  if (Div .eq. 1) then
    ! Add velocity vectors
    rcollection%IquickAccess(7) = 7
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(1),1)   ! u1,x,y
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(2),1)   ! u2,x,y       
       
    ! initializing the maximum norm value
    rcollection%DquickAccess(11) = 1.0E-16_DP
    
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error Div(u)')
    call output_line ('----------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('L2divu:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG) 
 
    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^inf Error Div(u)')
    call output_line ('------------------')
    call output_line (trim(sys_sdEP(rcollection%DquickAccess(11),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('Linfdivu:'//&
    trim(sys_sdEP(rcollection%DquickAccess(11),15,6)), & 
    coutputMode=OU_MODE_BENCHLOG) 
 
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
      call output_line ('flux input: -0.16666667')
      call output_line ('-----------------------')
      call output_line (trim(sys_sdEP(Dfluxi,17,10)))
      
      ! Better to use the exact value of the inflow fluxes, rather than
      !  calculating it numerically
      ! Flow Around Cylinder
      ! Dfluxi = -0.082_DP
      ! Poiseuelle Flow
      Dfluxi = -1.0_DP/6.0_DP
!      Dfluxi = -0.082_DP
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
      call output_line (trim(sys_sdEP(Dgmc,16,6)))     
      

      ! Output line flux
      Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
      Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
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
      call output_line (trim(sys_sdEP(Dgmc,16,6)))    
    
    else
      ! Calculate the number of elements from nlevels
      ! Input line coordinates
      call parlst_getvalue_string (rparams, 'POST', 'inputGMC',sstring)
      read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
                 Dcoords(1,2), Dcoords(2,2)
      ! Input line flux
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                     Dcoords(1:2,2),Dfluxi,nlevels=nlevels)
      
      ! Output line coordinates
      call parlst_getvalue_string (rparams, 'POST', 'outputGMC',sstring)
      read (sstring,*)   Dcoords(1,1), Dcoords(2,1), &
                 Dcoords(1,2), Dcoords(2,2)
      ! Output line flux
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                    Dcoords(1:2,2),Dfluxo,nlevels=nlevels)

      ! Print out the calculated inflow flux
      call output_lbrk()
!      call output_line ('flux input: -0.1666667')
!      call output_line ('----------------------')      
      call output_line ('flux input: -0.082')
      call output_line ('------------------')
      call output_line (trim(sys_sdEP(Dfluxi,17,10)))
      
      ! Exact value of the inflow flux
      ! ***Flow Around Cylinder*** !
      Dfluxi = -0.082_DP
      ! ***Poiseuille Flow*** !
      ! Dfluxi = -1.0_DP/6.0_DP
      
      ! Take the absolute flux values instead
      Dfluxi = abs(Dfluxi)
      Dfluxo = abs(Dfluxo)
            
      ! The GMC is then calculated as the normalized flux value
      Dgmc = 100.0_DP*(Dfluxi - Dfluxo) / Dfluxi
      
      ! Print the GMC value
      call output_lbrk()
      call output_line ('Global Mass Conservation(%)')
      call output_line (&
      '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
      call output_line (trim(sys_sdEP(Dgmc,16,6)))     


      ! Modify the coordinates for a new cross-section
      ! Output line flux will be calculated on:
      Dcoords(1,1) = Dcoords(1,1) + 0.45_DP
      Dcoords(1,2) = Dcoords(1,2) + 0.45_DP
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
      
      ! Take the absolute flux value instead
      Dfluxo5 = abs(Dfluxo5)
      
      ! The GMC is then calculated as
      Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                           
      
      ! Print the GMC value
      call output_lbrk()
      call output_line ('Global Mass Conservation(%)')
      call output_line (&
      '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
      call output_line (trim(sys_sdEP(Dgmc,16,6)))

      ! Modify the coordinates for a new cross-section
      ! Output line flux will be calculated on:
      Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
      Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
      
      ! Take the absolute flux value instead
      Dfluxo5 = abs(Dfluxo5)
      
      ! The GMC is then calculated as
      Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                           
      
      ! Print the GMC value
      call output_lbrk()
      call output_line ('Global Mass Conservation(%)')
      call output_line (&
      '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
      call output_line (trim(sys_sdEP(Dgmc,16,6)))


      ! Modify the coordinates for a new cross-section
      ! Output line flux will be calculated on:
      Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
      Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
      
      ! Take the absolute flux value instead
      Dfluxo5 = abs(Dfluxo5)
      
      ! The GMC is then calculated as
      Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                           
      
      ! Print the GMC value
      call output_lbrk()
      call output_line ('Global Mass Conservation(%)')
      call output_line (&
      '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
      call output_line (trim(sys_sdEP(Dgmc,16,6)))
      
      ! Modify the coordinates for a new cross-section
      ! Output line flux will be calculated on:
      Dcoords(1,1) = Dcoords(1,1) + 0.5_DP
      Dcoords(1,2) = Dcoords(1,2) + 0.5_DP
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
      
      ! Take the absolute flux value instead
      Dfluxo5 = abs(Dfluxo5)
      
      ! The GMC is then calculated as
      Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                           
      
      ! Print the GMC value
      call output_lbrk()
      call output_line ('Global Mass Conservation(%)')
      call output_line (&
      '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
      call output_line (trim(sys_sdEP(Dgmc,16,6)))
      
            ! Modify the coordinates for a new cross-section
      ! Output line flux will be calculated on:
      Dcoords(1,1) = Dcoords(1,1) + 0.2_DP
      Dcoords(1,2) = Dcoords(1,2) + 0.2_DP
      call ppns2D_calcFluxThroughLine (rvector,Dcoords(1:2,1),&
                     Dcoords(1:2,2),Dfluxo5,nlevels=nlevels)
      
      ! Take the absolute flux value instead
      Dfluxo5 = abs(Dfluxo5)
      
      ! The GMC is then calculated as
      Dgmc = 100.0_DP*(Dfluxi - Dfluxo5) / Dfluxi                           
      
      ! Print the GMC value
      call output_lbrk()
      call output_line ('Global Mass Conservation(%)')
      call output_line (&
      '--------at x='//trim(sys_sdp(Dcoords(1,1),5,2))//'-------------')
      call output_line (trim(sys_sdEP(Dgmc,16,6)))      
     
   end if   
  
  end if


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate The L^2 and H^1 norms for analytical solutions
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call parlst_getvalue_int (rparams, 'POST', 'L2U', L2U, 0)
  call parlst_getvalue_int (rparams, 'POST', 'L2P', L2P, 0)
  call parlst_getvalue_int (rparams, 'POST', 'L2S', L2S, 0)
  call parlst_getvalue_int (rparams, 'POST', 'H1U', H1U, 0)
  call parlst_getvalue_int (rparams, 'POST', 'H1P', H1P, 0)
  call parlst_getvalue_int (rparams, 'POST', 'StbblError', StbblError, 0)
  
  if (L2U == 1) then

    ! L^2 Norm velocity
    ! Add x-velocity vector
    rcollection%IquickAccess(7) = 0
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(1),0)   ! u1
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)    
    
    ! Add y-velocity vector
    rcollection%IquickAccess(7) = 1      
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(2),0)   ! u2
    call bma_buildIntegral (dintvalue1,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)
    
    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error velocity')
    call output_line ('--------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)))  
    
    call output_line ('L2velocity:'//&
    trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)), coutputMode=OU_MODE_BENCHLOG)

  end if

  if (H1U == 1) then    
    ! H^1 Norm velocity
    ! Add x-velocity vector
    rcollection%IquickAccess(7) = 0
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(1),1)   ! u1,x,y
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_H1_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)    
    
    ! Add y-velocity vector
    rcollection%IquickAccess(7) = 1      
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(2),1)   ! u2,x,y
    call bma_buildIntegral (dintvalue1,BMA_CALC_STANDARD,&
    ls_H1_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)
    
    ! Print the Norm value
    call output_lbrk()
    call output_line ('H^1 Error velocity')
    call output_line ('------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)))  

    call output_line ('H1velocity:'//&
    trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)), coutputMode=OU_MODE_BENCHLOG)

  end if
  
  if (L2P == 1) then  
    ! L^2 Norm pressure
    ! Add pressure vector
    rcollection%IquickAccess(7) = 2
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(3),0)   ! p
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error pressure')
    call output_line ('------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('L2pressure:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG)

  end if

  if (H1P == 1) then
    ! H^1 Norm pressure
    ! Add pressure vector
    rcollection%IquickAccess(7) = 2
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(3),1)   ! p,x,y
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_H1_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)    
    
    ! Print the Norm value
    call output_lbrk()
    call output_line ('H^1 Error pressure')
    call output_line ('------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('H1pressure:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG)    

  end if

    
  if (L2S == 1) then    
    ! L^2 Norm Sxx
    ! Add the vector
    rcollection%IquickAccess(7) = 3
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(4),0)   ! sxx
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error Sxx')
    call output_line ('-------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('L2Sxx:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG) 


    ! L^2 Norm Sxy
    ! Add the vector
    rcollection%IquickAccess(7) = 4
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(5),0)   ! sxy
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error Sxy')
    call output_line ('-------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('L2Sxy:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG)     
  
    ! L^2 Norm Syy
    ! Add the vector
    rcollection%IquickAccess(7) = 5
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(6),0)   ! syy
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error Syy')
    call output_line ('-------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)  

    call output_line ('L2Syy:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG) 

  end if
     
  if (StbblError == 1) then

    !!!! Calculate the errors 
    ! Add the pressure vector
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(3),0)   ! p
       
    ! Gauss 1-pt rule = 1 Point per element in the center.
    call spdiscr_createDefCubStructure(&  
    rdiscretisation%RspatialDiscr(3),rcubatureInfo2,CUB_G1_2D)

    ! Initialize the minimum pressure
    rcollection%DquickAccess(16) = 1.0E12_DP
    ! Calculate the real minimum of the pressure field
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_calc_min,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo2)
    
    ! Tell that we want to calculate the static bubble errors
    rcollection%IquickAccess(7) = 6
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo2)

    !!!! Print the errors to standard output and the log file
    call output_lbrk()
    call output_line ('Mean pressure inside bubble')
    call output_line ('---------------------------')
    call output_line (trim(sys_sdEP(&
              rcollection%DquickAccess(13)/0.1963495408493621_DP,15,6)))  
!    call output_lbrk()
!    call output_line ('Numerical area inside bubble')
!    call output_line ('----------------------------')
!    call output_line (trim(sys_sdEP(rcollection%DquickAccess(12),15,6)))

    call output_lbrk()
    call output_line ('Mean pressure outside bubble')
    call output_line ('----------------------------')
    call output_line (trim(sys_sdEP(&
              rcollection%DquickAccess(15)/3.8036504591506379_DP,15,6)))  
!    call output_lbrk()
!    call output_line ('Numerical area outside bubble')
!    call output_line ('-----------------------------')
!    call output_line (trim(sys_sdEP(rcollection%DquickAccess(14),15,6)))

    call output_lbrk()
    call output_line ('Abs. Error Young-Laplace Eq.')
    call output_line ('----------------------------')
    call output_line ( trim( sys_sdEP(abs(&
              rcollection%DquickAccess(13)/0.1963495408493621_DP-&
              rcollection%DquickAccess(15)/3.8036504591506379_DP-&
              4.0_DP),15,6)))
    call output_lbrk()
    call output_line ('Percent Error(%) Young-Laplace Eq.')
    call output_line ('----------------------------------')
    call output_line (trim(sys_sdEP(abs(&
              (rcollection%DquickAccess(13)/0.1963495408493621_DP-&
              rcollection%DquickAccess(15)/3.8036504591506379_DP-&
              4.0_DP)*25.0_DP),15,6)))
    
    !!!! Release the temp vector and the cubature structure
    call fev2_releaseVectorList(revalVectors)
    call spdiscr_releaseCubStructure(rcubatureInfo2)
    
  end if           
     
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_cleanup(rvector,rvector_old,rrhs,Rlevels,p_rtempVectorSc,&
                        rboundary,NLMAX,NLMIN,rparams)
 
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
  type(t_vectorBlock) :: rvector,rvector_old,rrhs
  
  ! Temporary scalar vector; used for calculating the nonlinear matrix
  ! on lower levels / projecting the solution from higher to lower levels.
  type(t_vectorScalar), pointer :: p_rtempVectorSc
  !</inputoutput>

!</subroutine>

  ! Local variables
  integer :: i
  
  ! Determine if zero mean pressure constraint was active 
  integer :: detZMV
  
  ! Zero mean pressure active or not
  call parlst_getvalue_int (rparams, 'ZMV', 'detZMV', detZMV, 0)
  
  ! Now, clean up so that all the memory is available again.  
  ! Release the block matrix/vectors
  call lsysbl_releaseVector (rvector)
  call lsysbl_releaseVector (rvector_old)
  call lsysbl_releaseVector (rrhs)
  
  do i = NLMAX, NLMIN, -1
    call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
    if (detZMV .eq. 1) call lsysbl_releaseMatrix (Rlevels(i)%rmatrixMass)
    if (i .ne. NLMAX) then
      call lsysbl_releaseVector (Rlevels(i)%rtempVector)
    end if
  end do
  if (detZMV .eq. 4) call lsysbl_releaseMatrix (Rlevels(NLMIN)%rmatrixMass)
  
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
  subroutine ls_svp2D_Matrix(RmatrixData,rassemblyData,rmatrixAssembly,&
             npointsPerElement,nelements,revalVectors,rcollection)
!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
  ! This is where the (BIG-BANG)**2 starts :D
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
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA12
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA13,p_DlocalMatrixA14
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA15,p_DlocalMatrixA16
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA21,p_DlocalMatrixA22
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA23,p_DlocalMatrixA24
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA25,p_DlocalMatrixA26  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA33,p_DlocalMatrixA34  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA36   
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA41
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA43,p_DlocalMatrixA44  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA45  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA51,p_DlocalMatrixA52
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA54  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA55,p_DlocalMatrixA56  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA62
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA63 
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA65,p_DlocalMatrixA66  
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA33,p_DbasTestA33
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA44,p_DbasTestA44
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA14,p_DbasTestA14
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA34,p_DbasTestA34

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11,p_rmatrixDataA33,p_rmatrixDataA44
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA13,p_rmatrixDataA14,p_rmatrixDataA34
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

  ! Scaling factors
  real(DP) :: beta, gama
  
  ! Continuity equation scaling
  real(DP) :: alpha
  
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dUx, dUy, dVx, dVy
  
  ! Viscosity and its 1st derivative
  ! 2nd invariant and a frequently used value
  real(DP) :: nu, dnu, Dii, F
  
  ! Continuity equation scaling factor ---> \alpha
  !  default = 1.0_DP  no scaling
  alpha = rcollection%DquickAccess(3)

  ! Linearization Scheme ---> \beta
  beta = rcollection%DquickAccess(4)

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
  p_DlocalMatrixA15 => RmatrixData(1,5)%p_Dentry
  p_DlocalMatrixA16 => RmatrixData(1,6)%p_Dentry

  p_DlocalMatrixA21 => RmatrixData(2,1)%p_Dentry
  p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
  p_DlocalMatrixA23 => RmatrixData(2,3)%p_Dentry
  p_DlocalMatrixA24 => RmatrixData(2,4)%p_Dentry
  p_DlocalMatrixA25 => RmatrixData(2,5)%p_Dentry
  p_DlocalMatrixA26 => RmatrixData(2,6)%p_Dentry
  
  p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
  p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
  p_DlocalMatrixA33 => RmatrixData(3,3)%p_Dentry
  p_DlocalMatrixA34 => RmatrixData(3,4)%p_Dentry
  p_DlocalMatrixA36 => RmatrixData(3,6)%p_Dentry   

  p_DlocalMatrixA41 => RmatrixData(4,1)%p_Dentry
  p_DlocalMatrixA43 => RmatrixData(4,3)%p_Dentry
  p_DlocalMatrixA44 => RmatrixData(4,4)%p_Dentry
  p_DlocalMatrixA45 => RmatrixData(4,5)%p_Dentry   
  
  p_DlocalMatrixA51 => RmatrixData(5,1)%p_Dentry
  p_DlocalMatrixA52 => RmatrixData(5,2)%p_Dentry
  p_DlocalMatrixA54 => RmatrixData(5,4)%p_Dentry
  p_DlocalMatrixA55 => RmatrixData(5,5)%p_Dentry
  p_DlocalMatrixA56 => RmatrixData(5,6)%p_Dentry   
  
  p_DlocalMatrixA62 => RmatrixData(6,2)%p_Dentry
  p_DlocalMatrixA63 => RmatrixData(6,3)%p_Dentry
  p_DlocalMatrixA65 => RmatrixData(6,5)%p_Dentry
  p_DlocalMatrixA66 => RmatrixData(6,6)%p_Dentry  
  
  ! Get the velocity field from the parameters
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
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
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    call ls_nonlinear_weight(Dii,rcollection,gama)
    ! And a frequently used combination
    F = dUy+dVx
    
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
      dval = p_DcubWeight(icubp,iel)*beta*(   dU*dUx*dbasJx*dbasI+&
             dV*dUx*dbasJy*dbasI + 4.0_DP*gama*&
             nu*dnu*(4.0_DP*dUx**2*dbasJx*dbasIx+ &
             2.0_DP*dUx*F*dbasJx*dbasIy+2.0_DP*F*dUx*dbasJy*dbasIx+ &
             F**2*dbasJy*dbasIy)   )

      dval = dval + p_DcubWeight(icubp,iel)*(   dU**2*dbasJx*dbasIx+dU*dV*dbasJx*dbasIy + &
             dV*dU*dbasJy*dbasIx + dV**2*dbasJy*dbasIy + &
             gama*nu**2*(4.0_DP*dbasJx*dbasIx + 2.0_DP*dbasJy*dbasIy) + & 
             alpha*(dbasJx*dbasIx)  )           
             
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
      ! A22
      dval = p_DcubWeight(icubp,iel)*beta*(   dU*dVy*dbasJx*dbasI+&
             dV*dVy*dbasJy*dbasI + 4.0_DP*gama* &
             nu*dnu*(4.0_DP*dVy**2*dbasJy*dbasIy+ &
             2.0_DP*dVy*F*dbasJy*dbasIx+2.0_DP*F*dVy*dbasJx*dbasIy+ &
             F**2*dbasJx*dbasIx)   )
      
      dval = dval + p_DcubWeight(icubp,iel)*(   dU**2*dbasJx*dbasIx+dU*dV*dbasJy*dbasIx + &
             dV*dU*dbasJx*dbasIy + dV**2*dbasJy*dbasIy + &
             gama*nu**2*(4.0_DP*dbasJy*dbasIy + 2.0_DP*dbasJx*dbasIx) + &
             alpha*(dbasJy*dbasIy)  )         
              
      p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
      
      ! A12
      dval = p_DcubWeight(icubp,iel)*beta*(   dU*dVx*dbasJx*dbasI+&
             dV*dVx*dbasJy*dbasI + & 
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy)   ) 
             
      dval = dval + p_DcubWeight(icubp,iel)*(   gama*nu**2*(2.0_DP*dbasJx*dbasIy)+& 
             alpha*(dbasJy*dbasIx)  ) 
            
      p_DlocalMatrixA12(jdofe,idofe,iel) = p_DlocalMatrixA12(jdofe,idofe,iel) + dval
      
      ! A21
      dval = p_DcubWeight(icubp,iel)*beta*(   dUy*dU*dbasJ*dbasIx + & 
             dUy*dV*dbasJ*dbasIy + & 
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy)   ) 
             
      dval = dval + p_DcubWeight(icubp,iel)*(   gama*nu**2*(2.0_DP*dbasJx*dbasIy)+& 
             alpha*(dbasJy*dbasIx)  )
                   
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
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    call ls_nonlinear_weight(Dii,rcollection,gama)    
    ! And a frequently used combination
    F = dUy+dVx
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA13%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasIx = p_DbasTestA13(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA13(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA13%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA13(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A13
      dval = -p_DcubWeight(icubp,iel)*gama*beta*4.0_DP*dnu* &
             (   (2.0_DP*dUx**2+2.0_DP*dVy*dUx)*dbasJ*dbasIx + &
             (dUx*F+dVy*F)*dbasJ*dbasIy    )
             
      dval = dval - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIx   )
      
      p_DlocalMatrixA13(jdofe,idofe,iel) = p_DlocalMatrixA13(jdofe,idofe,iel) + dval
      
      ! A31
      dval = - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIx   )      
      p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
           
      ! A23
      dval = -p_DcubWeight(icubp,iel)*gama*beta*4.0_DP*dnu* &
             (   (2.0_DP*dVy**2+2.0_DP*dVy*dUx)*dbasJ*dbasIy + &
             (dUx*F+dVy*F)*dbasJ*dbasIx   )

      dval = dval - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIy   )           
             
      p_DlocalMatrixA23(jdofe,idofe,iel) = p_DlocalMatrixA23(jdofe,idofe,iel) + dval
      
      ! A32
      dval = - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIy   )       
      p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + dval
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  
  ! ++++++++++++++++++++++++++++++++++++++
  ! Calculate blocks:   A14, A24, A41, A42
  ! A15, A25, A51, A52, A16, A26, A61, A62
  ! ++++++++++++++++++++++++++++++++++++++
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

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    call ls_nonlinear_weight(Dii,rcollection,gama)    
    ! And a frequently used combination
    F = dUy+dVx
    
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
      dval = -p_DcubWeight(icubp,iel) * beta* (   dUx*dbasJx*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dUx**2*dbasJ*dbasIx+dUx*F*dbasJ*dbasIy)  )
             
      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )             
             
      p_DlocalMatrixA14(jdofe,idofe,iel) = p_DlocalMatrixA14(jdofe,idofe,iel) + dval
      
      ! A41             
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )
                   
      p_DlocalMatrixA41(idofe,jdofe,iel) = p_DlocalMatrixA41(idofe,jdofe,iel) + dval

      ! A24
      dval = -p_DcubWeight(icubp,iel) * beta*(   dUy*dbasJx*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dUx*dVy*dbasJ*dbasIy+dUx*F*dbasJ*dbasIx)*gama   )
      p_DlocalMatrixA24(jdofe,idofe,iel) = p_DlocalMatrixA24(jdofe,idofe,iel) + dval
      
      ! A42

      ! A15
      dval = -p_DcubWeight(icubp,iel)*beta*(   dUx*dbasJy*dbasI+dVx*dbasJx*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dUx*F*dbasJ*dbasIx+F**2*dbasJ*dbasIy)   )
      
      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )
             
      p_DlocalMatrixA15(jdofe,idofe,iel) = p_DlocalMatrixA15(jdofe,idofe,iel) + dval
      
      ! A51
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )
      p_DlocalMatrixA51(idofe,jdofe,iel) = p_DlocalMatrixA51(idofe,jdofe,iel) + dval
                      
      ! A25
      dval = -p_DcubWeight(icubp,iel)*beta*(   dUy*dbasJy*dbasI+dVy*dbasJx*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dVy*F*dbasJ*dbasIy+F**2*dbasJ*dbasIx)   )

      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )           
             
      p_DlocalMatrixA25(jdofe,idofe,iel) = p_DlocalMatrixA25(jdofe,idofe,iel) + dval
      
      ! A52
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )  
      p_DlocalMatrixA52(idofe,jdofe,iel) = p_DlocalMatrixA52(idofe,jdofe,iel) + dval
      
      ! A16
      dval = -p_DcubWeight(icubp,iel) * beta*(   dVx*dbasJy*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dVy*dUx*dbasJ*dbasIx+dVy*F*dbasJ*dbasIy)*gama   )
      p_DlocalMatrixA16(jdofe,idofe,iel) = p_DlocalMatrixA16(jdofe,idofe,iel) + dval
      
      ! A61   
      
      ! A26
      dval = -p_DcubWeight(icubp,iel)*beta*(   dVy*dbasJy*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dVy**2*dbasJ*dbasIy+dVy*F*dbasJ*dbasIx)   )
             
      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )            
             
      p_DlocalMatrixA26(jdofe,idofe,iel) = p_DlocalMatrixA26(jdofe,idofe,iel) + dval
      
      ! A62  
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )       
      p_DlocalMatrixA62(idofe,jdofe,iel) = p_DlocalMatrixA62(idofe,jdofe,iel) + dval      
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

    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_nonlinear_weight(Dii,rcollection,gama) 

    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA33%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA33(idofe,DER_FUNC,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA33%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA33(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A33
      dval = p_DcubWeight(icubp,iel) * 2.0_DP*gama*dbasJ*dbasI
      p_DlocalMatrixA33(jdofe,idofe,iel) = p_DlocalMatrixA33(jdofe,idofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  

  ! +++++++++++++++++++++++++
  ! Calculate blocks A34, A43
  ! A35, A53, A36, A63 
  ! +++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_nonlinear_weight(Dii,rcollection,gama) 
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA34%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA34(idofe,DER_FUNC,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA34%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA34(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A34
      dval = p_DcubWeight(icubp,iel) * gama*dbasJ*dbasI
      p_DlocalMatrixA34(jdofe,idofe,iel) = p_DlocalMatrixA34(jdofe,idofe,iel) + dval

      ! A43
      p_DlocalMatrixA43(idofe,jdofe,iel) = p_DlocalMatrixA43(idofe,jdofe,iel) + dval

      ! A35
      ! A53

      ! A36
      p_DlocalMatrixA36(jdofe,idofe,iel) = p_DlocalMatrixA36(jdofe,idofe,iel) + dval

      ! A63
      p_DlocalMatrixA63(idofe,jdofe,iel) = p_DlocalMatrixA63(idofe,jdofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  
  
  ! +++++++++++++++++++++++++++++
  ! Calculate block A44, A55, A66
  !  A45, A54, A46, A64, A56, A65
  ! +++++++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_nonlinear_weight(Dii,rcollection,gama) 

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
      dval = p_DcubWeight(icubp,iel) * (  dbasJx*dbasIx + gama*dbasJ*dbasI  )
      p_DlocalMatrixA44(jdofe,idofe,iel) = p_DlocalMatrixA44(jdofe,idofe,iel) + dval

      ! A55
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIy + dbasJx*dbasIx + &
            2.0_DP*gama*dbasJ*dbasI  )
      p_DlocalMatrixA55(jdofe,idofe,iel) = p_DlocalMatrixA55(jdofe,idofe,iel) + dval

      ! A66
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIy + gama*dbasJ*dbasI  )
      p_DlocalMatrixA66(jdofe,idofe,iel) = p_DlocalMatrixA66(jdofe,idofe,iel) + dval
      
      ! A45
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIx  )
      p_DlocalMatrixA45(jdofe,idofe,iel) = p_DlocalMatrixA45(jdofe,idofe,iel) + dval  

      ! A54
      p_DlocalMatrixA54(idofe,jdofe,iel) = p_DlocalMatrixA54(idofe,jdofe,iel) + dval

      ! A46
      ! A64

      ! A56
      p_DlocalMatrixA56(jdofe,idofe,iel) = p_DlocalMatrixA56(jdofe,idofe,iel) + dval  

      ! A65
      p_DlocalMatrixA65(idofe,jdofe,iel) = p_DlocalMatrixA65(idofe,jdofe,iel) + dval

      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  end subroutine


 
  ! ***************************************************************************

!<subroutine>
  subroutine ls_svp2D_NewtonMatrix(RmatrixData,rassemblyData,rmatrixAssembly,&
             npointsPerElement,nelements,revalVectors,rcollection)
!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
  ! This is where the (BIG-BANG)**2 starts :D
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
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA12
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA21,p_DlocalMatrixA22
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA41,p_DlocalMatrixA42
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA51,p_DlocalMatrixA52
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA61,p_DlocalMatrixA62
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA14,p_DbasTestA14

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA13,p_rmatrixDataA14
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

  ! Scaling factors
  real(DP) :: gama
  
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dUx, dUy, dVx, dVy
  
  ! Viscosity and its 1st derivative
  ! 2nd invariant and a frequently used value
  real(DP) :: nu, dnu, Dii, F
  

  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
  p_rmatrixDataA13 => RmatrixData(1,3)
  p_rmatrixDataA14 => RmatrixData(1,4)


  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
  p_DbasTrialA13 => RmatrixData(1,3)%p_DbasTrial
  p_DbasTestA13 => RmatrixData(1,3)%p_DbasTest
  p_DbasTrialA14 => RmatrixData(1,4)%p_DbasTrial
  p_DbasTestA14 => RmatrixData(1,4)%p_DbasTest
  
  
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
  p_DlocalMatrixA12 => RmatrixData(1,2)%p_Dentry

  p_DlocalMatrixA21 => RmatrixData(2,1)%p_Dentry
  p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
  
  p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
  p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
 
  p_DlocalMatrixA41 => RmatrixData(4,1)%p_Dentry
  p_DlocalMatrixA42 => RmatrixData(4,2)%p_Dentry
  
  p_DlocalMatrixA51 => RmatrixData(5,1)%p_Dentry
  p_DlocalMatrixA52 => RmatrixData(5,2)%p_Dentry
  
  p_DlocalMatrixA61 => RmatrixData(6,1)%p_Dentry
  p_DlocalMatrixA62 => RmatrixData(6,2)%p_Dentry

  
  ! Get the velocity field from the parameters
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
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
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    call ls_nonlinear_weight(Dii,rcollection,gama)
    ! And a frequently used combination
    F = dUy+dVx
    
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
      dval = p_DcubWeight(icubp,iel)*(  dUx*dU*dbasJ*dbasIx+&
             dUx*dV*dbasJ*dbasIy + (dUx**2+dVx**2)*dbasJ*dbasI +&
             4.0_DP*gama*nu*dnu*(4.0_DP*dUx**2*dbasJx*dbasIx+ &
             2.0_DP*dUx*F*dbasJx*dbasIy+2.0_DP*F*dUx*dbasJy*dbasIx+ &
             F**2*dbasJy*dbasIy) + Dii*32.0_DP*dnu**2*(dUx**2*dbasJx*dbasIx+ &
             0.5_DP*dUx*F*dbasJx*dbasIy+0.5_DP*F*dUx*dbasJy*dbasIx+ &
             0.25_DP*F**2*dbasJy*dbasIy)*gama    )         
             
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
      ! A22
      dval = p_DcubWeight(icubp,iel)*(  dVy*dU*dbasJ*dbasIx+&
             dVy*dV*dbasJ*dbasIy + (dUy**2+dVy**2)*dbasJ*dbasI + &
             4.0_DP*gama*nu*dnu*(4.0_DP*dVy**2*dbasJy*dbasIy+ &
             2.0_DP*dVy*F*dbasJy*dbasIx+2.0_DP*F*dVy*dbasJx*dbasIy+ &
             F**2*dbasJx*dbasIx) + Dii*32.0_DP*dnu**2*(dVy**2*dbasJy*dbasIy+ &
             0.5_DP*dVy*F*dbasJy*dbasIx+0.5_DP*F*dVy*dbasJx*dbasIy+ &
             0.25_DP*F**2*dbasJx*dbasIx)*gama   )      
              
      p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
      
      ! A12
      dval = p_DcubWeight(icubp,iel)*(  dUy*dU*dbasJ*dbasIx+&
             dUy*dV*dbasJ*dbasIy + (dUy*dUx+dVy*dVx)*dbasJ*dbasI +&
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy) + Dii*32.0_DP*dnu**2*(dVy*dUx*dbasJy*dbasIx+ &
             0.5_DP*dVy*F*dbasJy*dbasIy+0.5_DP*F*dUx*dbasJx*dbasIx+ &
             0.25_DP*F**2*dbasJx*dbasIy)*gama   ) 
                         
      p_DlocalMatrixA12(jdofe,idofe,iel) = p_DlocalMatrixA12(jdofe,idofe,iel) + dval
      
      ! A21
      dval = p_DcubWeight(icubp,iel)*(  dU*dVx*dbasJx*dbasI+&
             dV*dVx*dbasJy*dbasI+(dUy*dUx+dVy*dVx)*dbasJ*dbasI+&
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy) + Dii*32.0_DP*dnu**2*(dVy*dUx*dbasJy*dbasIx+ &
             0.5_DP*dVy*F*dbasJy*dbasIy+0.5_DP*F*dUx*dbasJx*dbasIx+ &
             0.25_DP*F**2*dbasJx*dbasIy)*gama   )      
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
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    call ls_nonlinear_weight(Dii,rcollection,gama)    
    ! And a frequently used combination
    F = dUy+dVx
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA13%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasIx = p_DbasTestA13(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA13(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA13%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA13(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A13
      
      ! A31
      dval = -p_DcubWeight(icubp,iel)*gama*4.0_DP*dnu* &
             (   (2.0_DP*dUx**2+2.0_DP*dVy*dUx)*dbasJ*dbasIx + &
             (dUx*F+dVy*F)*dbasJ*dbasIy    )
      p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
           
      ! A23
      
      ! A32
      dval = -p_DcubWeight(icubp,iel)*gama*4.0_DP*dnu* &
             (   (2.0_DP*dVy**2+2.0_DP*dVy*dUx)*dbasJ*dbasIy + &
             (dUx*F+dVy*F)*dbasJ*dbasIx   )      
      p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + dval
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  
  ! ++++++++++++++++++++++++++++++++++++++
  ! Calculate blocks:   A14, A24, A41, A42
  ! A15, A25, A51, A52, A16, A26, A61, A62
  ! ++++++++++++++++++++++++++++++++++++++
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

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    call ls_nonlinear_weight(Dii,rcollection,gama)    
    ! And a frequently used combination
    F = dUy+dVx
    
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
                   
      ! A41
      dval = -p_DcubWeight(icubp,iel)*(  dUx*dbasJx*dbasI + 4.0_DP*gama*dnu* &
             (2.0_DP*dUx**2*dbasJ*dbasIx+dUx*F*dbasJ*dbasIy)   )
      p_DlocalMatrixA41(idofe,jdofe,iel) = p_DlocalMatrixA41(idofe,jdofe,iel) + dval

      ! A24
      
      ! A42
      dval = -p_DcubWeight(icubp,iel)*(  dUy*dbasJx*dbasI + 4.0_DP*gama*dnu* &
             (2.0_DP*dUx*dVy*dbasJ*dbasIy+dUx*F*dbasJ*dbasIx)   )      
      p_DlocalMatrixA42(idofe,jdofe,iel) = p_DlocalMatrixA42(idofe,jdofe,iel) + dval

      ! A15
      
      ! A51
      dval = -p_DcubWeight(icubp,iel)*(  dUx*dbasJy*dbasI+dVx*dbasJx*dbasI +&
             4.0_DP*gama*dnu* &
             (2.0_DP*dUx*F*dbasJ*dbasIx+F**2*dbasJ*dbasIy)   )  
      p_DlocalMatrixA51(idofe,jdofe,iel) = p_DlocalMatrixA51(idofe,jdofe,iel) + dval
                      
      ! A25
      
      ! A52
      dval = -p_DcubWeight(icubp,iel)*(  dUy*dbasJy*dbasI+dVy*dbasJx*dbasI + &
             4.0_DP*gama*dnu* &
             (2.0_DP*dVy*F*dbasJ*dbasIy+F**2*dbasJ*dbasIx)   ) 
      p_DlocalMatrixA52(idofe,jdofe,iel) = p_DlocalMatrixA52(idofe,jdofe,iel) + dval
      
      ! A16
      
      ! A61
      dval = -p_DcubWeight(icubp,iel)*(  dVx*dbasJy*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dVy*dUx*dbasJ*dbasIx+dVy*F*dbasJ*dbasIy)*gama   )       
      p_DlocalMatrixA61(idofe,jdofe,iel) = p_DlocalMatrixA61(idofe,jdofe,iel) + dval      
      
      ! A26
      
      ! A62
      dval = -p_DcubWeight(icubp,iel)*(  dVy*dbasJy*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dVy**2*dbasJ*dbasIy+dVy*F*dbasJ*dbasIx)*gama   )      
      p_DlocalMatrixA62(idofe,jdofe,iel) = p_DlocalMatrixA62(idofe,jdofe,iel) + dval
            
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  

 
  ! +++++++++++++++++++
  ! Calculate block A33
  ! +++++++++++++++++++
  
  

  ! +++++++++++++++++++++++++
  ! Calculate blocks A34, A43
  ! A35, A53, A36, A63 
  ! +++++++++++++++++++++++++
  
  
  
  ! +++++++++++++++++++++++++++++
  ! Calculate block A44, A55, A66
  !  A45, A54, A46, A64, A56, A65
  ! +++++++++++++++++++++++++++++

  
  end subroutine

 
  ! ***************************************************************************

!<subroutine>
  subroutine ls_svp2D_Matrix_un(RmatrixData,rassemblyData,rmatrixAssembly,&
             npointsPerElement,nelements,revalVectors,rcollection)
!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
  ! This is where the (BIG-BANG)**2 starts :D
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
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA12
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA13,p_DlocalMatrixA14
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA15,p_DlocalMatrixA16
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA21,p_DlocalMatrixA22
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA23,p_DlocalMatrixA24
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA25,p_DlocalMatrixA26  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA33,p_DlocalMatrixA34  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA36   
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA41
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA43,p_DlocalMatrixA44  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA45  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA51,p_DlocalMatrixA52
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA54  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA55,p_DlocalMatrixA56  
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA62
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA63 
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA65,p_DlocalMatrixA66  
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA33,p_DbasTestA33
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA44,p_DbasTestA44
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA14,p_DbasTestA14
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA34,p_DbasTestA34

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11,p_rmatrixDataA33,p_rmatrixDataA44
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA13,p_rmatrixDataA14,p_rmatrixDataA34
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

  ! Scaling factors
  real(DP) :: beta, gama
  
  ! Continuity equation scaling
  real(DP) :: alpha
  
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dUx, dUy, dVx, dVy
  
  ! Viscosity and its 1st derivative
  ! 2nd invariant and a frequently used value
  real(DP) :: nu, dnu, Dii, F
  
  ! Continuity equation scaling factor ---> \alpha
  !  default = 1.0_DP  no scaling
  alpha = rcollection%DquickAccess(3)

  ! Linearization Scheme ---> \beta
  beta = rcollection%DquickAccess(4)

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
  p_DlocalMatrixA15 => RmatrixData(1,5)%p_Dentry
  p_DlocalMatrixA16 => RmatrixData(1,6)%p_Dentry

  p_DlocalMatrixA21 => RmatrixData(2,1)%p_Dentry
  p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
  p_DlocalMatrixA23 => RmatrixData(2,3)%p_Dentry
  p_DlocalMatrixA24 => RmatrixData(2,4)%p_Dentry
  p_DlocalMatrixA25 => RmatrixData(2,5)%p_Dentry
  p_DlocalMatrixA26 => RmatrixData(2,6)%p_Dentry
  
  p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
  p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
  p_DlocalMatrixA33 => RmatrixData(3,3)%p_Dentry
  p_DlocalMatrixA34 => RmatrixData(3,4)%p_Dentry
  p_DlocalMatrixA36 => RmatrixData(3,6)%p_Dentry   

  p_DlocalMatrixA41 => RmatrixData(4,1)%p_Dentry
  p_DlocalMatrixA43 => RmatrixData(4,3)%p_Dentry
  p_DlocalMatrixA44 => RmatrixData(4,4)%p_Dentry
  p_DlocalMatrixA45 => RmatrixData(4,5)%p_Dentry   
  
  p_DlocalMatrixA51 => RmatrixData(5,1)%p_Dentry
  p_DlocalMatrixA52 => RmatrixData(5,2)%p_Dentry
  p_DlocalMatrixA54 => RmatrixData(5,4)%p_Dentry
  p_DlocalMatrixA55 => RmatrixData(5,5)%p_Dentry
  p_DlocalMatrixA56 => RmatrixData(5,6)%p_Dentry   
  
  p_DlocalMatrixA62 => RmatrixData(6,2)%p_Dentry
  p_DlocalMatrixA63 => RmatrixData(6,3)%p_Dentry
  p_DlocalMatrixA65 => RmatrixData(6,5)%p_Dentry
  p_DlocalMatrixA66 => RmatrixData(6,6)%p_Dentry  
  
  ! Get the velocity field from the parameters
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! We do the calculation in a block-by-block manner. All the
  ! relevant blocks will be calculated in the same loop over the
  ! elements.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Un-weighted formulation
  gama = 1.0_DP
  
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
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    ! And a frequently used combination
    F = dUy+dVx
    
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
      dval = p_DcubWeight(icubp,iel)*beta*(   dU*dUx*dbasJx*dbasI+&
             dV*dUx*dbasJy*dbasI + 4.0_DP*gama*&
             nu*dnu*(4.0_DP*dUx**2*dbasJx*dbasIx+ &
             2.0_DP*dUx*F*dbasJx*dbasIy+2.0_DP*F*dUx*dbasJy*dbasIx+ &
             F**2*dbasJy*dbasIy)   )

      dval = dval + p_DcubWeight(icubp,iel)*(   dU**2*dbasJx*dbasIx+dU*dV*dbasJx*dbasIy + &
             dV*dU*dbasJy*dbasIx + dV**2*dbasJy*dbasIy + &
             gama*nu**2*(4.0_DP*dbasJx*dbasIx + 2.0_DP*dbasJy*dbasIy) + & 
             alpha*(dbasJx*dbasIx)  )           
             
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
      ! A22
      dval = p_DcubWeight(icubp,iel)*beta*(   dU*dVy*dbasJx*dbasI+&
             dV*dVy*dbasJy*dbasI + 4.0_DP*gama* &
             nu*dnu*(4.0_DP*dVy**2*dbasJy*dbasIy+ &
             2.0_DP*dVy*F*dbasJy*dbasIx+2.0_DP*F*dVy*dbasJx*dbasIy+ &
             F**2*dbasJx*dbasIx)   )
      
      dval = dval + p_DcubWeight(icubp,iel)*(   dU**2*dbasJx*dbasIx+dU*dV*dbasJy*dbasIx + &
             dV*dU*dbasJx*dbasIy + dV**2*dbasJy*dbasIy + &
             gama*nu**2*(4.0_DP*dbasJy*dbasIy + 2.0_DP*dbasJx*dbasIx) + &
             alpha*(dbasJy*dbasIy)  )         
              
      p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
      
      ! A12
      dval = p_DcubWeight(icubp,iel)*beta*(   dU*dVx*dbasJx*dbasI+&
             dV*dVx*dbasJy*dbasI + & 
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy)   ) 
             
      dval = dval + p_DcubWeight(icubp,iel)*(   gama*nu**2*(2.0_DP*dbasJx*dbasIy)+& 
             alpha*(dbasJy*dbasIx)  ) 
            
      p_DlocalMatrixA12(jdofe,idofe,iel) = p_DlocalMatrixA12(jdofe,idofe,iel) + dval
      
      ! A21
      dval = p_DcubWeight(icubp,iel)*beta*(   dUy*dU*dbasJ*dbasIx + & 
             dUy*dV*dbasJ*dbasIy + & 
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy)   ) 
             
      dval = dval + p_DcubWeight(icubp,iel)*(   gama*nu**2*(2.0_DP*dbasJx*dbasIy)+& 
             alpha*(dbasJy*dbasIx)  )
                   
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
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu) 
    ! And a frequently used combination
    F = dUy+dVx
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA13%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasIx = p_DbasTestA13(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA13(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA13%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA13(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A13
      dval = -p_DcubWeight(icubp,iel)*gama*beta*4.0_DP*dnu* &
             (   (2.0_DP*dUx**2+2.0_DP*dVy*dUx)*dbasJ*dbasIx + &
             (dUx*F+dVy*F)*dbasJ*dbasIy    )
             
      dval = dval - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIx   )
      
      p_DlocalMatrixA13(jdofe,idofe,iel) = p_DlocalMatrixA13(jdofe,idofe,iel) + dval
      
      ! A31
      dval = - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIx   )      
      p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
           
      ! A23
      dval = -p_DcubWeight(icubp,iel)*gama*beta*4.0_DP*dnu* &
             (   (2.0_DP*dVy**2+2.0_DP*dVy*dUx)*dbasJ*dbasIy + &
             (dUx*F+dVy*F)*dbasJ*dbasIx   )

      dval = dval - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIy   )           
             
      p_DlocalMatrixA23(jdofe,idofe,iel) = p_DlocalMatrixA23(jdofe,idofe,iel) + dval
      
      ! A32
      dval = - p_DcubWeight(icubp,iel)*(   2.0_DP*gama*nu*dbasJ*dbasIy   )       
      p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + dval
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  
  ! ++++++++++++++++++++++++++++++++++++++
  ! Calculate blocks:   A14, A24, A41, A42
  ! A15, A25, A51, A52, A16, A26, A61, A62
  ! ++++++++++++++++++++++++++++++++++++++
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

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)  
    ! And a frequently used combination
    F = dUy+dVx
    
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
      dval = -p_DcubWeight(icubp,iel) * beta* (   dUx*dbasJx*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dUx**2*dbasJ*dbasIx+dUx*F*dbasJ*dbasIy)  )
             
      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )             
             
      p_DlocalMatrixA14(jdofe,idofe,iel) = p_DlocalMatrixA14(jdofe,idofe,iel) + dval
      
      ! A41             
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )
                   
      p_DlocalMatrixA41(idofe,jdofe,iel) = p_DlocalMatrixA41(idofe,jdofe,iel) + dval

      ! A24
      dval = -p_DcubWeight(icubp,iel) * beta*(   dUy*dbasJx*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dUx*dVy*dbasJ*dbasIy+dUx*F*dbasJ*dbasIx)*gama   )
      p_DlocalMatrixA24(jdofe,idofe,iel) = p_DlocalMatrixA24(jdofe,idofe,iel) + dval
      
      ! A42

      ! A15
      dval = -p_DcubWeight(icubp,iel)*beta*(   dUx*dbasJy*dbasI+dVx*dbasJx*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dUx*F*dbasJ*dbasIx+F**2*dbasJ*dbasIy)   )
      
      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )
             
      p_DlocalMatrixA15(jdofe,idofe,iel) = p_DlocalMatrixA15(jdofe,idofe,iel) + dval
      
      ! A51
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )
      p_DlocalMatrixA51(idofe,jdofe,iel) = p_DlocalMatrixA51(idofe,jdofe,iel) + dval
                      
      ! A25
      dval = -p_DcubWeight(icubp,iel)*beta*(   dUy*dbasJy*dbasI+dVy*dbasJx*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dVy*F*dbasJ*dbasIy+F**2*dbasJ*dbasIx)   )

      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )           
             
      p_DlocalMatrixA25(jdofe,idofe,iel) = p_DlocalMatrixA25(jdofe,idofe,iel) + dval
      
      ! A52
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJx*dbasIx+dV*dbasJx*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIx   )  
      p_DlocalMatrixA52(idofe,jdofe,iel) = p_DlocalMatrixA52(idofe,jdofe,iel) + dval
      
      ! A16
      dval = -p_DcubWeight(icubp,iel) * beta*(   dVx*dbasJy*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dVy*dUx*dbasJ*dbasIx+dVy*F*dbasJ*dbasIy)*gama   )
      p_DlocalMatrixA16(jdofe,idofe,iel) = p_DlocalMatrixA16(jdofe,idofe,iel) + dval
      
      ! A61   
      
      ! A26
      dval = -p_DcubWeight(icubp,iel)*beta*(   dVy*dbasJy*dbasI + &
             4.0_DP*gama*dnu*(2.0_DP*dVy**2*dbasJ*dbasIy+dVy*F*dbasJ*dbasIx)   )
             
      dval = dval - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )            
             
      p_DlocalMatrixA26(jdofe,idofe,iel) = p_DlocalMatrixA26(jdofe,idofe,iel) + dval
      
      ! A62  
      dval = - p_DcubWeight(icubp,iel)*(   dU*dbasJy*dbasIx+dV*dbasJy*dbasIy + &
             2.0_DP*gama*nu*dbasJ*dbasIy   )       
      p_DlocalMatrixA62(idofe,jdofe,iel) = p_DlocalMatrixA62(idofe,jdofe,iel) + dval      
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

    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)

    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA33%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA33(idofe,DER_FUNC,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA33%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA33(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A33
      dval = p_DcubWeight(icubp,iel) * 2.0_DP*gama*dbasJ*dbasI
      p_DlocalMatrixA33(jdofe,idofe,iel) = p_DlocalMatrixA33(jdofe,idofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  

  ! +++++++++++++++++++++++++
  ! Calculate blocks A34, A43
  ! A35, A53, A36, A63 
  ! +++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA34%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasI = p_DbasTestA34(idofe,DER_FUNC,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA34%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA34(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A34
      dval = p_DcubWeight(icubp,iel) * gama*dbasJ*dbasI
      p_DlocalMatrixA34(jdofe,idofe,iel) = p_DlocalMatrixA34(jdofe,idofe,iel) + dval

      ! A43
      p_DlocalMatrixA43(idofe,jdofe,iel) = p_DlocalMatrixA43(idofe,jdofe,iel) + dval

      ! A35
      ! A53

      ! A36
      p_DlocalMatrixA36(jdofe,idofe,iel) = p_DlocalMatrixA36(jdofe,idofe,iel) + dval

      ! A63
      p_DlocalMatrixA63(idofe,jdofe,iel) = p_DlocalMatrixA63(idofe,jdofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  
  
  ! +++++++++++++++++++++++++++++
  ! Calculate block A44, A55, A66
  !  A45, A54, A46, A64, A56, A65
  ! +++++++++++++++++++++++++++++
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)

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
      dval = p_DcubWeight(icubp,iel) * (  dbasJx*dbasIx + gama*dbasJ*dbasI  )
      p_DlocalMatrixA44(jdofe,idofe,iel) = p_DlocalMatrixA44(jdofe,idofe,iel) + dval

      ! A55
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIy + dbasJx*dbasIx + &
            2.0_DP*gama*dbasJ*dbasI  )
      p_DlocalMatrixA55(jdofe,idofe,iel) = p_DlocalMatrixA55(jdofe,idofe,iel) + dval

      ! A66
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIy + gama*dbasJ*dbasI  )
      p_DlocalMatrixA66(jdofe,idofe,iel) = p_DlocalMatrixA66(jdofe,idofe,iel) + dval
      
      ! A45
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIx  )
      p_DlocalMatrixA45(jdofe,idofe,iel) = p_DlocalMatrixA45(jdofe,idofe,iel) + dval  

      ! A54
      p_DlocalMatrixA54(idofe,jdofe,iel) = p_DlocalMatrixA54(idofe,jdofe,iel) + dval

      ! A46
      ! A64

      ! A56
      p_DlocalMatrixA56(jdofe,idofe,iel) = p_DlocalMatrixA56(jdofe,idofe,iel) + dval  

      ! A65
      p_DlocalMatrixA65(idofe,jdofe,iel) = p_DlocalMatrixA65(idofe,jdofe,iel) + dval

      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  end subroutine

 
  ! ***************************************************************************

!<subroutine>
  subroutine ls_svp2D_NewtonMatrix_un(RmatrixData,rassemblyData,rmatrixAssembly,&
             npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Assemble the solution matrix in a block-by-block procedures.
  ! This is where the (BIG-BANG)**2 starts :D
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
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA12
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA21,p_DlocalMatrixA22
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA41,p_DlocalMatrixA42
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA51,p_DlocalMatrixA52
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA61,p_DlocalMatrixA62
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA14,p_DbasTestA14

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA13,p_rmatrixDataA14
  
  ! Known velocity data
  real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2

  ! Scaling factors
  real(DP) :: gama
  
  ! Velocity values/derivatives in cubature points 
  real(DP) :: dU, dV, dUx, dUy, dVx, dVy
  
  ! Viscosity and its 1st derivative
  ! 2nd invariant and a frequently used value
  real(DP) :: nu, dnu, Dii, F
  

  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
  p_rmatrixDataA13 => RmatrixData(1,3)
  p_rmatrixDataA14 => RmatrixData(1,4)


  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
  p_DbasTrialA13 => RmatrixData(1,3)%p_DbasTrial
  p_DbasTestA13 => RmatrixData(1,3)%p_DbasTest
  p_DbasTrialA14 => RmatrixData(1,4)%p_DbasTrial
  p_DbasTestA14 => RmatrixData(1,4)%p_DbasTest
  
  
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
  p_DlocalMatrixA12 => RmatrixData(1,2)%p_Dentry

  p_DlocalMatrixA21 => RmatrixData(2,1)%p_Dentry
  p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
  
  p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
  p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
 
  p_DlocalMatrixA41 => RmatrixData(4,1)%p_Dentry
  p_DlocalMatrixA42 => RmatrixData(4,2)%p_Dentry
  
  p_DlocalMatrixA51 => RmatrixData(5,1)%p_Dentry
  p_DlocalMatrixA52 => RmatrixData(5,2)%p_Dentry
  
  p_DlocalMatrixA61 => RmatrixData(6,1)%p_Dentry
  p_DlocalMatrixA62 => RmatrixData(6,2)%p_Dentry

  
  ! Get the velocity field from the parameters
  p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
  p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! We do the calculation in a block-by-block manner. All the
  ! relevant blocks will be calculated in the same loop over the
  ! elements.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Un-weighted formulation
  gama = 1.0_DP
  
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
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)
    call ls_viscosity_model_der(Dii,rcollection,dnu)
    
    ! And a frequently used combination
    F = dUy+dVx
    
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
      dval = p_DcubWeight(icubp,iel)*(  dUx*dU*dbasJ*dbasIx+&
             dUx*dV*dbasJ*dbasIy + (dUx**2+dVx**2)*dbasJ*dbasI +&
             4.0_DP*gama*nu*dnu*(4.0_DP*dUx**2*dbasJx*dbasIx+ &
             2.0_DP*dUx*F*dbasJx*dbasIy+2.0_DP*F*dUx*dbasJy*dbasIx+ &
             F**2*dbasJy*dbasIy) + Dii*32.0_DP*dnu**2*(dUx**2*dbasJx*dbasIx+ &
             0.5_DP*dUx*F*dbasJx*dbasIy+0.5_DP*F*dUx*dbasJy*dbasIx+ &
             0.25_DP*F**2*dbasJy*dbasIy)*gama    )         
             
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
      
      ! A22
      dval = p_DcubWeight(icubp,iel)*(  dVy*dU*dbasJ*dbasIx+&
             dVy*dV*dbasJ*dbasIy + (dUy**2+dVy**2)*dbasJ*dbasI + &
             4.0_DP*gama*nu*dnu*(4.0_DP*dVy**2*dbasJy*dbasIy+ &
             2.0_DP*dVy*F*dbasJy*dbasIx+2.0_DP*F*dVy*dbasJx*dbasIy+ &
             F**2*dbasJx*dbasIx) + Dii*32.0_DP*dnu**2*(dVy**2*dbasJy*dbasIy+ &
             0.5_DP*dVy*F*dbasJy*dbasIx+0.5_DP*F*dVy*dbasJx*dbasIy+ &
             0.25_DP*F**2*dbasJx*dbasIx)*gama   )      
              
      p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
      
      ! A12
      dval = p_DcubWeight(icubp,iel)*(  dUy*dU*dbasJ*dbasIx+&
             dUy*dV*dbasJ*dbasIy + (dUy*dUx+dVy*dVx)*dbasJ*dbasI +&
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy) + Dii*32.0_DP*dnu**2*(dVy*dUx*dbasJy*dbasIx+ &
             0.5_DP*dVy*F*dbasJy*dbasIy+0.5_DP*F*dUx*dbasJx*dbasIx+ &
             0.25_DP*F**2*dbasJx*dbasIy)*gama   ) 
                         
      p_DlocalMatrixA12(jdofe,idofe,iel) = p_DlocalMatrixA12(jdofe,idofe,iel) + dval
      
      ! A21
      dval = p_DcubWeight(icubp,iel)*(  dU*dVx*dbasJx*dbasI+&
             dV*dVx*dbasJy*dbasI+(dUy*dUx+dVy*dVx)*dbasJ*dbasI+&
             4.0_DP*nu*dnu*gama*(4.0_DP*dVy*dUx*dbasJy*dbasIx+ &
             2.0_DP*dVy*F*dbasJy*dbasIy+2.0_DP*F*dUx*dbasJx*dbasIx+ &
             F**2*dbasJx*dbasIy) + Dii*32.0_DP*dnu**2*(dVy*dUx*dbasJy*dbasIx+ &
             0.5_DP*dVy*F*dbasJy*dbasIy+0.5_DP*F*dUx*dbasJx*dbasIx+ &
             0.25_DP*F**2*dbasJx*dbasIy)*gama   )      
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
    dUx = p_Du1(icubp,iel,DER_DERIV2D_X)
    dVx = p_Du2(icubp,iel,DER_DERIV2D_X)

    dUy = p_Du1(icubp,iel,DER_DERIV2D_Y)
    dVy = p_Du2(icubp,iel,DER_DERIV2D_Y)
    
    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)
  
    ! And a frequently used combination
    F = dUy+dVx
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rmatrixDataA13%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      dbasIx = p_DbasTestA13(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA13(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA13%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJ = p_DbasTrialA13(jdofe,DER_FUNC,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A13
      
      ! A31
      dval = -p_DcubWeight(icubp,iel)*gama*4.0_DP*dnu* &
             (   (2.0_DP*dUx**2+2.0_DP*dVy*dUx)*dbasJ*dbasIx + &
             (dUx*F+dVy*F)*dbasJ*dbasIy    )
      p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
           
      ! A23
      
      ! A32
      dval = -p_DcubWeight(icubp,iel)*gama*4.0_DP*dnu* &
             (   (2.0_DP*dVy**2+2.0_DP*dVy*dUx)*dbasJ*dbasIy + &
             (dUx*F+dVy*F)*dbasJ*dbasIx   )      
      p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + dval
                      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  
  ! ++++++++++++++++++++++++++++++++++++++
  ! Calculate blocks:   A14, A24, A41, A42
  ! A15, A25, A51, A52, A16, A26, A61, A62
  ! ++++++++++++++++++++++++++++++++++++++
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

    ! Calculate viscosity and its 1st derivative
    ! we need to evaluate the 2nd invariant of the 
    ! deformation rate tensor Dii(u)
    Dii = 1.0_DP/2.0_DP * &
              (4.0_DP*dUx**2 + 2.0_DP*(dUy+dVx)**2 + 4.0_DP*dVy**2)
    call ls_viscosity_model(Dii,rcollection,nu)              
    call ls_viscosity_model_der(Dii,rcollection,dnu)
  
    ! And a frequently used combination
    F = dUy+dVx
    
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
                   
      ! A41
      dval = -p_DcubWeight(icubp,iel)*(  dUx*dbasJx*dbasI + 4.0_DP*gama*dnu* &
             (2.0_DP*dUx**2*dbasJ*dbasIx+dUx*F*dbasJ*dbasIy)   )
      p_DlocalMatrixA41(idofe,jdofe,iel) = p_DlocalMatrixA41(idofe,jdofe,iel) + dval

      ! A24
      
      ! A42
      dval = -p_DcubWeight(icubp,iel)*(  dUy*dbasJx*dbasI + 4.0_DP*gama*dnu* &
             (2.0_DP*dUx*dVy*dbasJ*dbasIy+dUx*F*dbasJ*dbasIx)   )      
      p_DlocalMatrixA42(idofe,jdofe,iel) = p_DlocalMatrixA42(idofe,jdofe,iel) + dval

      ! A15
      
      ! A51
      dval = -p_DcubWeight(icubp,iel)*(  dUx*dbasJy*dbasI+dVx*dbasJx*dbasI +&
             4.0_DP*gama*dnu* &
             (2.0_DP*dUx*F*dbasJ*dbasIx+F**2*dbasJ*dbasIy)   )  
      p_DlocalMatrixA51(idofe,jdofe,iel) = p_DlocalMatrixA51(idofe,jdofe,iel) + dval
                      
      ! A25
      
      ! A52
      dval = -p_DcubWeight(icubp,iel)*(  dUy*dbasJy*dbasI+dVy*dbasJx*dbasI + &
             4.0_DP*gama*dnu* &
             (2.0_DP*dVy*F*dbasJ*dbasIy+F**2*dbasJ*dbasIx)   ) 
      p_DlocalMatrixA52(idofe,jdofe,iel) = p_DlocalMatrixA52(idofe,jdofe,iel) + dval
      
      ! A16
      
      ! A61
      dval = -p_DcubWeight(icubp,iel)*(  dVx*dbasJy*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dVy*dUx*dbasJ*dbasIx+dVy*F*dbasJ*dbasIy)*gama   )       
      p_DlocalMatrixA61(idofe,jdofe,iel) = p_DlocalMatrixA61(idofe,jdofe,iel) + dval      
      
      ! A26
      
      ! A62
      dval = -p_DcubWeight(icubp,iel)*(  dVy*dbasJy*dbasI + 4.0_DP*dnu* &
             (2.0_DP*dVy**2*dbasJ*dbasIy+dVy*F*dbasJ*dbasIx)*gama   )      
      p_DlocalMatrixA62(idofe,jdofe,iel) = p_DlocalMatrixA62(idofe,jdofe,iel) + dval
            
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  

 
  ! +++++++++++++++++++
  ! Calculate block A33
  ! +++++++++++++++++++
  
  

  ! +++++++++++++++++++++++++
  ! Calculate blocks A34, A43
  ! A35, A53, A36, A63 
  ! +++++++++++++++++++++++++
  
  
  
  ! +++++++++++++++++++++++++++++
  ! Calculate block A44, A55, A66
  !  A45, A54, A46, A64, A56, A65
  ! +++++++++++++++++++++++++++++

  
  end subroutine



  !****************************************************************************
  
!<subroutine>

  subroutine ls_Defect(rmatrix, rvector, rrhs, rdefect)
                
 !<description>  
  ! This subroutine calculates the defect on the finest level.
  ! The defect is calculated based on the rmatrix, rrhs and rvector
  !  into the vector rdefect.
 !</description>                

  !<input>
  ! Solution vector
  type(t_vectorBlock), intent(in) :: rvector
  
  ! RHS vector
  type(t_vectorBlock), intent(in) :: rrhs  

  ! The system block matrix.
  type(t_matrixBlock), intent(in) :: rmatrix
  !</input>


  !<output>    
  ! Defect vector
  type(t_vectorBlock), intent(out) :: rdefect  
  !</output>

!</subroutine>

  ! The defect is calculated simply!
  call lsysbl_copyVector (rrhs,rdefect)
  call lsysbl_blockMatVec (rmatrix, rvector, rdefect, 1.0_DP, -1.0_DP)    
    
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



  !****************************************************************************

!<subroutine>

  subroutine ls_L2_Norm(dintvalue,rassemblyData,rintegralAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Calculates the (squared) L2 error of an arbitrary FEM function v
  ! to the analytical function u:
  !
  !  <tex>  || v - u ||^2_{L2}  </tex>
  !
  ! The FEM function must be provided in revalVectors.
  ! The routine only supports non-interleaved vectors.
!</description>

!<input>
  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

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

!<output>
  ! Returns the value of the integral
  real(DP), dimension(:), intent(out) :: dintvalue
!</output>  

!<subroutine>

  ! Local variables
  real(DP) :: dval1,dval2,dx,dy, dC
  integer :: iel, icubp
  real(DP), dimension(:,:), pointer :: p_DcubWeight
  real(DP), dimension(:,:), pointer :: p_Dfunc,p_DderivX,p_DderivY
  real(DP), dimension(:,:,:), pointer :: p_Dpoints
  real(DP), dimension(:), pointer :: p_DelementArea
  real(DP) :: dpressure_in, darea_in, dpressure_o, darea_o,r
  
  ! Cancel if no FEM function is given.
  if (revalVectors%ncount .eq. 0) then
    call output_line ("FEM function missing.",&
      OU_CLASS_ERROR,OU_MODE_STD,"ls_L2_Norm")
    call sys_halt()
  end if

  ! Skip interleaved vectors.
  if (revalVectors%p_RvectorData(1)%bisInterleaved) return

  ! The constant in the pressure filed
  dC = rcollection%DquickAccess(9)

  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  
  ! Get the coordinates of the cubature points
  p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
  
  ! Get the data array with the values of the FEM function
  ! in the cubature points
  p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)

  dintvalue = 0.0_DP

  select case (rcollection%IquickAccess(7))
  
  case (0)
    !X-Velocity
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the analytic function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = 2.0_DP*dx**2*(1.0_DP-dx)**2*(dy*(1.0_DP-dy)**2 - dy**2*(1.0_DP-dy))
      dval1 = 0.0_DP
      
      ! Get the error of the FEM function to the analytic function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel
    
  case (1)
    ! Y-Velocity
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = -2.0_DP*dy**2*(1.0_DP-dy)**2*(dx*(1.0_DP-dx)**2 - dx**2*(1.0_DP-dx))
      dval1 = 0.0_DP

      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel
  
  case (2)
    ! Pressure
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = dC*(dx**3 - dy**3)

      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel

  case (3)
    ! Sxx
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = -dC*(dx**3 - dy**3) + 8.0_DP*dx*dy*(2.0_DP*dx - 1.0_DP) &
             *(2.0_DP*dy - 1.0_DP)*(dx - 1.0_DP)*(dy - 1.0_DP)

      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel    

  case (4)
    ! Sxy
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = 2.0_DP*dx**2*(dx-1.0_DP)**2*(6.0_DP*dy**2-6.0_DP*dy+1.0_DP) - &
          2.0_DP*dy**2*(dy-1.0_DP)**2*(6.0_DP*dx**2-6.0_DP*dx+1.0_DP)

      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel   

  case (5)
    ! Syy
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = -dC*(dx**3 - dy**3) - 8.0_DP*dx*dy*(2.0_DP*dx - 1.0_DP) &
             *(2.0_DP*dy - 1.0_DP)*(dx - 1.0_DP)*(dy - 1.0_DP)

      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel  

  case (6)
    ! Static bubble error analysis
    call storage_getbase_double (&
     rintegralAssembly%p_rtriangulation%h_DelementVolume, p_DelementArea)
      
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement
      
      ! Get the real coordinate of the cubarure points
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      ! Calculate the distance of the cubator point to origin at (0,0)
      r = sqrt(dx**2 + dy**2)
      
      ! Check if the cubature point lies inside of the bubble
      if (r .le. 0.3_DP) then
      
        ! Area of each element inside of the bubble
        darea_in = p_DelementArea(rassemblyData%P_IelementList(iel))
        ! Total area of the bubble (numerical)
        ! The exact value is: 0.1963495408493621_DP
        rcollection%DquickAccess(12) = rcollection%DquickAccess(12) + darea_in
        ! Pressure value of each element inside of the bubble
        !    (normalized with the min value)          
        dpressure_in = p_Dfunc(icubp,iel) - rcollection%DquickAccess(16)
        ! Pressure mutiplied by the area of each element inside of the bubble      
        rcollection%DquickAccess(13) = rcollection%DquickAccess(13) + &
                                                     darea_in*dpressure_in
                                                     
     else
     
        ! outside the bubble 
        ! Area of each element outside of the bubble
        darea_o = p_DelementArea(rassemblyData%P_IelementList(iel))
        ! Total area of the bubble (numerical)
        ! The exact value is: 0.1963495408493621_DP
        rcollection%DquickAccess(14) = rcollection%DquickAccess(14) + darea_o
        ! Pressure value of each element outside of the bubble
        !    (normalized with the min value)
        dpressure_o = p_Dfunc(icubp,iel) - rcollection%DquickAccess(16)
        ! Pressure mutiplied by the area of each element outside of the bubble
        rcollection%DquickAccess(15) = rcollection%DquickAccess(15) + &
                                                     darea_o*dpressure_o     
                                                           
     end if
      
    end do ! icubp
    
    end do ! iel 

  case (7)
    ! Velocity divergence L^2-error
    ! Loop over the elements in the current set.
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_DderivX => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_X)
    p_DderivY => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV2D_Y)
    
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! The analytic function is not required here
      dval1 = 0.0_DP

      ! Get the error of the FEM function to the analytic function
      dval2 = p_DderivX(icubp,iel) + p_DderivY(icubp,iel)
      
      ! Maximum norm of the divergence
      if (dval2 .gt. rcollection%DquickAccess(11)) then
        rcollection%DquickAccess(11) = dval2
      end if  
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel 
  
  end select   
    
  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine ls_H1_Norm(dintvalue,rassemblyData,rintegralAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Calculates the (squared) H1 error of an arbitrary FEM function v
  ! to a function u (based on the H1 semi-norm).
  !
  !  <tex>  | v - u |^2_{H1}  </tex>
  !
  ! The FEM function must be provided in revalVectors.
  ! The routine only supports non-interleaved vectors.
!</description>

!<input>
  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

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

!<output>
  ! Returns the value of the integral
  real(DP), dimension(:), intent(out) :: dintvalue
!</output>  

!<subroutine>

  ! Local variables
  real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy, dC
  integer :: iel, icubp
  real(DP), dimension(:,:), pointer :: p_DcubWeight
  real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY
  real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
  ! Calcel if no FEM function is given.
  if (revalVectors%ncount .eq. 0) then
    call output_line ("FEM function missing.",&
      OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleH1error")
    call sys_halt()
  end if

  ! Skip interleaved vectors.
  if (revalVectors%p_RvectorData(1)%bisInterleaved) return

  ! Skip vector-valued FE functions
  if (revalVectors%p_RvectorData(1)%ndimfe .ne. 1) return


  ! The constant in the pressure filed
  dC = rcollection%DquickAccess(9)

  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  
  ! Get the coordinates of the cubature points
  p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
  
  ! Get the data array with the values of the FEM function
  ! in the cubature points
  p_DderivX => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_X)
  p_DderivY => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_Y)

  dintvalue = 0.0_DP


  select case (rcollection%IquickAccess(7))
  
  case (0)
    !X-Velocity  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = 4.0_DP*dx*dy*(2*dx-1.0_DP)*(2.0_DP*dy-1.0_DP)*(dx-1.0_DP)*(dy-1.0_DP)
      dderivY1 = 2.0_DP*dx**2*(dx-1.0_DP)**2*(6.0_DP*dy**2-6.0_DP*dy+1.0_DP)
      dderivX1 = 0.0_DP
      dderivY1 = 0.0_DP
      ! Get the error of the FEM function derivatives of the bubble function
      ! in the cubature point
      dderivX2 = p_DderivX(icubp,iel)
      dderivY2 = p_DderivY(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) H1 error:
      dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
        ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
      
    end do ! icubp
    
    end do ! iel

  case (1)
    ! Y-Velocity
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = -2.0_DP*dy**2*(dy-1.0_DP)**2*(6.0_DP*dx**2-6.0_DP*dx+1.0_DP)
      dderivY1 = -4.0_DP*dx*dy*(2*dx-1.0_DP)*(2.0_DP*dy-1.0_DP)*(dx-1.0_DP)*(dy-1.0_DP)
      dderivX1 = 0.0_DP
      dderivY1 = 0.0_DP
      
      ! Get the error of the FEM function derivatives of the bubble function
      ! in the cubature point
      dderivX2 = p_DderivX(icubp,iel)
      dderivY2 = p_DderivY(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) H1 error:
      dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
        ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
      
    end do ! icubp
    
    end do ! iel
      
  case (2)
    ! Pressure    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = 3.0_DP*dC*dx**2
      dderivY1 = -3.0_DP*dC*dy**2

      ! Get the error of the FEM function derivatives of the bubble function
      ! in the cubature point
      dderivX2 = p_DderivX(icubp,iel)
      dderivY2 = p_DderivY(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) H1 error:
      dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
        ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
      
    end do ! icubp
    
    end do ! iel

  case (3)
    ! Vorticity
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = -4.0_DP*dx*(dx-1.0_DP)**2*(6.0_DP*dy**2-6.0_DP*dy+1.0_DP)-&
          2.0_DP*dy**2*(12.0_DP*dx-6.0_DP)*(dy-1.0_DP)**2-&
          2.0_DP*dx**2*(2.0_DP*dx-2.0_DP)*(6.0_DP*dy**2-6.0_DP*dy+1.0_DP)
      dderivY1 = -4.0_DP*dy*(dy-1.0_DP)**2*(6.0_DP*dx**2-6.0_DP*dx+1.0_DP)-&
          2.0_DP*dx**2*(12.0_DP*dy-6.0_DP)*(dx-1.0_DP)**2-&
          2.0_DP*dy**2*(2.0_DP*dy-2.0_DP)*(6.0_DP*dx**2-6.0_DP*dx+1.0_DP) 

      ! Get the error of the FEM function derivatives of the bubble function
      ! in the cubature point
      dderivX2 = p_DderivX(icubp,iel) 
      dderivY2 = p_DderivY(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) H1 error:
      dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
        ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
      
    end do ! icubp
    
    end do ! iel
      
  end select    
    
  end subroutine



  !****************************************************************************

!<subroutine>

  subroutine ls_calc_min(dintvalue,rassemblyData,rintegralAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
  ! Calculates the minimum value of a FEM function.
  ! The rcollection%DquickAccess(16) must be initiaized (high value)
  !  before calling this function. The final minimum value will 
  !  be stored in rcollection%DquickAccess(16) too.
  ! The FEM function must be provided in revalVectors.
!</description>

!<input>
  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

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

!<output>
  ! Returns the value of the integral
  real(DP), dimension(:), intent(out) :: dintvalue
!</output>  

!<subroutine>

  ! Local variables
  integer :: iel, icubp
  real(DP), dimension(:,:), pointer :: p_Dfunc
  
  ! Cancel if no FEM function is given.
  if (revalVectors%ncount .eq. 0) then
    call output_line ("FEM function missing.",&
      OU_CLASS_ERROR,OU_MODE_STD,"ls_L2_Norm")
    call sys_halt()
  end if
  
  ! Get the data array with the values of the FEM function
  ! in the cubature points
  p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)

  dintvalue = 0.0_DP

  ! Static bubble minimum pressure calculation
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement
        rcollection%DquickAccess(16) = &
               min(rcollection%DquickAccess(16),p_Dfunc(icubp,iel))
    end do ! icubp
  
  end do ! iel  
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine jstab_reacJumpStabil2d_mod ( &
      rmatrixScalar,dgamma,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the reactive jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \gamma \nu 1/|E| \int_E [u] [v] ds$$
  ! </tex>
  !
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
  !
  ! WARNING: For edge oriented stabilisation, the underlying matrix rmatrix
  !   must have an extended stencil! The matrix structure must be set up with
  !   the BILF_MATC_EDGEBASED switch!!!
!</description>

!<input>
  ! Stabilisation parameter. Standard=0.01
  real(DP), intent(in) :: dgamma

  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta

  ! 1D cubature formula to use for line integration
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType

  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dval,dedgelength,dedgeweight,dphi,dpsi,dcoeff

  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTempl
  integer, dimension(EL_MAXNBAS*2), target :: Idofs

  ! Arrays saving the local DOF numbers belonging to the global
  ! DOF numbers in Idofs
  integer, dimension(EL_MAXNBAS*2),target :: IlocalDofs

  ! Renumbering strategy for local DOF`s
  integer, dimension(EL_MAXNBAS), target :: IlocalDofRenum

  ! Number of local DOF`s on the patch
  integer :: ndof

  ! Number of local degees of freedom for trial and test functions
  integer :: indofPerElement

  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation

  ! Some triangulation arrays we need frequently
  integer, dimension(:,:), pointer :: p_IneighboursAtElement
  integer, dimension(:,:), pointer :: p_IelementsAtEdge
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  real(DP), dimension(:,:), pointer :: p_DvertexCoords

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! Underlying discretisation structure
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:), allocatable :: Kentry
  real(DP), dimension(:), allocatable :: Dentry

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_evalElementSet) :: revalElementSet

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable, target :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), pointer :: p_DcubPtsRef

  ! Cubature point weights
  real(DP), dimension(CUB_MAXCUBP) :: Domega

  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D

  ! number of cubature points on the reference element
  integer :: ncubp,icubp

  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder

  ! Whther the viscosity is constant or not
  logical :: bconstViscosity

    ! Currently we support only constant viscosity
    bconstViscosity = .true.

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rmatrixScalar%p_rspatialDiscrTest

    ! If a discretisation structure is present, take that one.
    if (present(rdiscretisation)) &
      p_rdiscretisation => rdiscretisation

    ! Get Kvert, Kadj, Kmid, Kmel,...
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
                                p_IneighboursAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IelementsAtEdge,&
                                p_IelementsAtEdge)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)

    NVT = p_rtriangulation%NVT

    ! Get KLD, KCol...
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar,p_Da)

    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! Assure that the element spaces are compatible
    if (elem_igetShape(p_relementDistribution%celement) .ne. &
        elem_igetShape(rmatrixScalar%p_rspatialDiscrTrial%RelementDistr(1)%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_uniDP')
      call sys_halt()
    end if

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference line [-1,1]
    call cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega)

    ! Allocate arrays accepting cubature point coordinates.
    ! We have ncubp cubature points on every egde.
    ! DcubPtsRef saves all possible combinations, which edge of
    ! one element might interact with one edge of one neighbour
    ! element.
    allocate(DcubPtsRefOnAllEdges(NDIM2D,ncubp,NVE))

    ! Put the 1D cubature points from to all of the edges on the
    ! reference element, so we can access them quickly.
    do iedge = 1,NVE
      call trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
      do i=1,ncubp
        DcubPtsRefOnAllEdges(1,i,iedge) = Dxi2D(i,1)
        DcubPtsRefOnAllEdges(2,i,iedge) = Dxi2D(i,2)
      end do
    end do

    ! We have always 2 elements in each patch...
    IELcount = 2

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,IELcount))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    allocate(Kentry(indofPerElement*2*indofPerElement*2))
    allocate(Dentry(indofPerElement*2*indofPerElement*2))

    ! Allocate memory for obtaining DOF`s:
    allocate(IdofsTempl(indofPerElement,2))

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    !
    ! We allocate space for 3 instead of 2 elements. The reason is that
    ! we later permute the values of the basis functions to get
    ! a local numbering on the patch. That is also the reason, we allocate
    ! not indofPerElement elements, but even indofPerElement,
    ! which is more than enough space to hold the values of the DOF`s of
    ! a whole element patch.

    allocate(Dbas(indofPerElement*2, &
            elem_getMaxDerivative(p_relementDistribution%celement),&
            ncubp,3))

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Do not calculate coordinates on the reference element -- we do this manually.
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    Bder(:) = .false.
    Bder(DER_FUNC) = .true.

      ! Fill the basis function arrays with 0. Essential, as only parts
      ! of them are overwritten later.
      Dbas = 0.0_DP

    ! We loop through all edges
    do IMT = 1,p_rtriangulation%NMT

      ! Check if we have 1 or 2 elements on the current edge
      if (p_IelementsAtEdge (2,IMT) .eq. 0) then

        ! This only happens if we have a boundary edge.
        ! The boundary handling is... doing nothing! So we can skip
        ! the handling of that edge completely!
        cycle

      end if

      ! On an example, we now show the relationship between all the DOF`s
      ! we have to consider in the current situation. We have two elements,
      ! let us say IEL1 and IEL2, with their local and global DOF`s
      ! (example for Q1~):
      !
      !   local DOF`s on the           corresponding global
      !   element and on the patch     DOF`s
      !
      !    +----4----+----7----+       +----20---+----50---+
      !    |    4    |    4    |       |         |         |
      !    |         |         |       |         |         |
      !    1 1     3 3 1     3 6       10        60        70
      !    |         |         |       |         |         |
      !    |    2    |    2    |       |         |         |
      !    +----2----+----3----+       +----40---+----30---+
      !        IEL1      IEL2              IEL1      IEL2
      !
      !
      ! On every element, we have 4 local DOF`s (1..4). On the other hand,
      ! we have "local DOF`s on the patch" (1..7), ehich we call "patch DOF`s"
      ! from now on. To every local DOF, there belongs a global DOF
      ! (10,20,30,... or whatever), which gives the coefficient of the basis
      ! function.
      !
      ! Numbering that way, the local DOF`s of IEL1 obviously coincide
      ! with the first couple of local DOF`s of the element patch! Only
      ! the local DOF`s of element IEL2 make trouble, as they have another
      ! numbering.

      ! Get the global DOF`s of the 1 or two elements
      call dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  IdofsTempl)

      ! Some of the DOF`s on element 2 may coincide with DOF`s on element 1.
      ! More precisely, some must coincide! Therefore, we now have to collect the
      ! DOF`s uniquely and to figure out, which local DOF`s of element 2
      ! must renumbered to the appropriate local patch-DOF (like in the
      ! example, where local DOF 1 of element IEL2 must be renumbered
      ! to local patch DOF 3!
      !
      ! As the first couple of local DOF`s of IEL1 coincide with the local
      ! DOF`s of the patch, we can simply copy them:

      ndof = indofPerElement
      Idofs(1:ndof) = IdofsTempl(1:ndof,1)

      ! Furthermore, we read IdofsTempl and store the DOF`s there in Idofs,
      ! skipping all DOF`s we already have and setting up the renumbering strategy
      ! of local DOF`s on element IEL2 to patch DOF`s.

      skiploop: do IDOFE = 1,indofPerElement

        ! Do we have the DOF?
        idof = IdofsTempl(IDOFE,IELcount)
        do JDOFE=1,ndof
          if (Idofs(JDOFE) .eq. idof) then
            ! Yes, we have it.
            ! That means, the local DOF idof has to be mapped to the
            ! local patch dof...
            IlocalDofRenum (IDOFE) = JDOFE

            ! Proceed with the next one.
            cycle skiploop
          end if
        end do

        ! We do not have that DOF! Append it to Idofs.
        ndof = ndof+1
        Idofs(ndof) = idof
        IlocalDofRenum (IDOFE) = ndof

        ! Save also the number of the local DOF.
        ! Note that the global DOF`s in IdofsTempl(1..indofPerElement)
        ! belong to the local DOF`s 1..indofPerElement -- in that order!
        IlocalDofs(ndof) = IDOFE

      end do skiploop

      ! Now we know: Our 'local' matrix (consisting of only these DOF`s we just
      ! calculated) is a ndofsTest*ndofsTrial matrix.
      !
      ! Now extract the corresponding entries from the matrix.
      ! Kentry is an index where the entries of the local matrix Dentry
      ! (which we build up later) can be found in the global matrix.
      !
      ! The test space gives us the rows; loop through them

      do IDOFE = 0,ndof-1

        ! The 'local' DOF IDOFE corresponds in the global matrix to line...
        irow = Idofs (1+IDOFE)

        ! Loop through that line to find the columns, indexed by the local
        ! DOF`s in the trial space.
        trialspaceloop: do JDOFE = 1,ndof

          do jcol = p_Kld(irow),p_Kld(irow+1)-1

            if (p_Kcol(jcol) .eq. Idofs(JDOFE)) then

              ! Found! Put as destination pointer to Kentry
              Kentry (IDOFE*ndof+JDOFE) = jcol

              ! Clear local matrix
              Dentry (IDOFE*ndof+JDOFE) = 0.0_DP

              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop

            end if

          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_uniDP')
          call sys_halt()

        end do trialspaceloop

      end do ! JDOFE

      ! Now we can set up the local matrix in Dentry. Later, we will plug it into
      ! the global matrix using the positions in Kentry.
      !
      ! The next step is to evaluate the basis functions in the cubature
      ! points on the edge. To compute the jump, this has to be done
      ! for both elements on the edge.
      !
      ! Figure out which edge on the current element is IMT.
      ! We need this local numbering later to determine on which edge we
      ! have to place the cubature points.
      do i = 1,IELcount

        IEL = p_IelementsAtEdge(i,IMT)
        do iedge = 1,NVE
          if (p_IedgesAtElement (iedge,IEL) .eq. IMT) exit
        end do

        ! Copy the coordinates of the corresponding cubature points
        ! to DcubPtsEval. We calculated the cubature points on the
        ! reference element in advance, so we can simply take them.

        p_DcubPtsRef(:,:,i) = DcubPtsRefOnAllEdges (:,:,iedge)

      end do

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementsAtEdge (1:IELcount,IMT), &
          ctrafoType, DpointsRef=p_DcubPtsRef, rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

      ! Apply the permutation of the local DOF`s on the test functions
      ! on element 2. The numbers of the local DOF`s on element 1
      ! coincides with the numbers of the local DOF`s on the patch.
      ! Those on element 2, we have to renumber according to the permutation
      ! so that they are in the correct order according to the DOF`s on the patch.
      !
      ! We copy the values of the basis functions to the space in
      ! p_DcubPtsTest which is reserved for the 3rd element!
      ! That way, Dbas has:
      !
      ! local DOF on element   :       1       2       3       4
      !
      ! Global DOF on element 1:      10      40      60      20
      ! Dbas(1..4,*,*,1):    d1psi10 d1psi40 d1psi60 d1psi20
      !
      ! Global DOF on elemenr 2:      60      30      70      50
      ! Dbas(1..4,*,*,2):    d2psi60 d2psi30 d2psi70 d2psi50
      !
      ! and will be transformed to:
      !
      ! local patch-DOF:            1       2       3       4       5      6        7
      ! global DOF:                10      40      60      20      30     70       50
      ! Dbas(1..7,*,*,1): d1psi10 d1psi40 d1psi60 d1psi20     0.0    0.0      0.0
      ! Dbas(1..7,*,*,2): --------------------------unused-----------------------
      ! Dbas(1..7,*,*,3):     0.0     0.0 d2psi60     0.0 d2psi30 d2psi70 d2psi50
      !
      ! ("d1psi10" = grad(psi_10) on element 1,
      !  "psi_10" = test function on global DOF #10)
      !
      ! That space Dbas(:,:,:,3) is unused up to now, initialise by 0.0 and
      ! partially overwrite with the reordered values of the basis functions.
      Dbas (:,:,:,3) = 0.0_DP
      Dbas (IlocalDofRenum(1:indofPerElement),:,:,3) &
        = Dbas (1:indofPerElement,:,:,2)

      ! What do we have now? When looking at the example, we have:
      !
      ! ndof = 7
      ! Idofs(1..7)             = global DOF`s in local numbering 1..7
      ! Dbas(1..7,*,1..ncubp,1) = values of basis functions on element 1
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF`s only appearing at element 2
      ! Dbas(1..7,*,1..ncubp,3) = values of basis functions on element 2
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF`s only appearing at element 1
      !
      ! Now we can start to integrate using this.

      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit interval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 &
            + (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2 )
      dedgeweight = dedgelength !* 0.5_DP

      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E gamma*nu/|E| int_E [u] [v] ds
      dcoeff = dgamma * dnu / dedgelength
!      dcoeff = dgamma * dnu * dedgelength
      
      ! Now we have the values of the basis functions in all the cubature
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      do IDOFE = 1,ndof

        ! Loop through the trial basis functions
        do JDOFE = 1,ndof

          dval = 0.0_DP

          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          do icubp = 1,ncubp

            ! [ phi ]   ( jump in the derivative of trial basis function)
            ! = [ phi  ,  phi ]
            dphi  = Dbas (JDOFE,DER_FUNC,icubp,1) &
                  - Dbas (JDOFE,DER_FUNC,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsi  = Dbas (IDOFE,DER_FUNC,icubp,1) &
                  - Dbas (IDOFE,DER_FUNC,ncubp-icubp+1,3)


            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dval = dval + Domega(icubp) * dedgeweight * (dphi*dpsi)

          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndof+JDOFE) = &
            Dentry ((IDOFE-1)*ndof+JDOFE) + dcoeff*dval

        end do

      end do

      ! Incorporate our "local" system matrix
      ! into the global matrix. The position of each entry DENTRY(X,Y)
      ! in the global matrix array A was saved in element Kentry(X,Y)
      ! before.
      ! Kentry gives the position of the additive contributions in Dentry.
      ! The entry is weighted by the current dtheta, which is usually
      ! the weighting parameter of the corresponding THETA-scheme of a
      ! nonstationary simulation. For stationary simulations, dtheta is typically
      ! 1.0 which includes the local matrix into the global one directly.)

      do IDOFE = 0,ndof-1
        do JDOFE = 1,ndof

          p_Da(Kentry(IDOFE*ndof+JDOFE)) = &
            p_Da(Kentry(IDOFE*ndof+JDOFE)) + &
            dtheta * Dentry (IDOFE*ndof+JDOFE)

        end do
      end do

      ! Proceed with next edge

    end do ! IMT

    ! Clean up allocated arrays and memory.
    deallocate(DcubPtsRefOnAllEdges)

    call elprep_releaseElementSet(revalElementSet)
    deallocate(p_DcubPtsRef)

    deallocate(Dbas)

    deallocate(Kentry)
    deallocate(Dentry)

    deallocate(IdofsTempl)

  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine ls_viscosity_model(Dii,rcollection,nu)
!<description>
  ! This subroutine calculates the kinematic viscosity based on different
  ! non-Newtonian fluid models, e.g. power law model.
!</description>

!<input>
  ! The second invariant of the deformation rate tensor, Dii(u)
  real(DP), intent(inout) :: Dii
  
  ! User defined collection structure
  type(t_collection), intent(in) :: rcollection
!</input>

!<output>
  ! Returns the value of the kinematic viscosity
  real(DP), intent(out) :: nu
!</output>  

!</subroutine>

  ! local variables
  ! Viscosity model constants
  real(DP) :: r, nu_0, nu_inf, landa
  
  ! Take care of the initial solutions
  ! where Dii(u) = 0
  if (Dii .eq. 0.0_DP) then
    Dii = 1.0_DP
  end if
  
  ! Switch on the viscosity model
  select case (rcollection%IquickAccess(1))
  
  case (1)
    ! Power Law model
    !  Viscosity defines with the POWER LAW:
    !                             r/2-1
    !   { \nu(Dii(u)) = \nu_0 Dii(u)^     
    !   { where  ( r>1 , \nu_0>0 ) 
    !
    !  and
    !
    !  Dii(u) = 1/2 ( 2D(u):2D(u) )  
    !  D(u) = 1/2 ( grad(u) + grad(u)^T ) 
    !
    ! Viscosity model constants
    nu_0 = rcollection%DquickAccess(17)
    r = rcollection%DquickAccess(18)
    
    ! Viscosity
    nu = nu_0 * Dii**(r/2.0_DP-1.0_DP)

  case (2)
    ! Carreau Law model
    !  Viscosity defines with the Carreau LAW:
    !                                                            r/2-1
    !   { \nu(Dii(u)) = \nu_inf + (\nu_0-\nu_inf) (1+\landa Dii(u))^     
    !   { where  ( \landa>0, r>1 , \nu_0>\nu_inf>=0 )
    !
    ! Viscosity model constants
    nu_0 = rcollection%DquickAccess(17)
    r = rcollection%DquickAccess(18)
    nu_inf = rcollection%DquickAccess(19)
    landa = rcollection%DquickAccess(20)    
    
    ! Viscosity    
    nu = nu_inf + (nu_0 - nu_inf) * (1 + landa*Dii)**(r/2.0_DP-1.0_DP)
  
  case default
    ! Not defined model!!
           
  end select
  
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine ls_viscosity_model_der(Dii,rcollection,dnu)
!<description>
  ! This subroutine calculates the derivative of the kinematic viscosity
  ! based on different non-Newtonian fluid models, e.g. power law model.
!</description>

!<input>
  ! The second invariant of the deformation rate tensor, Dii(u)
  real(DP), intent(inout) :: Dii
  
  ! User defined collection structure
  type(t_collection), intent(in) :: rcollection
!</input>

!<output>
  ! Returns the derivative of the kinematic viscosity
  real(DP), intent(out) :: dnu
!</output>  

!</subroutine>

  ! local variables
  ! Viscosity model constants
  real(DP) :: r, nu_0, nu_inf, landa

  ! Take care of the initial solutions
  ! where Dii(u) = 0
  if (Dii .eq. 0.0_DP) then
    Dii = 1.0_DP
  end if

  ! Switch on the viscosity model     
  select case (rcollection%IquickAccess(1))
  
  case (1)
    ! Power Law model
    !  Viscosity defines with the POWER LAW:
    !                             r/2-1
    !   { \nu(Dii(u)) = \nu_0 Dii(u)^     
    !   { where  ( r>1 , \nu_0>0 ) 
    !
    !  and
    !
    !  Dii(u) = 1/2 ( 2D(u):2D(u) )  
    !  D(u) = 1/2 ( grad(u) + grad(u)^T ) 
    !
    ! Viscosity model constants
    nu_0 = rcollection%DquickAccess(17)
    r = rcollection%DquickAccess(18)
    
    ! Derivative of the viscosity    
    dnu = nu_0 * (r/2.0_DP-1.0_DP) * Dii**(r/2.0_DP-2.0_DP)

  case (2)
    ! Carreau Law model
    !  Viscosity defines with the Carreau LAW:
    !                                                            r/2-1
    !   { \nu(Dii(u)) = \nu_inf + (\nu_0-\nu_inf) (1+\landa Dii(u))^     
    !   { where  ( \landa>0, r>1 , \nu_0>\nu_inf>=0 )
    !
    ! Viscosity model constants
    nu_0 = rcollection%DquickAccess(17)
    r = rcollection%DquickAccess(18)
    nu_inf = rcollection%DquickAccess(19)
    landa = rcollection%DquickAccess(20)    
    
    ! Derivative of the viscosity    
    dnu = (nu_0 - nu_inf) * (r/2.0_DP-1.0_DP) * landa * &
          (1 + landa*Dii)**(r/2.0_DP-2.0_DP)
  
  case default
    ! Not defined model!!
           
  end select
    
  end subroutine



  !****************************************************************************

!<subroutine>

  subroutine ls_nonlinear_weight(Dii,rcollection,gama)
!<description>
  ! This subroutine calculates the kinematic viscosity based on different
  ! non-Newtonian fluid models, e.g. power law model.
!</description>

!<input>
  ! The second invariant of the deformation rate tensor, Dii(u)
  real(DP), intent(inout) :: Dii
  
  ! User defined collection structure
  type(t_collection), intent(in) :: rcollection
!</input>

!<output>
  ! Returns the value of the kinematic viscosity
  real(DP), intent(out) :: gama
!</output>  

!</subroutine>

  ! local variables
  ! Viscosity model constants
  real(DP) :: r, nu_0, nu_inf, landa
  
  ! Take care of the initial solutions
  ! where Dii(u) = 0
  if (Dii .eq. 0.0_DP) then
    Dii = 1.0_DP
  end if
  
  ! Switch on the viscosity model
  select case (rcollection%IquickAccess(1))
  
  case (1)
    ! Power Law model
    !  Viscosity defines with the POWER LAW:
    !                             r/2-1
    !   { \nu(Dii(u)) = \nu_0 Dii(u)^     
    !   { where  ( r>1 , \nu_0>0 ) 
    !
    !  and
    !
    !  Dii(u) = 1/2 ( 2D(u):2D(u) )  
    !  D(u) = 1/2 ( grad(u) + grad(u)^T ) 
    !
    ! Viscosity model constants
    nu_0 = rcollection%DquickAccess(17)
    r = rcollection%DquickAccess(18)
    
    ! Viscosity
    gama = 1.0_DP/(nu_0 * Dii**(r/2.0_DP-1.0_DP))

  case (2)
    ! Carreau Law model
    !  Viscosity defines with the Carreau LAW:
    !                                                            r/2-1
    !   { \nu(Dii(u)) = \nu_inf + (\nu_0-\nu_inf) (1+\landa Dii(u))^     
    !   { where  ( \landa>0, r>1 , \nu_0>\nu_inf>=0 )
    !
    ! Viscosity model constants
    nu_0 = rcollection%DquickAccess(17)
    r = rcollection%DquickAccess(18)
    nu_inf = rcollection%DquickAccess(19)
    landa = rcollection%DquickAccess(20)    
    
    ! Viscosity    
    gama = 1.0_DP/(nu_inf + (nu_0 - nu_inf) * (1 + landa*Dii)**(r/2.0_DP-1.0_DP))
  
  case default
    ! Not defined model!!
           
  end select
  
  end subroutine

    
end

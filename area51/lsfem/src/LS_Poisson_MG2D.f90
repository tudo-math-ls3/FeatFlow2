!##############################################################################
!# ****************************************************************************
!# <name> LS_Poisson_MG2D </name>
!# ****************************************************************************
!# <purpose>   
!# This module solves the 2D Poisson eq. using LSFEM. The equation reads:
!#    { -\triangle p = f   in \Omega        
!#    {            p = 0   on \Gamma_D,    (1)
!#  where f is a body force.
!#
!# The second-order elliptic equation (1) is reformulated into first-order
!# system of equations using the first derivatives of the primary variable:
!#      flux: \bu = -\nabla(p).
!# The resulting system is the 'div-grad' system which reads:
!#    {        \nabla p = f   in \Omega        
!#    { \bu + \nabla(p) = 0   in \Omega          (2)
!#    {               p = 0   on \Gamma_D,       
!# The problem is solved in a coupled manner for the solution of:
!#   - primary variable  p
!#   - the fluxes        \bu.
!#
!# In an attempt to make the least-squares formulation "norm-equivalent",
!#  the div-grad system will be augmented by additional curl equation
!#       \nabla \times \bu = 0,
!#  which yields to the div-grad-curl system. Additional curl operator
!#  can be added to the boundary conitions, i.e.
!#      \bn \times \bu = 0  on \Gamma_D.
!#  So, the 'div-grad-curl' system reads
!#    {          \nabla p = f   in \Omega        
!#    {   \bu + \nabla(p) = 0   in \Omega
!#    { \nabla \times \bu = 0   in \Omega         (3)
!#    {    \bn \times \bu = 0   on \Omega_D
!#    {                 p = 0   on \Gamma_D,     
!#
!# The LSFEM formulation then applied to system of equations (2) and (3)
!# which yiedls a symmetric, positive-definite linear system of equations.
!# This routine uses the Multigrid as linear solver. The discretisation uses
!# the block assembly method to evaluate the mattrix all-in-one.
!# </purpose>
!#
!# Author:    Masoud Nickaeen
!# First Version: Jan  14, 2011
!# Last Update:   Jun  25, 2013
!# 
!##############################################################################

module LS_Poisson_MG2D

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
  use pprocerror
  use meshmodification
  use paramlist
  
  use LS_Poisson_callback

  use blockmatassemblybase
  use blockmatassembly
  use collection
  use feevaluation2
  use feevaluation
  use vectorio
  use matrixmodification
  use dofmapping
  use meshmodification
  use matrixio
  use statistics
  
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

    ! A variable describing the discrete boundary conditions.    
    type(t_discreteBC) :: rdiscreteBC

    ! Cubature information structure which defines the cubature formula.
    type(t_scalarCubatureInfo) :: rcubatureInfo
  
  end type

!</types>


contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine LS_Poisson_mg
  
!<description>
  ! This is an all-in-one stokes solver for directly solving a stokes
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.)  Read in parametrisation
  ! 2.)  Read in triangulation
  ! 3.)  Set up the discretization structure
  ! 4.)  Initialization of BCs, solution/RHS vectors
  ! 5.)  Set up matrix
  ! 6.)  Set up RHS
  ! 7.)  Implement BCs
  ! 8.)  Solve the problem
  ! 9.)  Post processing
  ! 10.) Release all variables, finish.
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
  type(t_vectorBlock) :: rvector,rrhs
  
  ! Max. Min. level to solve the problem
  integer :: NLMAX, NLMIN
  
  ! Loop index
  integer :: i
  
  ! Extra curl equation and boundary condition
  integer :: icurl,icurl_boundary
    
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
  ! 1)- Read in parametrisation
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call ls_initParams(rparams,NLMAX,NLMIN,icurl_boundary,rcollection)
  

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 2)- Read in triangulation
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_grid(Rlevels,rboundary,rparams,rcollection,NLMAX,NLMIN)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 3)- Set up the discretization structure
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_structure(Rlevels,rboundary,rparams)

  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 4)- Initialization of Boundary conditions, and the solution/RHS vectors.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize the discrete boundary conditions
  do i = NLMIN, NLMAX
    call ls_BCs_Dirichlet_one(Rlevels(i)%rdiscretisation,rboundary,&
        Rlevels(i)%rDiscreteBC,icurl_boundary,rcollection)
  end do
  
  ! Initialize the RHS and fine level solution vector
  call ls_Init_RhsndSolution(Rlevels,rrhs,rvector,NLMAX)
  
              
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 5)- System Matrix Assembly
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_MatAssembly(Rlevels,rcollection,rparams,NLMIN,NLMAX)
    

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 6)- Right hand side Vector Assembly
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    
  call ls_RHS_Assembly(rrhs,rcollection,Rlevels(NLMAX)%rcubatureInfo) 
    
    
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 7)- Implement Boundary Conditions
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_BCs_Dirichlet(Rlevels,rrhs,rvector,NLMAX,NLMIN,rcollection)

     
  ! Start Solver Timer
  call stat_startTimer(rtimerSolver)   
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 8)- Solve the System
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_Solver_linear(Rlevels,rrhs,rvector,rparams)
  ! Stop Solver Timer
  call stat_stopTimer(rtimerSolver)


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 9)- Post-processing starts here. Writing the solution to file ...
  ! export GMV/VTK files
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_postprocess(rboundary,Rlevels(NLMAX)%rmatrix,&
       rvector,Rlevels(NLMAX)%rtriangulation,Rlevels(NLMAX)%rcubatureInfo,&
       Rlevels(NLMAX)%rdiscretisation,rparams,icurl_boundary,rcollection)
  
  
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! 10)- Clean up, free all the memory which is used for variables.
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call ls_cleanup(rvector,rrhs,Rlevels,rboundary,NLMAX,NLMIN,rparams)
  
  
  ! Stop the total-time Timer
  call stat_stopTimer(rtimerTotal)
  ! Print some statistics
  call output_lbrk ()
  call output_line ("Time for Total Solution: "//&
                    trim(sys_sdL(rtimerTotal%delapsedReal,10)))  
    
  call output_lbrk ()
  call output_line ("Time for Linear Solver: "//&
                    trim(sys_sdL(rtimerSolver%delapsedReal,10)))
    
  end subroutine


  !****************************************************************************
  
!<subroutine>
  subroutine ls_initParams(rparams,NLMAX,NLMIN,icurl_boundary,&
                           rcollection)

 !<description>  
  ! In this subroutine, the collection structrure is initialized.
 !</description>                

 !<output>
  ! All parameters in LSFEM solver
  type(t_parlist), intent(out) :: rparams
  
  ! Collection structure for callback routines  
  type(t_collection), intent(out) :: rcollection
  
  ! Multigrid levels
  integer, intent(out) :: NLMAX,NLMIN
  
  ! Adding extra curl boundary condition 
  integer, intent(out) :: icurl_boundary  
 !</output>

 !</subroutine>

  ! Local variables
  integer :: icurl
  
  ! reading data from the *.dat file 
  ! We want to solve our Poisson problem on level...
  call parlst_init(rparams)
  call parlst_readfromfile (rparams, "./data/ls_poisson.dat") 

  ! The multigrid coarse and fine levels
  call parlst_getvalue_int (rparams, 'MESH', 'NLMAX', NLMAX, 6)
  call parlst_getvalue_int (rparams, 'MESH', 'NLMIN', NLMIN, 4)

  ! Adding the extra curl equation to the formulation
  call parlst_getvalue_int (rparams, 'GENERAL', 'icurl', icurl, 1)    
  
  ! For later use in matrix assembly routines
  ! save the information in the rcollection
  if (icurl .eq. 1) then
    rcollection%DquickAccess(1) = 1.0_DP
  else
    rcollection%DquickAccess(1) = 0.0_DP
  end if
  
  ! Adding the extra curl equation to the formulation
  call parlst_getvalue_int (rparams, 'GENERAL', 'icurl_boundary', icurl_boundary, 1) 

  ! The problem to solve, by default sqare domain
  call parlst_getvalue_int (rparams, 'GENERAL', 'Problem', &
       rcollection%IquickAccess(1), 0) 
  
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
  ! Set up a discretisation and cubature info structure.
  ! Also, this subroutine tells the code which finite element to use.
  ! Finally, it creates a block matrix structure.
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
  integer(I32) :: Velm, Pelm

  ! Type of cubature rule for numerical integration
  integer(I32) :: ccubType
  
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
  ! a block discretisation structure that specifies the blocks in the
  ! solution vector. In this problem, we have three blocks.
  ! Do this for all levels
  do i = NLMIN, NLMAX
    call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation,3,&
                               Rlevels(i)%rtriangulation, rboundary)
  end do                      

  ! rdiscretisation%RspatialDiscr is a list of scalar 
  ! discretisation structures for every component of the solution vector.
  ! We have a solution vector with three components:
  !  Component 1 =  u   1st component of flux variable
  !  Component 2 =  v   2nd component of flux variable
  !  Component 3 =  p   the primary variable
  
  ! Read the finite element for fluxe variables
  call parlst_getvalue_string (rparams, 'MESH', 'Velm', sstring)
  Velm = elem_igetID(sstring)

  ! Here we set up one discretisation structure for the 
  ! 1st component of the vector variable
  do i = NLMIN, NLMAX
    call spdiscr_initDiscr_simple (&
      Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
        Velm, Rlevels(i)%rtriangulation, rboundary)
  end do
          
  ! ...and copy this structure also to the discretisation structure
  ! of the 2nd component (v). This needs no additional memory,
  ! as both structures will share the same dynamic information afterwards.
  do i = NLMIN, NLMAX
    call spdiscr_duplicateDiscrSc(&
         Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
           Rlevels(i)%rdiscretisation%RspatialDiscr(2))
  end do
  
  ! For the primary variable, we set up a separate discretisation
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
    ! The global system looks like this, a full 3*3 block matrix
    ! which is symmetric and positive definite.
    !
    !  ( A11 A12 A13 )
    !  ( A21 A22 A23 )
    !  ( A31 A32 A33 )
    !  
    ! Create the matrix structure of the u. Block A11,
    call bilf_createMatrixStructure (&
    Rlevels(i)%rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
    Rlevels(i)%rmatrix%RmatrixBlock(1,1))  
                     
    ! Block A12,
    call bilf_createMatrixStructure (&
        Rlevels(i)%rdiscretisation%RspatialDiscr(1), LSYSSC_MATRIX9, &
         Rlevels(i)%rmatrix%RmatrixBlock(1,2))  
    ! Block A13
    call bilf_createMatrixStructure (&
        Rlevels(i)%rdiscretisation%RspatialDiscr(3), LSYSSC_MATRIX9, &
        Rlevels(i)%rmatrix%RmatrixBlock(1,3),&
        Rlevels(i)%rdiscretisation%RspatialDiscr(1))                                    
                                    
    ! Use u structure for the v component. Block A22,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
      Rlevels(i)%rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
    ! Block A23,
    call lsyssc_duplicateMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,3),&
      Rlevels(i)%rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)

    ! Create the matrix structure of the primary variable p. Block A33,
    call bilf_createMatrixStructure (&
       Rlevels(i)%rdiscretisation%RspatialDiscr(3), LSYSSC_MATRIX9, &
        Rlevels(i)%rmatrix%RmatrixBlock(3,3))
    
      
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
                              icurl_boundary,rcollection)
                
 !<description>  
  ! In this subroutine we discretise the boundary conditions and
  ! prepare them to be applied to the matrix/RHS/sol later.
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
  
  ! The curl boundary condition
  integer, intent(in) :: icurl_boundary
 !</input>

!</subroutine>

  ! Local variables
  integer :: i
  type(t_boundaryRegion) :: rboundaryRegion 

  ! Create a t_discreteBC structure where we store all discretised boundary
  ! conditions.
  call bcasm_initDiscreteBC(rdiscreteBC)

  ! Switch on the problem we want to solve
  select case (rcollection%IquickAccess(1))
  
  case (0)
    
    ! The unit square problem
    ! -----------------------
    ! We first set up the boundary conditions for u
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
    !   We specify icomponent='3' to indicate that we set up the
    !   Dirichlet BC`s for the third (here: one and only) component in the 
    !   solution vector.
    ! - Discretise the boundary condition so that the BC`s can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,3,&
                                     rboundaryRegion,rdiscreteBC,&
                                     getBoundaryValues_2D)
                           
    ! Edge 2 of boundary component 1.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,3,&
                                     rboundaryRegion,rdiscreteBC,&
                                     getBoundaryValues_2D)
                           
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,3,&
                                     rboundaryRegion,rdiscreteBC,&
                                     getBoundaryValues_2D)
  
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,3,&
                                     rboundaryRegion,rdiscreteBC,&
                                     getBoundaryValues_2D)
			      			       				       

    if (icurl_boundary .eq. 1) then
      !  Edge 1 of boundary component 1. That is it.
  	  call boundary_createRegion(rboundary,1,1,rboundaryRegion)
  	  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
  			       
      ! Edge 3 of boundary component 1. That is it.
  	  call boundary_createRegion(rboundary,1,3,rboundaryRegion)
  	  call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)				       
  			
  			
      ! Edge 2 of boundary component 2. That is it.
  	  call boundary_createRegion(rboundary,1,2,rboundaryRegion)
  	  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
  			       
      ! Edge 4 of boundary component 2. That is it.
  	  call boundary_createRegion(rboundary,1,4,rboundaryRegion)
  	  call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)				
			       
    endif
    
  case default
    ! Un-known problem
    call output_line ("Unknown problem.", OU_CLASS_WARNING, OU_MODE_STD, &
            "ls_BCs_Dirichlet_One")
    call sys_halt()   
  end select

  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_Init_RhsndSolution(Rlevels,rrhs,rvector,NLMAX)
                
 !<description>  
  ! Initializing the RHS and the solution vector
 !</description>                

 !<output>
  ! Block vectors
  type(t_vectorBlock), intent(out) :: rvector,rrhs
 !</output>
 
 !<input>
  ! Max. level to solve the problem
  integer, intent(in) :: NLMAX
 !</input>

 !<inputoutput>
  ! An array of problem levels for the multigrid solver
  type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>

!</subroutine>

  ! local variables
  
  ! Let's start
  ! Create a RHS and a solution vector based on the discretisation.
  ! Fill with zero.
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                                rrhs,.true.)
  call lsysbl_createVectorBlock (Rlevels(NLMAX)%rdiscretisation,&
                               rvector,.true.)
  
  end subroutine

  !****************************************************************************

!<subroutine>
  subroutine ls_MatAssembly(Rlevels,rcollection,rparams,NLMIN,NLMAX)
                
 !<description>  
  ! Initializing the solution vectors on all levels by calculating the 
  ! memory required to the interlevel projections.
 !</description>                
 
 !<input>
  ! Level info.
  integer, intent(in) :: NLMAX,NLMIN
  
  ! Collection structure for callback routines  
  type(t_collection) :: rcollection
  
  ! All parameters in LSFEM solver
  type(t_parlist), intent(in) :: rparams  
 !</input>
  
 !<inputoutput>
  ! An array of problem levels for the multigrid solver
  type(t_level), dimension(:), pointer :: Rlevels
 !</inputoutput>

!</subroutine>

  ! Local variables
  integer :: i
    
  ! Temporary block matrix and vectors
  type(t_matrixBlock), pointer :: p_rmatrix

  ! On all levels, we have to set up the system matrix  
  nullify(p_rmatrix)

  do i=NLMAX,NLMIN,-1
  
    ! Get the matrix on the current level.
    p_rmatrix => Rlevels(i)%rmatrix

    ! Assemble the whole system matrix on each level  
    call bma_buildMatrix (p_rmatrix,BMA_CALC_STANDARD,ls_poisson_Matrix,&
       rcubatureInfo=Rlevels(i)%rcubatureInfo,rcollection=rcollection)

  end do


  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_RHS_Assembly(rrhs,rcollection,rcubatureInfo)
                
 !<description>  
  ! Initializing the solution vectors on all levels by calculating the 
  ! memory required to the interlevel projections.
 !</description>                
 
 !<input>
  ! Cubature information structure which defines the cubature formula.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo 
 !</input> 
 
 !<inputoutput>
  ! RHS block vector
  type(t_vectorBlock), intent(inout) :: rrhs
  
  ! Collection structure for callback routines  
  type(t_collection) :: rcollection  
 !</inputoutput>
  
!</subroutine>


  ! Local variables  
  
  ! Call the subroutine to assemble the RHS 
  call bma_buildVector (rrhs,BMA_CALC_STANDARD,ls_poisson_rhs,&
                        rcubatureInfo=rcubatureInfo,rcollection=rcollection)
                 
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_BCs_Dirichlet(Rlevels,rrhs,rvector,NLMAX,NLMIN,rcollection)
                
 !<description>  
  ! Implementing BCs to the matrix and solution/RHS vectors.
 !</description>                

 !<input>
  integer, intent(in) :: NLMAX,NLMIN
  
  ! Collection structure for callback routines  
  type(t_collection), intent(in) :: rcollection
 !</input>
 
 !<inputoutput>
  ! Block vectors 
  type(t_vectorBlock), intent(inout) :: rvector,rrhs
   
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
  
  ! Implement the filter  
  call vecfil_discreteBCrhs (rrhs)
  call vecfil_discreteBCsol (rvector)

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
    case (3,4)
    ! The Jacobi-type coarse grid solvers!!
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
    call parlst_getvalue_double (rparams, 'MULTI', 'DampSmoothMG', DampSmoothMG, 1.0_DP)  
    call linsol_initSSOR (p_rsmoother, DampSmoothMG, .true.)
    case (4)
    call parlst_getvalue_double (rparams, 'MULTI', 'DampSmoothMG', DampSmoothMG, 1.0_DP)  
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
    call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)
    
    ! Check for the linear solver divergence 
    if (p_rsolverNode%DFinalDefect .gt. 1E8) then
      call output_lbrk()
      call output_line ('Linear Solver is diverged :(')
      call output_line ('Try to modify linear solver properties.')
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
      call output_line ('Linear Solver is diverged :(')
      call output_line (&
       'Try to modify the linear solver properties.')
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
  subroutine ls_postprocess(rboundary,rmatrix,rvector,rtriangulation,&
              rcubatureInfo,rdiscretisation,rparams,icurl_boundary,rcollection)
                
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
  
  ! The curl boundary condition
  integer, intent(in) :: icurl_boundary
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
  !   to export GMV/VTK outputs
  !   to write the real/projected data
  integer ::  ExporType
  integer :: Vtild, Ptild
 
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
    
  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation) :: rblockDiscr  

  ! Norm Calculations
  integer :: L2U,L2P,H1U,H1P
  real(DP) :: dintvalue, dintvalue1
  type(t_fev2Vectors) :: revalVectors
  

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

  if ( (Vtild .eq. 1) .or. (Ptild .eq. 1) ) then
      
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
                              icurl_boundary,rcollection)

    ! Hang the pointer into the vector.
    rprjVector%p_rdiscreteBC => rprjDiscreteBC

    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector) 
    
    if (ExporType .eq. 0) then
    
      ! Start UCD export to VTK file:
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/poisson_lsfem.vtk')

      ! Write primary variable p
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

      ! Write flux variables u,v
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'flux',p_Ddata,p_Ddata2)
    
    else
      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/poisson_lsfem.gmv')  

      ! Write primary variable p
      call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
      
      ! Write flux variables u,v
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'flux',p_Ddata,p_Ddata2)
        
    
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
        trim(sucddir)//'/poisson_lsfem.vtk')

      ! Write primary variable p
      call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)

      ! Write flux variables u,v
      call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'flux',p_Ddata,p_Ddata2)
    
    else
      ! Start UCD export to GMV file:
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//'/poisson_lsfem.gmv')  

      ! Write primary variable p
      call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'p',UCD_VAR_STANDARD,p_Ddata)
      
      ! Write flux variables u,v    
      call lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
      call lsyssc_getbase_double (rvector%RvectorBlock(2),p_Ddata2)
      call ucd_addVarVertBasedVec(rexport,'flux',p_Ddata,p_Ddata2)    
    
    end if
  
  end if ! end of real or projected data condition
    
  ! Write the file to disc, that is it.
  call ucd_write (rexport)
  call ucd_release (rexport)

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate The L^2 and H^1 norms for analytical solutions
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  call parlst_getvalue_int (rparams, 'POST', 'L2U', L2U, 0)
  call parlst_getvalue_int (rparams, 'POST', 'L2P', L2P, 0)
  call parlst_getvalue_int (rparams, 'POST', 'H1U', H1U, 0)
  call parlst_getvalue_int (rparams, 'POST', 'H1P', H1P, 0)
  
  if (L2U == 1) then
    
    ! L^2 Norm flux
    ! Add x-component vector
    rcollection%IquickAccess(7) = 0
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(1),0)   ! u
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)    
    
    ! Add y-component vector
    rcollection%IquickAccess(7) = 1      
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(2),0)   ! v
    call bma_buildIntegral (dintvalue1,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)
    
    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error flux variables')
    call output_line ('------------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)))  
    
    call output_line ('L2flux:'//&
    trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)), coutputMode=OU_MODE_BENCHLOG)

  end if

  if (H1U == 1) then    
    ! H^1 Norm flux    
    ! Add x-component vector
    rcollection%IquickAccess(7) = 0
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(1),1)   ! u,x,y
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_H1_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)    
    
    ! Add y-component vector
    rcollection%IquickAccess(7) = 1      
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(2),1)   ! v,x,y
    call bma_buildIntegral (dintvalue1,BMA_CALC_STANDARD,&
    ls_H1_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    call fev2_releaseVectorList(revalVectors)
    
    ! Print the Norm value
    call output_lbrk()
    call output_line ('H^1 Error flux variables')
    call output_line ('------------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)))  

    call output_line ('H1flux:'//&
    trim(sys_sdEP(sqrt(dintvalue+dintvalue1),15,6)), coutputMode=OU_MODE_BENCHLOG)

  end if
  
  if (L2P == 1) then
    ! L^2 Norm primary variable
    ! Add p
    rcollection%IquickAccess(7) = 2
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(3),0)   ! p
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_L2_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('L^2 Error primary variable')
    call output_line ('--------------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('L2primary:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG)

  end if


  if (H1P == 1) then
    ! H^1 Norm primary variable
    ! Add p
    rcollection%IquickAccess(7) = 2
    call fev2_addVectorToEvalList(revalVectors,&
       rvector%RvectorBlock(3),1)   ! p,x,y
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
    ls_H1_Norm,rcollection=rcollection, &
    revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)

    ! Print the Norm value
    call output_lbrk()
    call output_line ('H^1 Error primary variable')
    call output_line ('--------------------------')
    call output_line (trim(sys_sdEP(sqrt(dintvalue),15,6)))  
    call fev2_releaseVectorList(revalVectors)

    call output_line ('H1primary:'//&
    trim(sys_sdEP(sqrt(dintvalue),15,6)), coutputMode=OU_MODE_BENCHLOG)
    
  end if
     
  end subroutine


  !****************************************************************************

!<subroutine>
  subroutine ls_cleanup(rvector,rrhs,Rlevels,rboundary,NLMAX,NLMIN,rparams)
 
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
  type(t_vectorBlock) :: rvector,rrhs
  !</inputoutput>

!</subroutine>

  ! Local variables
  ! Loop index, fictitious BC
  integer :: i, ifictitious
  
  ! Now, clean up so that all the memory is available again.  
  ! Release the block matrix/vectors
  call lsysbl_releaseVector (rvector)
  call lsysbl_releaseVector (rrhs)
  
  do i = NLMAX, NLMIN, -1
    call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
  end do

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
  subroutine ls_poisson_Matrix(RmatrixData,rassemblyData,rmatrixAssembly,&
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
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA12
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA13
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA21,p_DlocalMatrixA22
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA23
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
  real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA33 
  
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA33,p_DbasTestA33
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13

  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA11,p_rmatrixDataA33
  type(t_bmaMatrixData), pointer :: p_rmatrixDataA13
  
  ! Additional curl equation parameter
  real(DP) :: Dcurl
  
  ! If (the additional curl equation is augmented) then
  !   Dcurl = 1.0
  ! else 
  !  Dcurl = 1.0
  ! end if
  Dcurl = rcollection%DquickAccess(1)
  
  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rmatrixDataA11 => RmatrixData(1,1)
  p_rmatrixDataA33 => RmatrixData(3,3)
  p_rmatrixDataA13 => RmatrixData(1,3)
  
  p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
  p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
  p_DbasTrialA33 => RmatrixData(3,3)%p_DbasTrial
  p_DbasTestA33 => RmatrixData(3,3)%p_DbasTest
  p_DbasTrialA13 => RmatrixData(1,3)%p_DbasTrial
  p_DbasTestA13 => RmatrixData(1,3)%p_DbasTest  
  
  
  p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
  p_DlocalMatrixA12 => RmatrixData(1,2)%p_Dentry
  p_DlocalMatrixA13 => RmatrixData(1,3)%p_Dentry

  p_DlocalMatrixA21 => RmatrixData(2,1)%p_Dentry
  p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
  p_DlocalMatrixA23 => RmatrixData(2,3)%p_Dentry

  p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
  p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
  p_DlocalMatrixA33 => RmatrixData(3,3)%p_Dentry

  
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
      dval = p_DcubWeight(icubp,iel) * (  dbasJx*dbasIx + dbasJ*dbasI + &
                                          Dcurl*dbasJy*dbasIy  )
      p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
     
      ! A22
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIy + dbasJ*dbasI + &
                                          Dcurl*dbasJx*dbasIx  )      
      p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
         
      ! A12      
      dval = p_DcubWeight(icubp,iel) * (  dbasJy*dbasIx - Dcurl*dbasJx*dbasIy  )       
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
      dval = p_DcubWeight(icubp,iel) *(  dbasJx*dbasI  )
      p_DlocalMatrixA13(jdofe,idofe,iel) = p_DlocalMatrixA13(jdofe,idofe,iel) + dval
      
      ! A31
      p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + dval
           
      ! A23
      dval = p_DcubWeight(icubp,iel) *(  dbasJy*dbasI  )
      p_DlocalMatrixA23(jdofe,idofe,iel) = p_DlocalMatrixA23(jdofe,idofe,iel) + dval
      
      ! A32
      p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + dval
                      
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
      dbasIx = p_DbasTestA33(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTestA33(idofe,DER_DERIV2D_Y,icubp,iel)
      
      ! Inner loop over the DOF's j=1..ndof, which corresponds to
      ! the basis function Phi_j:
      do jdofe=1,p_rmatrixDataA33%ndofTrial
      
      ! Fetch the contributions of the (trial) basis function Phi_j
      dbasJx = p_DbasTrialA33(jdofe,DER_DERIV2D_X,icubp,iel)
      dbasJy = p_DbasTrialA33(jdofe,DER_DERIV2D_Y,icubp,iel)

      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local matrices.
      ! A33  
      dval = p_DcubWeight(icubp,iel) *(dbasJx*dbasIx + dbasJy*dbasIy)
      p_DlocalMatrixA33(jdofe,idofe,iel) = p_DlocalMatrixA33(jdofe,idofe,iel) + dval
      
      end do ! idofe
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel  
  
  
  end subroutine

  !****************************************************************************


!<subroutine>
  subroutine ls_poisson_rhs(rvectorData,rassemblyData,rvectorAssembly,&
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
  real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
  real(DP), dimension(:,:), pointer :: p_DlocalVector3
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTest1,p_DbasTest3
  real(DP), dimension(:,:), pointer :: p_DcubWeight
  type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData3
  real(DP), dimension(:,:,:), pointer :: p_Dpoints
  real(DP) :: dfx, dfy, dx, dy


  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  p_rvectorData1 => RvectorData(1)
  p_rvectorData3 => RvectorData(3)

  p_DlocalVector1 => RvectorData(1)%p_Dentry
  p_DlocalVector2 => RvectorData(2)%p_Dentry
  p_DlocalVector3 => RvectorData(3)%p_Dentry

  p_DbasTest1 => RvectorData(1)%p_DbasTest
  p_DbasTest3 => RvectorData(3)%p_DbasTest 
    
  ! Get the real coordinates of the cubature points
  p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
 
     
  ! Calculate the RHS of the velocities
  ! Loop over the elements in the current set.
  do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

    ! Get the coordinates of the cubature point.
    dx = p_Dpoints(1,icubp,iel)
    dy = p_Dpoints(2,icubp,iel)
    
    ! Outer loop over the DOF's i=1..ndof on our current element,
    ! which corresponds to the (test) basis functions Phi_i:
    do idofe=1,p_rvectorData1%ndofTest
    
      ! Fetch the contributions of the (test) basis functions Phi_i
      ! into dbasI
      dbasIx = p_DbasTest1(idofe,DER_DERIV2D_X,icubp,iel)
      dbasIy = p_DbasTest1(idofe,DER_DERIV2D_Y,icubp,iel)
                
      ! Values of the RHS functions for the u and v components
      dfx = 32.0_DP*(dy*(1.0_DP - dy) + dx*(1.0_DP - dx))
      dfy = dfx
      
      dval = dfx*dbasIx
      ! Multiply the values of the basis functions by
      ! the cubature weight and sum up into the local vectors.
      p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
        p_DcubWeight(icubp,iel) * dval
        
      dval = dfx*dbasIy
      p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
        p_DcubWeight(icubp,iel) * dval
      
    end do ! jdofe

    end do ! icubp
  
  end do ! iel
  
  
  ! Calculate the RHS of the primary variable
  ! there is no rhs for the third component. So, leave it zero!
  
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
  real(DP) :: dval1,dval2,dx,dy, pi
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

  ! Get cubature weights data
  p_DcubWeight => rassemblyData%p_DcubWeight
  
  ! Get the coordinates of the cubature points
  p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

  pi = 3.1415926535897932_DP
  dintvalue = 0.0_DP

  select case (rcollection%IquickAccess(7))
  
  case (0)
    !X-component (u)
    ! Loop over the elements in the current set.
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)   
     
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the analytic function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = -16.0_DP*dy*(1.0_DP - dy)*(1.0_DP -2.0_DP*dx)
       
      ! Get the error of the FEM function to the analytic function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel
    
  case (1)
    ! Y-Velocity (v)
    ! Loop over the elements in the current set.
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)   
        
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = -16.0_DP*dx*(1.0_DP - dx)*(1.0_DP -2.0_DP*dy)
            
      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
      ! Multiply the values by the cubature weight and sum up
      ! into the (squared) L2 error:
      dintvalue = dintvalue + &
        p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
      
    end do ! icubp
    
    end do ! iel
  
  case (2)
    ! Primary variable
    ! Loop over the elements in the current set.
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)   
        
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the value of the bubble function
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dval1 = 16.0_DP*dx*dy*(1.0_DP - dx)*(1.0_DP - dy)
      
      ! Get the error of the FEM function to the bubble function
      dval2 = p_Dfunc(icubp,iel)
      
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
  real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy
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
    !X-component (u)  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = 32.0_DP*dy*(1.0_DP - dy)
      dderivY1 = -16.0_DP*(1.0_DP - 2.0_DP*dx)*(1.0_DP - 2.0_DP*dy)
      
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
    ! Y-component (v)
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = -16.0_DP*(1.0_DP - 2.0_DP*dx)*(1.0_DP - 2.0_DP*dy)
      dderivY1 = 32.0_DP*dx*(1.0_DP - dx)

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
    ! Primary   
    ! Loop over the elements in the current set.
    do iel = 1,nelements

    ! Loop over all cubature points on the current element
    do icubp = 1,npointsPerElement

      ! Get the derivatives of the bubble function in the cubature point
      dx = p_Dpoints(1,icubp,iel)
      dy = p_Dpoints(2,icubp,iel)
      
      dderivX1 = 16.0_DP*dy*(1.0_DP - dy)*(1.0_DP -2.0_DP*dx)
      dderivY1 = 16.0_DP*dx*(1.0_DP - dx)*(1.0_DP -2.0_DP*dy)

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


end 

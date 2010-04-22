!#########################################################################################
!# ***************************************************************************************
!# <name> elasticity_2d_disp_smallDeform_static </name>
!# ***************************************************************************************
!#
!# <purpose>
!#   This module solves the basic elasticity problem:
!#     - 2D
!#     - pure displacement formulation
!#     - small deformation (i.e., a linear problem)
!#     - static
!#   There is also the possibility to compute a 2D Poisson problem.
!# </purpose>
!#########################################################################################

module elasticity_2d_disp_smallDef_stat

  implicit none

contains

! ****************************************************************************************

!<subroutine>
  subroutine elast_2d_disp_smallDef_stat

    use storage
    use genoutput
    use paramlist
    use derivatives
    use boundary
    use bilinearformevaluation
    use linearformevaluation
    use bcassembly
    use linearsolver
    use linearsolverautoinitialise
    use filtersupport
    use linearsystemscalar
    use matrixfilters
    use vectorfilters
    use scalarpde
    use ucd
    use collection
    use elasticity_basic
    use elasticity_callback
  
!<description>
  ! This routine realises the basic elasticity solver. It performs the following tasks:
  !
  ! 1.) read in parametrisation
  ! 2.) read in triangulation
  ! 3.) set up RHS
  ! 4.) set up matrix
  ! 5.) create solver structure
  ! 6.) solve the problem
  ! 7.) write solution to GMV file
  ! 8.) release all variables, finish
!</description>

!</subroutine>

    ! array of problem levels for the multigrid solver
    type(t_level), dimension(:), pointer :: Rlevels

    ! object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! bilinear form (matrix) and linear form (RHS) describing the problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! block vectors for solution, right hand side and temp. vector
    type(t_vectorBlock) :: rsol, rrhs, rtempBlock

    ! variable for selecting a specifig boundary region
    type(t_boundaryRegion) :: rboundaryRegion

    ! solver node
    type(t_linsolNode), pointer :: p_rsolver

    ! array for the block structured system matrix
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

    ! filter chain that describes how to filter the matrix/vector before/during the
    ! solution process (used for implementing Dirichlet boundary conditions)
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! error indicator during initialisation of the solver
    integer :: ierror
    
    ! output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata, p_Ddata2
    real(DP), Dimension(:,:), pointer :: p_DvertexCoords

    ! some auxiliary variables
    integer :: i, j, k, ilev, irow, jcol, nlevels

    ! collection structure to provide additional information to the coefficient routine.
    type(t_collection) :: rcollection

    ! parameter list read from the parameter file
    type (t_parlist) :: rparams

    
    ! +------------------------------------------------------------------------
    ! | READ PARAMETER FILE
    ! +------------------------------------------------------------------------

    ! call subroutine that reads in the parameter file
    call elast_readParameterFile(rprob, rparams)

    ! +------------------------------------------------------------------------
    ! | BOUNDARY AND TRIANGULATION
    ! +------------------------------------------------------------------------
  
    ! set number of levels
    nlevels = rprob%ilevelMax - rprob%ilevelMin + 1
  
    ! allocate memory for all levels
    allocate(Rlevels(rprob%ilevelMin:rprob%ilevelMax))

    call output_line('reading boundary parameterisation from file ' // &
                     trim(rprob%sgridFilePRM) // '...')
    ! read in the parameterisation of the boundary and save it to rboundary.
    call boundary_read_prm(rboundary, rprob%sgridFilePRM)
        
    ! read in the basic triangulation.
    call output_line('reading triangulation from file ' // &
                     trim(rprob%sgridFileTRI) // '...')
    call tria_readTriFile2D(Rlevels(rprob%ilevelMin)%rtriangulation, &
                            rprob%sgridFileTRI, rboundary)
     
    ! refine it once
    call tria_quickRefine2LevelOrdering(rprob%ilevelMin-1, &
           Rlevels(rprob%ilevelMin)%rtriangulation, rboundary)
    
    ! create information about adjacencies etc.
    call tria_initStandardMeshFromRaw(Rlevels(rprob%ilevelMin)%rtriangulation, rboundary)

    ! create all refinement levels
    do ilev = rprob%ilevelMin+1, rprob%ilevelMax

      ! refine the grid using the 2-Level-Ordering algorithm
      call tria_refine2LevelOrdering(Rlevels(ilev-1)%rtriangulation, &
             Rlevels(ilev)%rtriangulation,rboundary)
      
      ! Create a standard mesh
      call tria_initStandardMeshFromRaw(Rlevels(ilev)%rtriangulation, rboundary)
    
    end do

    ! +------------------------------------------------------------------------
    ! | DISCRETISATION
    ! +------------------------------------------------------------------------

    ! set up a block discretisation structure that specifies the blocks in the
    ! solution vector for all levels. For the scalar Poisson problem, we only have
    ! one block. Do this for all levels
    do ilev = rprob%ilevelMin, rprob%ilevelMax
      call spdiscr_initBlockDiscr(Rlevels(ilev)%rdiscretisation, rprob%nblocks, &
                                  Rlevels(ilev)%rtriangulation, rboundary)
    end do

    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation structures
    ! for every component of the solution vector.
    do ilev = rprob%ilevelMin, rprob%ilevelMax
      ! create a discretisation for the first component
      call spdiscr_initDiscr_simple(Rlevels(ilev)%rdiscretisation%RspatialDiscr(1),&
          rprob%celement, rprob%ccubature2D, Rlevels(ilev)%rtriangulation, rboundary)

      ! ...and copy this structure to the discretisation structure of the 2nd component
      ! (y-displacement) (no additional memory needed)
      do j = 2, rprob%nblocks
        call spdiscr_duplicateDiscrSc (&
            Rlevels(ilev)%rdiscretisation%RspatialDiscr(1),&
            Rlevels(ilev)%rdiscretisation%RspatialDiscr(j))
      enddo
    end do

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    do ilev = rprob%ilevelMin, rprob%ilevelMax

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr(Rlevels(ilev)%rdiscretisation, &
                                        Rlevels(ilev)%rmatrix)
      
      ! in case of the Poisson equation only one block has to be set up
      if (rprob%cequation .eq. EQ_POISSON) then 
  
        ! generate the structure of the system matrix
        call bilf_createMatrixStructure(Rlevels(ilev)%rdiscretisation%RspatialDiscr(1), &
               LSYSSC_MATRIX9, Rlevels(ilev)%rmatrix%RmatrixBlock(1,1))
        
        ! bilinear form (grad Psi_j, grad Phi_i)
        rform%itermCount = 2
        rform%Idescriptors(1,1) = DER_DERIV_X
        rform%Idescriptors(2,1) = DER_DERIV_X
        rform%Idescriptors(1,2) = DER_DERIV_Y
        rform%Idescriptors(2,2) = DER_DERIV_Y
  
        ! use constant coefficients
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.
        rform%Dcoefficients(1)  = 1.0
        rform%Dcoefficients(2)  = 1.0
  
        ! Build the matrix entries. The callback function elast_mat_Poisson_2D for the
        ! coefficients is only used if specifying ballCoeffConstant = .FALSE. and
        ! BconstantCoeff = .FALSE.
        call bilf_buildMatrixScalar(rform, .true., &
               Rlevels(ilev)%rmatrix%RmatrixBlock(1,1), elast_mat_Poisson_2D)
      else if (rprob%cequation .eq. EQ_ELASTICITY) then
        
        ! common information for all blocks
        rform%itermCount = 2
        rform%ballCoeffConstant = .true.
        rform%BconstantCoeff = .true.

        do irow = 1, rprob%nblocks
          do jcol = 1, rprob%nblocks
            call bilf_createMatrixStructure(&
                   Rlevels(ilev)%rdiscretisation%RspatialDiscr(irow), &
                   LSYSSC_MATRIX9, Rlevels(ilev)%rmatrix%RmatrixBlock(irow,jcol))
            if (irow .eq. 1 .and. jcol .eq. 1) then
              ! block (1,1)
              rform%Idescriptors(1,1) = DER_DERIV_X
              rform%Idescriptors(2,1) = DER_DERIV_X
              rform%Idescriptors(1,2) = DER_DERIV_Y
              rform%Idescriptors(2,2) = DER_DERIV_Y
          
              rform%Dcoefficients(1)  = 2*rprob%dmu + rprob%dlambda
              rform%Dcoefficients(2)  = rprob%dmu

            else if (irow .eq. 1 .and. jcol .eq. 2) then
              ! block (1,2)
              rform%Idescriptors(1,1) = DER_DERIV_Y
              rform%Idescriptors(2,1) = DER_DERIV_X
              rform%Idescriptors(1,2) = DER_DERIV_X
              rform%Idescriptors(2,2) = DER_DERIV_Y

              rform%Dcoefficients(1)  = rprob%dlambda
              rform%Dcoefficients(2)  = rprob%dmu

            else if (irow .eq. 2 .and. jcol .eq. 1) then
              ! block (2,1)
              rform%Idescriptors(1,1) = DER_DERIV_X
              rform%Idescriptors(2,1) = DER_DERIV_Y
              rform%Idescriptors(1,2) = DER_DERIV_Y
              rform%Idescriptors(2,2) = DER_DERIV_X
          
! BRAL: sicher?
              rform%Dcoefficients(1)  = rprob%dlambda
              rform%Dcoefficients(2)  = rprob%dmu

            else if (irow .eq. 2 .and. jcol .eq. 2) then
              ! block (2,2)
              rform%Idescriptors(1,1) = DER_DERIV_X
              rform%Idescriptors(2,1) = DER_DERIV_X
              rform%Idescriptors(1,2) = DER_DERIV_Y
              rform%Idescriptors(2,2) = DER_DERIV_Y

              rform%Dcoefficients(1)  = rprob%dmu
              rform%Dcoefficients(2)  = 2*rprob%dmu + rprob%dlambda
            endif

! BRAL:
! As soon as I have to use the callback function, the following call has to be put
! inside the above if-else-block.
            call bilf_buildMatrixScalar(rform, .true., &
                   Rlevels(ilev)%rmatrix%RmatrixBlock(irow,jcol), elast_mat_Poisson_2D)
          enddo
        enddo

      endif
    end do

    ! set up structure for RHS and solution vector by using the block matrix as template
    call lsysbl_createVecBlockIndMat(Rlevels(rprob%ilevelMax)%rmatrix, rrhs, .false.)
    call lsysbl_createVecBlockIndMat(Rlevels(rprob%ilevelMax)%rmatrix, rsol, .false.)

    ! clear RHS and solution vector on the finest level
    call lsysbl_clearVector(rrhs)
    call lsysbl_clearVector(rsol)

    ! assemble RHS vector corresponding to linear form (f,Phi_j)
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! build the RHS entries on the finest level using callback functions to compute
    ! the coefficients
    if (rprob%cequation .eq. EQ_POISSON) then 
      ! Poisson equation
      
      ! use the callback function elast_RHS_Poisson_2D_vol
      call linf_buildVectorScalar(Rlevels(rprob%ilevelMax)%rdiscretisation%RspatialDiscr(1),&
                    rlinform, .true., rrhs%RvectorBlock(1), elast_RHS_Poisson_2D_vol)

    else if (rprob%cequation .eq. EQ_ELASTICITY) then
      ! elasticity equation
      
      ! compute volumen forces using the callback routine elast_RHS_2D_vol
      ! (x-direction: rcollection%IquickAccess(1) = 1, 
      !  y-direction: rcollection%IquickAccess(1) = 2)
      do irow = 1, rprob%nblocks
        rcollection%IquickAccess(1) = irow
        call linf_buildVectorScalar(&
               Rlevels(rprob%ilevelMax)%rdiscretisation%RspatialDiscr(irow), &
               rlinform, .true., rrhs%RvectorBlock(irow), elast_RHS_2D_vol, rcollection)
      enddo
    endif

    ! print number of DOF
    call lsysbl_getbase_double(rsol, p_Ddata)
    call output_line('Number of DOF: ' // trim(sys_siL(size(p_Ddata),12)) )
  
    ! For implementing boundary conditions, we use a filter technique with discretised
    ! boundary conditions. This means, we first have to calculate a discrete version of
    ! the analytic BC, which we can implement into the solution/RHS vectors using the
    ! corresponding filter.

    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    do ilev = rprob%ilevelMin, rprob%ilevelMax
      call bcasm_initDiscreteBC(Rlevels(ilev)%rdiscreteBC)
    end do

    ! set up the boundary conditions per boundary and segment
    do i = 1, rprob%nboundaries
      do j = 1,rprob%NboundarySegments(i)
        ! Create 'boundary region',  which is simply a part of the boundary corresponding
        ! to a boundary segment. A boundary region roughly contains the type, the min/max
        ! parameter value and whether the endpoints are inside the region or not.
        call boundary_createRegion(rboundary,i,j,rboundaryRegion)
        ! mark start and end point as belonging to the region
        rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
        do k = 1, rprob%nblocks
          if (rprob%Cbc(k,j,i) .eq. BC_DIRICHLET) then
            do ilev = rprob%ilevelMin, rprob%ilevelMax
              ! We use this boundary region to specify that we want to have Dirichlet
              ! boundary there. The following call does the following:
              ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
              !   We specify icomponent='k' to indicate that we set up the Dirichlet BCs
              !   for the k-th component in the solution vector.
              ! - Discretise the boundary condition so that the BCs can be applied
              !   to matrices and vectors.
              ! - Add the calculated discrete BCs to rdiscreteBC for later use.
              call bcasm_newDirichletBConRealBD(Rlevels(ilev)%rdiscretisation, k, &
                     rboundaryRegion, Rlevels(ilev)%rdiscreteBC, elast_boundValue_2D)
            end do  
          else if (rprob%Cbc(k,j,i) .eq. BC_NEUMANN) then
            ! store the current segment number in rcollection%IquickAccess(2) to make it
            ! accessible in the callback routine
            rcollection%IquickAccess(2) = j
            ! The non-zero Neumann contributations on the current boundary region are
            ! added to the RHS vector which already contains the volumetric contributions.
            if (rprob%cequation .eq. EQ_POISSON) then 
              call linf_buildVectorScalarBdr2d(rlinform, rprob%ccubature1D, .false., &
                     rrhs%RvectorBlock(1), elast_RHS_Poisson_2D_bound, rboundaryRegion, &
                     rcollection)
            else if (rprob%cequation .eq. EQ_ELASTICITY) then
              ! store the current component in rcollection%IquickAccess(1), such that
              ! the callback routine elast_RHS_2D_surf knows whether to compute surface
              ! forces in x- or y-direction
              rcollection%IquickAccess(1) = k
              call linf_buildVectorScalarBdr2d(rlinform, rprob%ccubature1D, .false., &
                rrhs%RvectorBlock(k), elast_RHS_2D_surf, rboundaryRegion, rcollection)
            endif
          else
            call output_line('Invalid BC found!', OU_CLASS_ERROR, OU_MODE_STD, &
                             'elast_2d_disp_smallDef_stat')
            call sys_halt()
          end if
        enddo
      end do ! end segments
    end do ! end boundaries

    ! attach the boundary conditions to the matrix (on all levels)
    do ilev = rprob%ilevelMin, rprob%ilevelMax
      Rlevels(ilev)%rmatrix%p_rdiscreteBC => Rlevels(ilev)%rdiscreteBC
    end do
!BRAL: testen, ob x- und y-Dirichlet-BC wirklich unabhaengig sind!!! (Block)
  
    ! attach the boundary conditions to the solution and the RHS vector
    rrhs%p_rdiscreteBC => Rlevels(rprob%ilevelMax)%rdiscreteBC
    rsol%p_rdiscreteBC => Rlevels(rprob%ilevelMax)%rdiscreteBC
    
    ! Next step is to implement boundary conditions into the matrix, the solution and the
    ! RHS vector. This is done using a vector/matrix filter for discrete boundary
    ! conditions that modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs(rrhs)
    call vecfil_discreteBCsol(rsol)
    do ilev = rprob%ilevelMin, rprob%ilevelMax
      call matfil_discreteBC(Rlevels(ilev)%rmatrix)
    end do

    ! +------------------------------------------------------------------------
    ! | SOLVER
    ! +------------------------------------------------------------------------
      
    ! During the linear solve, the boundary conditions must frequently be imposed to
    ! vectors. This is done using a 'filter chain'. (As the linear solver does not work
    ! with the actual solution vectors but with defect vectors instead, a filter for
    ! implementing the real boundary conditions would be wrong.)
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! The above filter chain is attached to the solver, so that it automatically filters
    ! the vector during the solution process. Set the following pointer for this.
    p_RfilterChain => RfilterChain

    call linsolinit_initFromFile(p_rsolver, rparams, "SOLVER", nlevels, p_RfilterChain)

    ! Attach the system matrices to the solver.

    ! Copy all matrices to a big matrix array and transfer it to the setMatrices routines,
    ! which intitialise the matrices on all levels. Note that no new memory is allocated,
    ! Rmatrices(:) only contains 'links' to existing matrices.
    allocate(Rmatrices(rprob%ilevelMin:rprob%ilevelMax))
    do ilev = rprob%ilevelMin, rprob%ilevelMax
      call lsysbl_duplicateMatrix(Rlevels(ilev)%rmatrix, Rmatrices(ilev), &
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_rsolver, Rmatrices(rprob%ilevelMin:rprob%ilevelMax))

    ! Rmatrices can be released immediately
    do ilev = rprob%ilevelMin,rprob%ilevelMax
      call lsysbl_releaseMatrix(Rmatrices(ilev))
    end do
    deallocate(Rmatrices)
    
    ! initialise solver structure/data
    call linsol_initStructure(p_rsolver, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line('Error in initialisation of the solver structure: ' // &
                       sys_siL(ierror,8), OU_CLASS_ERROR, OU_MODE_STD, &
                       'elast_2d_disp_smallDef_stat')
      call sys_halt()
    endif
    call linsol_initData(p_rsolver, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line('Error in initialisation of the solver data: ' // &
                       sys_siL(ierror,8), OU_CLASS_ERROR, OU_MODE_STD, &
                       'elast_2d_disp_smallDef_stat')
      call sys_halt()
    endif

    ! Solve the system. We want to solve Ax=b with b being the 'real' RHS and x being
    ! the 'real' solution vector, so linsol_solveAdaptively is the appropriate routine.
    ! If b is a defect and x a defect correction to be added to a solution vector, we
    ! would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively(p_rsolver, rsol, rrhs, rtempBlock)

    call output_line('*********************************************************')
    ! number of iterations
    call output_line('Number of iterations: ' // sys_siL(p_rsolver%iiterations,8) )
    ! rate of convergence
    call output_line('Convergence rate: ' // sys_sdL(p_rsolver%dconvergenceRate,5) )
    ! rate of asymptotic convergence
    call output_line('Asymptotic convergence rate: ' // &
                      sys_sdL(p_rsolver%dasymptoticConvergenceRate,5) )
    call output_line('*********************************************************')
    call output_lbrk()

    ! calculate and display errors
    call elast_calcErrors(rprob, rsol)

    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) then
      sucddir = './gmv'
    endif

    ! get pointers to the arrays containing the displacement solutions
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_Ddata)
    if (rprob%cequation .eq. EQ_ELASTICITY) then
      call lsyssc_getbase_double(rsol%RvectorBlock(2),p_Ddata2)
    end if
  
    ! add displacements to the vertex coordinates when deformation is to be displayed
    ! in the gmv file
    if (rprob%cshowDeformation .eq. YES .and. rprob%cequation .eq. EQ_ELASTICITY) then
      ! get pointer to the arrays containing the vertex coordinates
      call storage_getbase_double2D( &
             Rlevels(rprob%ilevelMax)%rtriangulation%h_DvertexCoords, p_DvertexCoords)

      ! add displacements
      do i = 1,Rlevels(rprob%ilevelMax)%rtriangulation%NVT
        p_Dvertexcoords(1,i) = p_Dvertexcoords(1,i) + p_Ddata(i)
        p_Dvertexcoords(2,i) = p_dvertexCoords(2,i) + p_Ddata2(i)
      end do
    end if

    ! initialise UCD export to GMV file
    call ucd_startGMV(rexport, UCD_FLAG_STANDARD, &
                      Rlevels(rprob%ilevelMax)%rtriangulation, &
                      trim(sucddir) // '/elasticity_2d.gmv')

    ! add the displacement solution to the file
    if (rprob%cequation .eq. EQ_ELASTICITY) then
      call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata, p_Ddata2)
    else if (rprob%cequation .eq. EQ_POISSON) then
      call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata)
    end if
!BRAL: output of derivatives, strains and stresses is missing...
  
    ! write the file
    call ucd_write(rexport)
    call ucd_release(rexport)

    ! -----------------------
    ! simulation finished
    ! -----------------------
    
    ! clean up so that all the memory is available again.

    ! release solver data and structure and the solver node itself
    call linsol_doneData(p_rsolver)
    call linsol_doneStructure(p_rsolver)
    call linsol_releaseSolver(p_rsolver)

    ! release block matrix/vectors
    call lsysbl_releaseVector(rsol)
    call lsysbl_releaseVector(rtempBlock)
    call lsysbl_releaseVector(rrhs)
    do ilev = rprob%ilevelMax, rprob%ilevelMin, -1
      call lsysbl_releaseMatrix(Rlevels(ilev)%rmatrix)
    end do

    ! release discrete version of the boundary conditions
    do ilev = rprob%ilevelMax, rprob%ilevelMin, -1
      call bcasm_releaseDiscreteBC(Rlevels(ilev)%rdiscreteBC)
    end do

    ! release discretisation structure
    do ilev = rprob%ilevelMax, rprob%ilevelMin, -1
      call spdiscr_releaseBlockDiscr(Rlevels(ilev)%rdiscretisation)
    end do
    
    ! release triangulation
    do ilev = rprob%ilevelMax, rprob%ilevelMin, -1
      call tria_done(Rlevels(ilev)%rtriangulation)
    end do
    deallocate(Rlevels)
    
    ! release domain
    call boundary_release(rboundary)

    ! release all manually allocated arrays
    deallocate(rprob%NboundarySegments)
    deallocate(rprob%Cbc)
    if (rprob%csimulation .eq. SIMUL_REAL) then
      deallocate(rprob%DforceSurface)
    endif
    if (rprob%nevalPoints .gt. 0) then
      deallocate(rprob%DevalPoints, rprob%Dvalues)
      if (rprob%nrefSols .gt. 0) then
        deallocate(rprob%DrefSols)
      endif
    endif

  end subroutine elast_2d_disp_smallDef_stat

end module elasticity_2d_disp_smallDef_stat


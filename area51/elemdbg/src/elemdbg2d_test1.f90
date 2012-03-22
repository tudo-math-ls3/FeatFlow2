!##############################################################################
!# ****************************************************************************
!# <name> elemdbg2d_test1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module elemdbg2d_test1

  use fsystem
  use genoutput
  use paramlist
  use storage
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use bcassemblybase
  use linearalgebra
  use dofmapping
  use basicgeometry
  use meshregion
  use triangulation
  use element
  use spatialdiscretisation
  use transformation
  use elementpreprocessing
  use linearsystemscalar
  use linearsystemblock
  use discretebc
  use scalarpde
  use pprocerror
  use stdoperators
  use meshmodification
  use ucd
  use spdiscprojection
  use convection
  use collection
  use sortstrategy
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use pprocerror
    
  use elemdbg2d_callback
  
  use disto2d_aux
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine elemdbg2d_1(rparam, sConfigSection, itest)
  type(t_parlist), intent(INOUT) :: rparam
  character(LEN=*), intent(IN) :: sConfigSection
  integer, intent(IN) :: itest
  
!<description>
!</description>

!</subroutine>

  type(t_boundary) :: rboundary
  type(t_triangulation) :: rtriangulation
  type(t_blockDiscretisation) :: rdiscretisation
  type(t_linearForm) :: rlinform
  type(t_matrixBlock) :: rmatrix
  type(t_vectorBlock), target :: rvecSol, rvecRhs, rvectmp
  type(t_boundaryRegion) :: rboundaryRegion
  type(t_discreteBC), target :: rdiscreteBC
  type(t_linsolNode), pointer :: p_rsolver, p_rprecond
  type(t_matrixBlock), dimension(1) :: Rmatrices
  type(t_errorScVec) :: rerror
  integer :: NLMIN,NLMAX,ierror,ilvl
  real(DP), dimension(:,:), allocatable, target :: Derror
  integer, dimension(:,:), allocatable :: Istat
  integer :: isolver, ioutput, nmaxiter,ccubature,idistType,idistLevel,idistLvl
  integer(I32) :: celement, cshape
  real(DP) :: ddist, depsRel, depsAbs, drelax, daux1, daux2, ddist2
  character(LEN=64) :: selement,scubature
  type(t_bilinearForm) :: rform
  integer :: iwritemesh, h_Ipermute
  type(t_ucdexport) :: rexport
  real(DP), dimension(:), pointer :: p_Ddata
  real(DP) :: dnu,dbeta1,dbeta2,dupsam,dgamma
  integer :: istabil, isolution, ifillin, clocalh
  type(t_convStreamlineDiffusion) :: rconfigSD
  type(t_jumpStabilisation) :: rconfigEOJ
  type(t_collection) :: rcollect
  integer, dimension(:), allocatable :: p_Iedges, p_Itemp
  integer :: iedge,ibc,nedgecount
  type(t_vectorBlock) :: rvelocityVector
  type(t_blockDiscretisation) :: rdiscretisationVel
  
    h_Ipermute = ST_NOHANDLE

    ! Fetch minimum and maximum levels
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMIN', NLMIN, -1)
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMAX', NLMAX, -1)
    
    if((NLMIN .lt. 1) .or. (NLMAX .lt. NLMIN)) then
      call output_line('Invalid NLMIN/NLMAX parameters',&
           OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end if

    ! Fetch mesh distortion type
    call parlst_getvalue_int(rparam, sConfigSection, 'IDISTTYPE', idistType, 0)

    ! Fetch mesh distortion level
    call parlst_getvalue_int(rparam, sConfigSection, 'IDISTLEVEL', idistLevel, 0)

    ! Fetch mesh distortion parameters
    call parlst_getvalue_double(rparam, sConfigSection, 'DMESHDISTORT', ddist, 0.1_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DMESHDISTORT2', ddist2, 0.0_DP)
    
    ! Fetch element and cubature rule
    call parlst_getvalue_string(rparam, sConfigSection, 'SELEMENT', selement, '')
    call parlst_getvalue_string(rparam, sConfigSection, 'SCUBATURE', scubature, '')
    
    ! Fetch desired solution
    call parlst_getvalue_int(rparam, sConfigSection, 'ISOLUTION', isolution, 0)
    
    ! Fetch solver type
    call parlst_getvalue_int(rparam, sConfigSection, 'ISOLVER', isolver, 0)
    
    ! Fetch solver output level
    call parlst_getvalue_int(rparam, sConfigSection, 'IOUTPUT', ioutput, 0)
    
    ! Fetch maximum number of iterations
    call parlst_getvalue_int(rparam, sConfigSection, 'NMAXITER', nmaxiter, 100)
    
    ! Fetch solver tolerances
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSABS', depsAbs, 1E-11_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSREL', depsRel, 1E-8_DP)

    ! Fetch problem parameters
    call parlst_getvalue_double(rparam, sConfigSection, 'DNU', dnu, 1.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DBETA1', dbeta1, 0.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DBETA2', dbeta2, 0.0_DP)
    
    ! Fetch stabilisation parameters
    call parlst_getvalue_int(rparam, sConfigSection, 'ISTABIL', istabil, 0)
    call parlst_getvalue_int(rparam, sConfigSection, 'CLOCALH', clocalh, 0)
    call parlst_getvalue_double(rparam, sConfigSection, 'DUPSAM', dupsam, 1.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DGAMMA', dgamma, 0.01_DP)
    
    ! Fetch relaxation parameter for CG-SSOR
    call parlst_getvalue_double(rparam, sConfigSection, 'DRELAX', drelax, 1.2_DP)

    ! Fetch fill-in level for ILU(k) preconditioner
    call parlst_getvalue_int(rparam, sConfigSection, 'IFILLIN', ifillin, 0)

    ! Writing of the mesh
    call parlst_getvalue_int(rparam, sConfigSection, 'IWRITEMESH', iwritemesh, 0)
    
    ! Parse element and cubature
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)
    
    ! Get the shape of the element
    cshape = elem_igetShape(celement)
    if(cshape .ne. cub_igetShape(ccubature)) then
      call output_line('Element and cubature formula incompatible', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end if
    
    ! Allocate arrays
    allocate(Derror(2,NLMIN:NLMAX))
    allocate(Istat(5,NLMIN:NLMAX))
    
    call output_separator(OU_SEP_STAR)
    call output_line('ELEMENT-DEBUGGER: 2D TEST #1')
    call output_line('============================')
    
    ! Print out that we are going to do:
    select case(itest)
    case(201)
      call output_line('System.............: L2-projection')
    case(202)
      call output_line('System.............: Poisson')
    case(203)
      call output_line('System.............: Convection-Diffusion')
      call output_line('DNU................: ' // trim(sys_sdEP(dnu,20,12)))
      call output_line('DBETA1.............: ' // trim(sys_sdEP(dbeta1,20,12)))
      call output_line('DBETA2.............: ' // trim(sys_sdEP(dbeta2,20,12)))
      select case(istabil)
      case(0)
        call output_line('Stabilisation......: none')
      case(1)
        call output_line('Stabilisation......: Streamline-Diffusion')
        call output_line('DUPSAM.............: ' // trim(sys_sdEP(dupsam,20,12)))
      case(2)
        call output_line('Stabilisation......: Jump-Stabilisation')
        call output_line('DGAMMA.............: ' // trim(sys_sdEP(dgamma,20,12)))
      case(3)
        call output_line('Stabilisation......: Streamline-Diffusion (kernel)')
        call output_line('DUPSAM.............: ' // trim(sys_sdEP(dupsam,20,12)))
        call output_line('CLOCALH............: ' // trim(sys_si(clocalh,4)))
      case default
        call output_line('Invalid ISTABIL parameter', &
          OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
        call sys_halt()
      end select
    case default
      call output_line('Invalid ITEST parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    select case(isolution)
    case(0)
      call output_line('Solution...........: u(x,y) = sin(pi*x) * sin(pi*y)')
    case(1)
      call output_line('Solution...........: u(x,y) = x')
    case(2)
      call output_line('Solution...........: u(x,y) = y')
    case(3)
      call output_line('Solution...........: u(x,y) = 16*x*(1-x)*y*(1-y)')
    case default
      call output_line('Invalid ISOLUTION parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    select case(cshape)
    case(BGEOM_SHAPE_TRIA)
      call output_line('Coarse Mesh........: TRIA.tri')
    case(BGEOM_SHAPE_QUAD)
      call output_line('Coarse Mesh........: QUAD.tri')
    case default
      call output_line('Element is not a valid 2D element!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    call output_line('NLMIN..............: ' // trim(sys_siL(NLMIN,4)))
    call output_line('NLMAX..............: ' // trim(sys_siL(NLMAX,4)))
    select case(idistType)
    case(0)
      ! no distortion => nothing else to print
    case(1)
      call output_line('Distortion Type....: index-based')
    case(2)
      call output_line('Distortion Type....: X-line-wise index-based')
    case(3)
      call output_line('Distortion Type....: Y-line-wise index-based')
    case default
      call output_line('Invalid IDISTTYPE parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    if(idistType .ne. 0) then
      call output_line('Distortion Level...: ' // trim(sys_siL(idistLevel,4)))
      call output_line('Mesh Distortion....: ' // trim(sys_sdL(ddist,8)))
    end if
    if((idistType .ge. 2) .and. (idistType .le. 3)) then
      call output_line('Mesh Distortion 2..: ' // trim(sys_sdL(ddist2,8)))
    end if
    call output_line('Element............: ' // trim(selement))
    call output_line('Cubature rule......: ' // trim(scubature))
    select case(isolver)
    case(0)
      call output_line('Solver.............: UMFPACK4')
    case(1)
      call output_line('Solver.............: CG-SSOR')
      call output_line('Relaxation.........: ' // trim(sys_sdL(drelax,8)))
    case(2)
      call output_line('Solver.............: BiCGStab-ILU(k)')
      if(ifillin .lt. 0) then
        call output_line('IFILLIN must not be less than 0', &
          OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
        call sys_halt()
      end if
      call output_line('Allowed Fill-In....: ' // trim(sys_siL(ifillin,4)))
    case default
      call output_line('Invalid ISOLVER parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    call output_line('Output Level.......: ' // trim(sys_siL(ioutput,4)))
    call output_line('Maximum Iter.......: ' // trim(sys_siL(nmaxiter,12)))
    call output_line('Absolute EPS.......: ' // trim(sys_sdEP(depsAbs,20,12)))
    call output_line('Relative EPS.......: ' // trim(sys_sdEP(depsRel,20,12)))
    call output_lbrk()
    
    ! Copy important parameters into quick-access arrays of the collection,
    ! these ones are needed by the callback functions.
    rcollect%IquickAccess(1) = itest
    rcollect%IquickAccess(2) = isolution
    rcollect%DquickAccess(1) = dnu
    rcollect%DquickAccess(2) = dbeta1
    rcollect%DquickAccess(3) = dbeta2

    ! Read in parametrisation
    select case(cshape)
    case(BGEOM_SHAPE_TRIA)
      call boundary_read_prm(rboundary, './pre/TRIA.prm')
    case(BGEOM_SHAPE_QUAD)
      call boundary_read_prm(rboundary, './pre/QUAD.prm')
    end select
    
    ! Loop over all levels
    do ilvl = NLMIN, NLMAX
    
      call output_line('Processing Level ' // trim(sys_siL(ilvl,4)) // '...')

      ! Now read in the basic triangulation.
      select case(cshape)
      case(BGEOM_SHAPE_TRIA)
        call tria_readTriFile2D (rtriangulation, './pre/TRIA.tri', rboundary)
      case(BGEOM_SHAPE_QUAD)
        call tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
      end select
       
      if(idistLevel .le. 0) then
        idistLvl = max(1,ilvl-idistLevel)
      else
        idistLvl = idistLevel
      end if
       
      if((idistLvl .le. ilvl) .and. (idistType .ne. 0)) then
      
        ! Refine up to distortion level
        call tria_quickRefine2LevelOrdering (idistLvl-1, rtriangulation,rboundary)
        
        ! Distort the mesh
        select case(idistType)
        case (1)
          ! index-based distortion
          call meshmod_disturbMesh (rtriangulation, ddist)
        
        case (2)
          ! X-line-wise index-based distortion
          call disto2d_distortQuadLineX(rtriangulation, ddist, ddist2)
        
        case (3)
          ! Y-line-wise index-based distortion
          call disto2d_distortQuadLineY(rtriangulation, ddist, ddist2)
        
        end select
        
        ! Refine up to current level
        if(idistLvl .lt. ilvl) then
          call tria_quickRefine2LevelOrdering (ilvl-idistLvl, rtriangulation,rboundary)
        end if
      
      else
      
        ! Refine up to current level
        call tria_quickRefine2LevelOrdering (ilvl-1, rtriangulation,rboundary)
        
      end if
      
      ! And create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

      ! Set up discretisation
      call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          celement, ccubature, rtriangulation, rboundary)
                   
      ! Create matrix structure
      call lsysbl_createMatBlockByDiscr(rdiscretisation, rmatrix)
      if ((itest .eq. 203) .and. (istabil .eq. 2)) then
        ! Extended matrix stencil
        call bilf_createMatrixStructure(rdiscretisation%RspatialDiscr(1),&
                                 LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1),&
                                 cconstrType=BILF_MATC_EDGEBASED)
      else
        call bilf_createMatrixStructure(rdiscretisation%RspatialDiscr(1),&
                                 LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1))
      end if
            
      ! Create vectors
      call lsysbl_createVecBlockIndMat(rmatrix, rvecSol, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecRhs, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecTmp, .true.)
      
      ! Set up discrete BCs
      if(itest .ne. 201) then
        call bcasm_initDiscreteBC(rdiscreteBC)
        call boundary_createRegion(rboundary,1,1,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
             rboundaryRegion,rdiscreteBC,getBoundaryValues2D,rcollect)
        call boundary_createRegion(rboundary,1,2,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
             rboundaryRegion,rdiscreteBC,getBoundaryValues2D,rcollect)
        call boundary_createRegion(rboundary,1,3,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
             rboundaryRegion,rdiscreteBC,getBoundaryValues2D,rcollect)
        call boundary_createRegion(rboundary,1,4,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
             rboundaryRegion,rdiscreteBC,getBoundaryValues2D,rcollect)
                               
        ! Assign BCs
        rmatrix%p_rdiscreteBC => rdiscreteBC
        rvecSol%p_rdiscreteBC => rdiscreteBC
        rvecRhs%p_rdiscreteBC => rdiscreteBC
      end if

      ! Store statistics
      Istat(1,ilvl) = rmatrix%NEQ
      Istat(2,ilvl) = rmatrix%rmatrixBlock(1,1)%NA
      Istat(3,ilvl) = rtriangulation%NVT
      Istat(4,ilvl) = rtriangulation%NMT
      Istat(5,ilvl) = rtriangulation%NEL

      ! Assemble system
      select case(itest)
      case(201)
        ! Assemble Mass matrix
        call stdop_assembleSimpleMatrix(rmatrix%RmatrixBlock(1,1), &
                                        DER_FUNC2D, DER_FUNC2D)
      
      case(202)
        ! Assemble Laplace matrix
        call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1))

      case(203)
        ! Assemble Laplace matrix
        call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1),.true.,dnu)
        
        ! Stabilisation?
        select case (istabil)
        case (0)
          ! Assemble convective operator, sum up to the Laplace operator
          rform%itermcount = 2
          rform%Idescriptors(1,1) = DER_DERIV_X
          rform%Idescriptors(2,1) = DER_FUNC
          rform%Idescriptors(1,2) = DER_DERIV_Y
          rform%Idescriptors(2,2) = DER_FUNC
          rform%Dcoefficients(1) = dbeta1
          rform%Dcoefficients(2) = dbeta2
          call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(1,1))
        
        case (1)
          ! Assemble the convection directly including th´e stabilisation
          rconfigSD%dupsam = dupsam
          rconfigSD%dnu = dnu
          call elemdbg2d_sd (dbeta1,dbeta2,rconfigSD,rmatrix%RmatrixBlock(1,1))
          
        case (2)
          ! Assemble convective operator, sum up to the Laplace operator
          rform%itermcount = 2
          rform%Idescriptors(1,1) = DER_DERIV_X
          rform%Idescriptors(2,1) = DER_FUNC
          rform%Idescriptors(1,2) = DER_DERIV_Y
          rform%Idescriptors(2,2) = DER_FUNC
          rform%Dcoefficients(1) = dbeta1
          rform%Dcoefficients(2) = dbeta2
          call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(1,1))
          
          ! Assemble the jump stabilisation
          rconfigEOJ%dgamma = dgamma
          rconfigEOJ%dgammastar = dgamma
          rconfigEOJ%dnu = dnu
          call conv_JumpStabilisation2d (rconfigEOJ, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1))
          
          ! Subtract the boundary edges from the matrix.
          rconfigEOJ%dtheta = -1.0_DP

          allocate (p_Iedges(rtriangulation%NMT))
          allocate (p_Itemp(rtriangulation%NMT))

          ! Loop over the edges and boundary components
          do ibc = 1,boundary_igetNBoundComp(rboundary)
            do iedge = 1,boundary_igetNsegments(rboundary, ibc)

              ! Get a region identifying that boundary edge
              call boundary_createRegion(rboundary,ibc,iedge,rboundaryRegion)

              ! Get triangulation edges on that boundary edge
              call bcasm_getEdgesInBdRegion (rtriangulation,rboundaryRegion, &
                  nedgecount, p_Iedges, p_Itemp)

              ! Subtract the jump
              call conv_JumpStabilisation2d(rconfigEOJ,CONV_MODMATRIX,rmatrix%RmatrixBlock(1,1), &
                  InodeList=p_Iedges(1:nedgecount))

            end do
          end do

          deallocate (p_Iedges,p_Itemp)

          rconfigEOJ%dtheta = 1.0_DP
          
        case (3)
          ! Prepare a velocity field resembling our convection
          call spdiscr_initBlockDiscr (rdiscretisationVel,2,rtriangulation, rboundary)
          call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
              rdiscretisationVel%RspatialDiscr(1), .true.)
          call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
              rdiscretisationVel%RspatialDiscr(2), .true.)
          
          ! Create the velocity vector
          call lsysbl_createVectorBlock (rdiscretisationVel,rvelocityVector)
          call lsyssc_clearVector (rvelocityVector%RvectorBlock(1),dbeta1)
          call lsyssc_clearVector (rvelocityVector%RvectorBlock(2),dbeta2)
        
          ! Assemble the convection using the kernel stabilisation
          rconfigSD%dupsam = dupsam
          rconfigSD%dnu = dnu
          rconfigSD%clocalh = clocalh
          call conv_streamlineDiffusionBlk2d ( &
              rvelocityVector, rvelocityVector, 1.0_DP, 0.0_DP,&
              rconfigSD, CONV_MODMATRIX,rmatrix)
          
          ! Release all that temporary stuff.
          call lsysbl_releaseVector (rvelocityVector)
          call spdiscr_releaseBlockDiscr (rdiscretisationVel)
          
        end select
      
      end select
      
      ! Assemble RHS vector
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC2D
      call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
            rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS2D,rcollect)

      ! In any case except for L2-projection, filter the system
      if(itest .ne. 201) then
        ! Implement BCs
        call vecfil_discreteBCrhs (rvecRhs)
        call vecfil_discreteBCsol (rvecsol)
        call matfil_discreteBC (rmatrix)
      end if
      
      ! Set up the solver
      select case(isolver)
      case (0)
        ! UMFPACK solver
        call linsol_initUMFPACK4(p_rsolver)
      
      case (1)
        ! CG-SSOR[1.2] solver
        nullify(p_rprecond)
        call linsol_initSSOR(p_rprecond,drelax)
        call linsol_initCG(p_rsolver, p_rprecond)
      
      case (2)
        ! BiCGStab-ILU(k) solver with RCMK
        nullify(p_rprecond)
        call linsol_initMILUs1x1(p_rprecond,ifillin,0.0_DP)
        call linsol_initBiCGStab(p_rsolver, p_rprecond)
        
        ! Calculate a RCMK permutation
        call sstrat_calcRevCuthillMcKee(rmatrix%RmatrixBlock(1,1), h_Ipermute)
        
        ! Permute matrix and vectors
        call lsyssc_sortMatrix(rmatrix%RmatrixBlock(1,1), .true., &
                               SSTRAT_RCM, h_Ipermute)
        call lsyssc_sortVectorInSitu(rvecSol%RvectorBlock(1), &
            rvecTmp%RvectorBlock(1), SSTRAT_RCM, h_Ipermute)
        call lsyssc_sortVectorInSitu(rvecRhs%RvectorBlock(1), &
            rvecTmp%RvectorBlock(1), SSTRAT_RCM, h_Ipermute)
      
      end select
      
      p_rsolver%ioutputLevel = ioutput
      p_rsolver%depsRel = depsRel
      p_rsolver%depsAbs = depsAbs
      p_rsolver%nmaxiterations = nmaxiter

      Rmatrices = (/rmatrix/)
      call linsol_setMatrices(p_rsolver,Rmatrices)
      call linsol_initStructure (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      call linsol_initData (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Solve the system
      call linsol_solveAdaptively (p_rsolver,rvecSol,rvecRhs,rvecTmp)
      
      ! If we have a permutation, unsort the solution vector
      if(h_Ipermute .ne. ST_NOHANDLE) then
        call lsyssc_sortVectorInSitu(rvecSol%RvectorBlock(1), &
            rvecTmp%RvectorBlock(1), -SSTRAT_RCM)
      end if
      
      ! Calculate the errors to the reference function
      rerror%p_RvecCoeff => rvecSol%RvectorBlock(1:1)
      rerror%p_DerrorL2 => Derror(1:1,ilvl)
      rerror%p_DerrorH1 => Derror(2:2,ilvl)
      call pperr_scalarVec(rerror, getReferenceFunction2D, rcollect)
      
      ! Print the errors
      call output_line('Errors (L2/H1): ' // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)))

      ! Probably write the mesh to disc
      if (iwritemesh .eq. 1) then
        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                          'ucd/sol2d_'//trim(sys_siL(ilvl,5))//'.gmv')

        ! Project the solution to the vertices
        allocate (p_Ddata(rtriangulation%NVT))
        call spdp_projectToVertices (rvecSol%RvectorBlock(1), p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)

        call ucd_write (rexport)
        call ucd_release (rexport)
        deallocate(p_Ddata)
      else if (iwritemesh .eq. 2) then
        call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                          'ucd/sol2d_'//trim(sys_siL(ilvl,5))//'.vtk')

        ! Project the solution to the vertices
        allocate (p_Ddata(rtriangulation%NVT))
        call spdp_projectToVertices (rvecSol%RvectorBlock(1), p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)

        call ucd_write (rexport)
        call ucd_release (rexport)
        deallocate(p_Ddata)
      end if

      ! Clean up this level
      call linsol_doneData (p_rsolver)
      call linsol_doneStructure (p_rsolver)
      call linsol_releaseSolver (p_rsolver)
      call lsysbl_releaseVector (rvecTmp)
      call lsysbl_releaseVector (rvecSol)
      call lsysbl_releaseVector (rvecRhs)
      call lsysbl_releaseMatrix (rmatrix)
      if(h_Ipermute .ne. ST_NOHANDLE) call storage_free(h_Ipermute)
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      call spdiscr_releaseBlockDiscr(rdiscretisation)
      call tria_done (rtriangulation)
    
    end do ! ilvl
    
    ! Release the domain
    call boundary_release (rboundary)
    
    ! Print some statistics
    call output_separator(OU_SEP_MINUS)
    call output_line('Level         NEQ        NNZE         NVT' // &
                     '         NMT         NEL')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // &
          trim(sys_si(Istat(1,ilvl),12)) // &
          trim(sys_si(Istat(2,ilvl),12)) // &
          trim(sys_si(Istat(3,ilvl),12)) // &
          trim(sys_si(Istat(4,ilvl),12)) // &
          trim(sys_si(Istat(5,ilvl),12)))
    end do ! ilvl

    ! Print out the L2- and H1-errors on each level
    call output_separator(OU_SEP_MINUS)
    call output_line('Level     L2-error            H1-error')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // '   ' // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)))
    end do ! ilvl
    
    ! Print out the L2- and H1- factors for each level pair
    if(NLMAX .gt. NLMIN) then
      call output_separator(OU_SEP_MINUS)
      call output_line('Level     L2-factor           H1-factor')
      do ilvl = NLMIN+1, NLMAX
        
        ! avoid division by zero here
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        if(abs(Derror(1,ilvl)) .gt. SYS_EPSREAL_DP) &
          daux1 = Derror(1,ilvl-1)/Derror(1,ilvl)
        if(abs(Derror(2,ilvl)) .gt. SYS_EPSREAL_DP) &
          daux2 = Derror(2,ilvl-1)/Derror(2,ilvl)
        
        ! print out the factors
        call output_line(trim(sys_si(ilvl,5)) // '   ' // &
            trim(sys_sdEP(daux1,20,12)) // &
            trim(sys_sdEP(daux2,20,12)))
      end do ! ilvl
    end if

    ! Deallocate arrays
    deallocate(Istat)
    deallocate(Derror)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine elemdbg2d_sd (dbeta1,dbeta2,rconfig,rmatrix)

!<description>
  ! Standard streamline diffusion for the convective operator (beta1,beta2).
  ! Sets up the matrix rmatrix.
!</description>

!<input>

  ! Constant convection field.
  real(DP) :: dbeta1,dbeta2
  
  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamlineDiffusion), intent(IN) :: rconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: celement
    
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      print *,'SD: Unsupported matrix format'
      call sys_halt()
    end if

    celement = rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement
    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'SD: Unsupported discretisation.'
      call sys_halt()
    end if

    if (.not. rconfig%bconstViscosity) then
      print *,'SD: Only constant viscosity supported at the moment!'
      call sys_halt()
    end if
    
    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      print *,'SD: Viscosity parameter nu not initialised!'
      call sys_halt()
    end if
    
    ! Calculate only the matrix. "The" matrix is used by the caller
    ! on all diagonal blocks of a block matrix!
    call elemdbg2d_strdiff2d_double ( &
                  dbeta1,dbeta2,rmatrix, rconfig%dupsam, rconfig%dnu, &
                  rconfig%clocalh)

    !!! DEBUG:
    !call matio_writeMatrixHR (rmatrix, 'matrix',&
    !                          .TRUE., 0, 'matrixL.txt', '(D10.3)')

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine elemdbg2d_strdiff2d_double (dbeta1,dbeta2, &
                  rmatrix,dupsam,dnu,clocalh)
!<description>
  ! Standard streamline diffusion method to set up the operator beta*grad(u).

!</description>

!<input>
  ! Convection direction.
  real(DP), intent(in) :: dbeta1,dbeta2
  
  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(IN) :: dupsam
  
  ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
  real(DP), intent(IN) :: dnu
  
  ! Method how to compute the local h
  integer, intent(in) :: clocalh
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  type(t_matrixScalar), intent(INOUT), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IDOFE,JDOFE,icubp
  integer :: JCOL0,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,OM,AH,dre
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
  integer :: NVE
  
  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation
  
  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet

  ! Arrays for saving Jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj
  
  ! An allocateable array accepting the DOF's of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs
  
  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity
  
  ! An array with local DELTA's, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
  
    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! For ALE we don't even need so much
    BderALE = .false.
    BderALE(DER_FUNC) = .true.
    
    ! Shortcut to the spatial discretisation
    p_rdiscretisation => rmatrix%p_rspatialDiscrTest
    
    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)
    
    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    
    ! Get the number of local DOF's for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF's Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(p_relementDistribution%celement)
    
    ! Number of local DOF's
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(BILF_NELEMSIM,p_rtriangulation%NEL)
    
    ! Get matrix arrays
    call lsyssc_getbase_double (rmatrix,p_Da)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeBilForm)
    
    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(p_relementDistribution%ccubTypeBilForm,p_DcubPtsRef, Domega)
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA's
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    allocate(Kentry(indof,indof,nelementsPerBlock))
    allocate(Dentry(indof,indof))
    
    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))
    
    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu
    else
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if

    call lalg_clearVectorDble (DlocalDelta)
    
    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it's too small...

    du1loc = dbeta1
    du2loc = dbeta2
    dunorm = sqrt(du1loc**2+du2loc**2)
    dumax = dunorm

    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)
    
    ! Loop over the elements - blockwise.
    do IELset = 1, size(p_IelementList), BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(p_IelementList),IELset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)
                                  
      ! Calculate local DELTA's for streamline diffusion method.
      ! (cf. p. 121 in Turek's CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we don't do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      call getLocalDeltaQuad (clocalh,dbeta1,dbeta2,&
                    p_IelementList(IELset:IELmax),&
                    duMaxR,DlocalDelta,p_rtriangulation,Idofs,dupsam,dre)
                                   
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j)
      ! "active" (i.e. have common support) on our current element, each
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position
      !   in the global system matrix, where the corresponding value
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately,
      !  but we directly add them to the global matrix in this
      !  approach.)
      !
      ! We build local matrices for all our elements
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof
        
          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do JDOFE=1,indof
            
            ! Get the global DOF of the "X" which interacts with
            ! our "O".
            
            JDFG=Idofs(JDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.
            
            do JCOL=JCOL0,rmatrix%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

            ! Because columns in the global matrix are sorted
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          end do ! IDOFE
          
        end do ! JDOFE
        
      end do ! IEL
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      if (IELset .eq. 1) then
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

      ! We want to set up the nonlinear part of the matrix
      !
      !   n~_h (u_h, u_h, v_h)
      !
      ! = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
      !   ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !  standard nonlin. part                  stabilization
      !
      ! More precisely, as we want to assemble the matrix which is
      ! later multiplied with coefficient vectors, we have to insert
      ! basis functions in the above terms instead of u_h and v_h.
      ! Assuming the representation u_h=sum_j(u_j*Phi_j) and
      ! v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
      ! DOF's as:
      !
      !   n_h (u_h, Phi_j, Phi_i)
      ! + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
      !
      ! In nonstationary simulations, the system matrix typically
      ! contains a mass matrix to respect the time derivative.
      ! The matrix has the form
      !
      ! [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
      !
      ! In a first step, we calculate the velocity field in all
      ! cubature points on all elements of the current block.
      ! If we only have a primary velocity field
      ! (dweight2=0), we can calculate that only by summing up the
      ! velocities in U1Lx, otherwise we have to sum up
      ! dweight1*u1vel + dweight2*u2vel
      
      ! Loop over all elements in the current set
      do IEL=1,IELmax-IELset+1
      
        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp
        
          ! Save the computed velocity
          Dvelocity(1,ICUBP,IEL) = dbeta1
          Dvelocity(2,ICUBP,IEL) = dbeta2
        
        end do ! ICUBP
        
      end do ! IEL
          
      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a biliinear form!
      !
      ! Loop over the elements in the current set.

      do IEL=1,IELmax-IELset+1
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        
        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

          ! Current velocity in this cubature point:
          du1loc = Dvelocity (1,ICUBP,IEL)
          du2loc = Dvelocity (2,ICUBP,IEL)
          
          ! We take a more detailed look onto the last scalar product
          ! of n~_h (u_h, u_h, v_h) what we want to calculate here.
          !
          ! The vector u_h=(DU1,DU2) contains both velocity components,
          ! for the X as well as for the Y velocity. On the other hand
          ! the system matrix we want to build here will be designed for
          ! one velocity component only! Therefore, Phi_i and Phi_j
          ! are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
          ! with two components. Therefore, the last scalar product is more
          ! in detail:
          !
          !     ( u_h*grad Phi_j, u_h*grad Phi_i )_T
          !
          ! =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
          !         (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)
          !
          ! =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
          !       (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
          !
          ! =   HSUMJ * HSUMI
          !
          ! i.e. a product of two scalar values!
          !
          ! Summing up over all pairs of multiindices.
          !
          ! Outer loop over the DOF's i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do IDOFE=1,indof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:
          
            HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
            HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
            HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)
           
            ! Calculate
            !
            !     U * grad(Phi_i)  =  < grad(Phi_i), U >
            !
            !   = ( grad(Phi_i)_1 , (DU1) )
            !     ( grad(Phi_i)_2   (DU2) )
            !
            ! Remember: DU1MV=DU2MV=0 in this case.
            !
            ! If ALE is active, use v=mesh velocity and calculate
            !
            !     (U-v) * grad(Phi_i)  =  < grad(Phi_i), U-v >
            !
            !   = ( grad(Phi_i)_1 , (DU1-DU1MV) )
            !     ( grad(Phi_i)_2   (DU2-DU2MV) )

            HSUMI = HBASI2*du1loc + HBASI3*du2loc

            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do JDOFE=1,indof
              
              if (IDOFE.eq.JDOFE) then
              
                ! Short version of the evaluation of the matrix
                ! contribution - see below for a more detailed
                ! description what is added together here!
              
                AH = HSUMI*(DlocalDelta(IEL)*HSUMI+HBASI1)
    
              else
              
                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:
              
                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)

                ! Calculate
                !
                !     U * grad(Phi_j)  =  < grad(Phi_j), U >
                !
                !   = ( grad(Phi_j)_1 , (DU1) )
                !     ( grad(Phi_j)_2   (DU2) )
                !
                ! Remember: DU1MV=DU2MV=0 in this case.
                !
                ! If ALE is active, use v=mesh velocity and calculate
                !
                !     (U-v) * grad(Phi_j)  =  < grad(Phi_j), U-v >
                !
                !   = ( grad(Phi_j)_1 , (DU1-DU1MV) )
                !     ( grad(Phi_j)_2   (DU2-DU2MV) )
                !
                ! But as v is already incorporated into DVelocity,
                ! we don't have to worry about that.

                HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc
    
                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of DNU,
                ! dalpha,ddelta,... this decomposes into three
                ! different parts:
                !
                ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
                !    + dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                !    + dalpha*(phi_j*phi_i)         | Mass matrix
                !
                ! The last two parts are probably not added to the
                ! matrix by setting DNY or CT0 to 0, respectively.
                !
                ! For saving some numerical operations, we write:
                !
                !     HSUMJ * (Delta * HSUMI + HBASI1)
                !
                ! =   Delta * HSUMJ * HSUMI
                !   + HSUMJ * HBASI1
                !
                ! =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
                !   + (U*grad(Phi_j),Phi_i)
                !
                ! <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
                !     + n_h (u_h, Phi_j, Phi_i)
                !
                ! plus the terms for the Stokes and Mass matrix,
                ! if their coefficient is <> 0.
                
                AH = HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1)
    
              end if ! (IDOFE.EQ.JDOFE)

              ! Weighten the calculated value AH by the cubature
              ! weight OM and add it to the local matrix. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+OM*AH
              
            end do ! IDOFE
            
          end do ! JDOFE

        end do ! ICUBP
        
        ! Now we have set up a "local" system matrix. We can either
        ! include it into the real matrix or we can use it to simply
        ! modify the RHS vector to create a defect vector (throwing
        ! away the information about the matrix afterwards, which would
        ! result in a matrix free modification of the RHS vector).
        !
        ! For cdef= containing CONV_MODMATRIX, incorporate our "local" system matrix
        ! into the global matrix. The position of each entry DENTRY(X,Y)
        ! in the global matrix array A was saved in element Kentry(X,Y)
        ! before.
        ! Kentry gives the position of the additive contributions in Dentry.
        ! The entry is weighted by the current dtheta, which is usually
        ! the weighting parameter of the corresponding THETA-scheme of a
        ! nonstationary simulation. For stationary simulations, dtheta is typically
        ! 1.0 which includes the local matrix into the global one directly.)
        
        do IDOFE=1,indof
          do JDOFE=1,indof
            p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + Dentry(JDOFE,IDOFE)
          end do
        end do
        
      end do ! IEL

    end do ! IELset
    
    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(p_DcubPtsRef)
    deallocate(Domega)
    deallocate(DlocalDelta)
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(Idofs)
    deallocate(Dbas)

  end subroutine

  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaQuad (clocalh,&
                      dbeta1,dbeta2,Ielements,&
                      duMaxR,Ddelta,rtriangulation,Idofs,dupsam,dnurec)

  ! This routine calculates a local ddelta=DELTA_T for a set of finite
  ! elements Ielements. This can be used by the streamline diffusion
  ! stabilisation technique as a multiplier of the (local) bilinear form.
  !
  ! The effective velocity that is used for calculating the ddelta
  ! is combined by a weighted mean of the two velocity fields U1,U2
  ! by:
  !                   du = (beta1,beta2)
  ! The coefficients A1,A2 allow the caller to take influence on which
  ! velocity field to weight more.
  
  ! Method how to compute the local h.
  ! =0: Use the root of the area of the element as local H
  ! =1: Use the length of the way that a particle travels through
  !     the element in direction of the flow
  integer, intent(in) :: clocalH
  
  ! Velocity direction
  real(DP) :: dbeta1,dbeta2
  
  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(IN) :: duMaxR
  
  ! Reciprocal value 1/NU of coefficient NU in front of the
  ! Laplacian term of the Navier-Stokes equation
  !   NU * Laplace(u) + u*grad(u) + ...
  real(DP), intent(IN) :: dnuRec
  
  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(IN) :: dupsam
  
  ! List of elements where the Ddelta should be calculated
  integer, dimension(:), intent(IN) :: Ielements
  
  ! Array with global degrees of freedom on the elements
  integer, dimension(:,:), intent(IN) :: Idofs
  
  ! Triangulation that defines the mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Out: local Ddelta on all elements
  real(DP), dimension(:), intent(OUT) :: ddelta

  ! local variables
  real(DP) :: dlocalH,du1,du2,dunorm,dreLoc
  integer :: iel,ielidx
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(DP), dimension(:), pointer :: p_DelementVolume

    ! Get some crucial data
    if (clocalh .eq. 0) then
      call storage_getbase_double (rtriangulation%h_DelementVolume,p_DelementVolume)
      
      ! Loop through all elements
      do ielidx = 1,size(Ielements)
      
        iel = Ielements(ielidx)

        ! Loop through the local degrees of freedom on element IEL.
        ! Sum up the velocities on these DOF's. This will result
        ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
        ! through element IEL.

        ! For elements whose DOF's represent directly the velocity, U1/U2
        ! represent the mean velocity
        ! along an egde/on the midpoint of each edge, so U1/U2 is
        ! clearly an approximation to the velocity in element T.

        du1=dbeta1
        du2=dbeta2

        ! Calculate the norm of that local velocity:

        dunorm = sqrt(du1**2+du2**2) / real(ubound(Idofs,1),DP)
        
        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then
        
          Ddelta(ielidx) = 0.0_DP

        else

          ! Calculate the local h from the area of the element
          dlocalH = sqrt(p_DelementVolume(iel))

          ! Calculate ddelta... (cf. p. 121 in Turek's CFD book)

          if (dupsam .lt. 0.0_DP) then

            ! For UPSAM<0, we use simple calculation of ddelta:
          
            Ddelta(ielidx) = abs(dupsam)*dlocalH
            
          else
          
            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU
            
            dreLoc = dunorm*dlocalH*dnuRec
            
            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
            
            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))
            
          end if ! (UPSAM.LT.0.0)
          
        end if ! (dunorm.LE.1D-8)

      end do
      
    else
    
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

      ! Loop through all elements
      do ielidx = 1,size(Ielements)
      
        iel = Ielements(ielidx)

        ! Loop through the local degrees of freedom on element IEL.
        ! Sum up the velocities on these DOF's. This will result
        ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
        ! through element IEL.

        ! For elements whose DOF's represent directly the velocity, U1/U2
        ! represent the mean velocity
        ! along an egde/on the midpoint of each edge, so U1/U2 is
        ! clearly an approximation to the velocity in element T.

        du1=dbeta1
        du2=dbeta2

        ! Calculate the norm of that local velocity:

        dunorm = sqrt(du1**2+du2**2) / real(ubound(Idofs,1),DP)
        
        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then
        
          Ddelta(ielidx) = 0.0_DP

        else

          ! u_T defines the "slope" of the velocity through
          ! the element T. At next, calculate the local mesh width
          ! dlocalH = h = h_T on our element T=IEL:

          call getLocalMeshWidthQuad (dlocalH,dunorm, du1, du2, iel, &
              p_IverticesAtElement,p_DvertexCoords)

          ! Calculate ddelta... (cf. p. 121 in Turek's CFD book)

          if (dupsam .lt. 0.0_DP) then

            ! For UPSAM<0, we use simple calculation of ddelta:
          
            Ddelta(ielidx) = abs(dupsam)*dlocalH
            
          else
          
            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU
            
            dreLoc = dunorm*dlocalH*dnuRec
            
            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
            
            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))
            
          end if ! (UPSAM.LT.0.0)
          
        end if ! (dunorm.LE.1D-8)

      end do

    end if

  end subroutine

  ! ----------------------------------------------------------------------

  pure subroutine getLocalMeshWidthQuad (dlocalH, dunorm,  XBETA1, &
                      XBETA2, JEL,Kvert,Dcorvg)
  
  ! Determine the local mesh width for an element JEL of a
  ! triangulation.
  
  ! Element where the local h should be calculated
  integer, intent(IN)               :: JEL
  
  integer, dimension(TRIA_MAXNVE2D,*), intent(IN) :: Kvert
  real(DP), dimension(NDIM2D,*), intent(IN)          :: Dcorvg
  
  ! norm ||u||_T = mean velocity through element T=JEL
  real(DP), intent(IN)  :: dunorm
  
  ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
  real(DP), intent(IN)  :: XBETA1, XBETA2
  
  ! local mesh width
  real(DP), intent(OUT) :: dlocalH
  
  ! local variables
  real(DP) :: dlambda
  integer :: NECK1,NECK2,NECK3,NECK4
  real(DP) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
  real(DP) :: dalphaMax, dalpha

    ! Fetch the numbers of the four corners of element JEL

    neck1=Kvert(1,JEL)
    neck2=Kvert(2,JEL)
    neck3=Kvert(3,JEL)
    neck4=Kvert(4,JEL)

    ! Fetch the coordinates of these corners

    x1=Dcorvg(1,neck1)
    y1=Dcorvg(2,neck1)
    x2=Dcorvg(1,neck2)
    y2=Dcorvg(2,neck2)
    x3=Dcorvg(1,neck3)
    y3=Dcorvg(2,neck3)
    x4=Dcorvg(1,neck4)
    y4=Dcorvg(2,neck4)

    ! Scale: (deactivated)

    !  dsp=max(xbeta1,xbeta2)

    !  xbeta1=xbeta1
    !  xbeta2=xbeta2

    dalphaMax=0.0_DP
    
    ! In the next step, we calculate the 'maximum possible mesh with
    ! in direction of the flow'; this is the maximum possible length
    ! that a particle can cross in the current element.
    ! The picture in mind is the following:
    !
    !          G3
    !   +-------------X-------+
    !   |            /        |
    !   |           /         |
    !   |          /          |
    !   |         /           |
    !   |        /            |
    ! G4|       /             | G2
    !   |      ^ (beta1,beta2)|
    !   |     /               |
    !   |    /                |
    !   |   /                 |
    !   |  /                  |
    !   | /                   |
    !   |/                    |
    !   O---------------------+
    !            G1
    !
    ! The vector (beta1,beta2) gives the direction of the flow.
    ! A particle starting in point O and moves at most up to point X.
    ! The length of the line (O,X) is the local mesh with h.
    !
    ! Loop through the four corners of element JEL and check
    ! of a line with slope BETA=(xbeta1,xbeta2) starting in this
    ! corner really intersects with one of the edges of the element.
    ! Remark that we only have to check the two opposite edges
    ! to the current corner!

    ! -----------------------------------------------------------------
    ! Check the first corner:

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The second one...
    
    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The third one...
    
    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! And the fourth=last one...
    
    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)

    ! -----------------------------------------------------------------
    ! finally determine the local h=h_T
    !
    ! dalphaMax is the maximum alpha, normalised as 'parameter value',
    ! i.e. dalphaMax=1.0 corresponds to a vector 1.0*(dbeta1,dbeta2).
    ! We multiply with dunorm=|(dbeta1,dbeta2)| to get the actual length
    ! of the vector which can be placed inside of the element.
    !
    ! Furthermore, we multiply with an additional weight 4. (why ?!?)

    dlocalH=dalphaMax*4.0_DP*dunorm

  end subroutine
  
  ! ----------------------------------------------------------------------

  pure subroutine intersectLines2D (XO,YO,dalpha,BETA1,BETA2, &
                      XA,YA,dlambda,XB,YB)

  ! Intersect two lines in R^2

  ! Origin of line 1
  real(DP), intent(IN) :: XO,YO
  
  ! Direction of line 1
  real(DP), intent(IN) :: BETA1,BETA2
  
  ! One point on the second line
  real(DP), intent(IN) :: XA,YA
  
  ! Another point on the second line
  real(DP), intent(IN) :: XB,YB
  
  ! Parameter value of the intersection point on line 1.
  ! =0.0, if there is no intersection point
  real(DP), intent(OUT) :: dalpha
  
  real(DP), intent(OUT) :: dlambda
  
  ! local variables
  double precision :: dsp

    ! Scalar product of the line (xa,ya)->(xb,yb) with the
    ! counterclockwise normal n1 of (beta1,beta2)
    dsp=BETA2*(XB-XA)-BETA1*(YB-YA)
    
    if (dsp.eq.0.0_DP) then
    
      ! beta and the vector are parallel
      dalpha=0.0_DP
      
    else

      ! Scalar product of (beta1,beta2) with the (inner) normal vector n2
      ! of the line (xo,yo)->(xa,ya).
      dlambda=(BETA1*(YA-YO)-BETA2*(XA-XO))/dsp

      !                    (xb,yb)
      !   +-----------------+
      !   |                 |
      !   |                 |
      !   ^ n2              |
      !   !                 |
      !   !  (beta1,beta2)  |    (beta1,beta2)
      !   !    ^            |    ^
      !   !   /  ^__ n1     |   /
      !   !  /      \__     |  /
      !   ! /          \__  | /
      !   !/              \_|/
      !   +-----------------+
      ! (xo,yo)            (xa,ya)
      !
      ! (What is this? Documentation incomplete. Has someone a good
      ! reference?)

      ! is the intersection point inside of the element?
      if ((dlambda.ge.-1E-1_DP).and.(dlambda.le.1.11E0_DP)) then
        if (BETA1 .ne. 0.0_DP) then
          dalpha=((XA-XO)+dlambda*(XB-XA))/BETA1
        else
          if (BETA2 .ne. 0.0_DP) then
            dalpha=((YA-YO)+dlambda*(YB-YA))/BETA2
          else
            dalpha=0.0_DP
          end if
        end if
      else
        dalpha=0.0_DP
      end if
      
    end if

  end subroutine

end module
